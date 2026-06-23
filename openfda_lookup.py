# -*- coding: utf-8 -*-
"""
openFDA Drug Event API (FAERS) から、薬剤ペアの併用報告における
死亡・重篤転帰の件数と割合を取得するモジュール。

データソース: https://api.fda.gov/drug/event.json
  - 登録・APIキー不要、無料、REST/JSON（project_interaction_checker.md で検証済み）
  - 米国データ(FAERS)であり日本国内の報告傾向とは異なる点に留意

注意: FAERSは医薬品名を英語(INN/商品名)で記録しているため、
日本語の一般名のままでは検索できない。英語名への変換は
drug_name_map.json の対応表を使う（未登録の薬剤は lookup_pair が
None を返すので、呼び出し側が英語名を直接 co_report_stats に
渡すか、対応表に追記して使う）。
"""
import concurrent.futures
import json
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import severity

BASE = "https://api.fda.gov/drug/event.json"
NAME_MAP_FILE = Path(__file__).parent / "drug_name_map.json"
CACHE_DIR = Path(__file__).parent / "cache" / "openfda"
# FAERS全体の報告総数。ROR(報告オッズ比)の2×2表で「いずれの薬剤も含まない報告数」を
# 求める分母に使う。実数はAPIから取得してキャッシュするが、取得失敗時のフォールバック値
# として直近の概数(約2,033万件/2026-06)を置く。総数は四半期ごとに緩やかに増えるだけで、
# 2,000万規模に対する数%の差はRORの桁を変えないため、フォールバックでも実害は小さい。
_DB_TOTAL_FALLBACK = 20_328_575
CACHE_DIR.mkdir(parents=True, exist_ok=True)


def _load_name_map() -> dict:
    if NAME_MAP_FILE.exists():
        return json.loads(NAME_MAP_FILE.read_text(encoding="utf-8"))
    return {}


def to_english(jp_name: str):
    """日本語の一般名(例: アムロジピン)を openFDA検索用の英語名(例: amlodipine)に変換。
    対応表に無ければ None を返す。"""
    return _load_name_map().get(jp_name)


def _escape(name: str) -> str:
    """openFDAクエリのフレーズ検索に埋め込めるよう薬剤名をエスケープ＆URLエンコードする。

    "cefcapene pivoxil"のようにスペースを含む英語名をそのままURLに渡すと
    urlopenが`InvalidURL: URL can't contain control characters`で拒否するため、
    quote()でパーセントエンコードする（openFDA側でデコードされ、フレーズとして
    一致するため検索結果には影響しない）。
    """
    return urllib.parse.quote(name.replace('"', ""), safe="")


def _request(query_string: str):
    url = f"{BASE}?{query_string}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            return json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # openFDAは0件ヒット時に404を返す
            return {"meta": {"results": {"total": 0}}, "results": []}
        raise


def _total(search_query: str) -> int:
    data = _request(f"search={search_query}&limit=1")
    return data.get("meta", {}).get("results", {}).get("total", 0) if data else 0


def _top_reactions(search_query: str, n: int = 5):
    data = _request(f"search={search_query}&count=patient.reaction.reactionmeddrapt.exact&limit={n}")
    if not data:
        return []
    return [(d["term"], d["count"]) for d in data.get("results", [])]


def _db_total() -> int:
    """FAERS全体の報告総数を取得する（RORの分母用）。日次で変わるほどの値ではないため
    cache/openfda/_db_total.json に丸ごとキャッシュし、取得失敗時はフォールバック概数を返す。"""
    cache_file = CACHE_DIR / "_db_total.json"
    if cache_file.exists():
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))["total"]
        except Exception:
            pass
    try:
        data = _request("limit=1")
        total = data.get("meta", {}).get("results", {}).get("total", 0)
    except Exception:
        total = 0
    total = total or _DB_TOTAL_FALLBACK
    cache_file.write_text(json.dumps({"total": total}, ensure_ascii=False), encoding="utf-8")
    return total


def _ror(co: int, solo_a: int, solo_b: int, n_total: int) -> dict:
    """併用報告の不均衡を報告オッズ比(ROR)で評価する。

    2×2表（各セル＝報告件数）:
        a = 両剤を含む報告           = co
        b = A を含むが B を含まない  = solo_a - co
        c = B を含むが A を含まない  = solo_b - co
        d = どちらも含まない         = n_total - solo_a - solo_b + co
      ROR = (a·d)/(b·c)、95%信頼区間は ln(ROR) ± 1.96·√(1/a+1/b+1/c+1/d)。

    生の死亡件数ではなく ROR を使う理由: 頻用薬は相互作用が無くても併用報告が
    大量に出る（＝絶対数の閾値はよく使われる薬ほど無条件に超える）。ROR は各剤
    単独の報告頻度から期待される併用数を基準に「偶然より何倍多く併用報告されたか」
    を測るため、薬剤の使用頻度の影響を補正できる。

    注意: FAERS上の「併用報告が多い」は相互作用だけでなく“併用処方が多い”ことでも
    起こり得る（同じ疾患に併用される2剤など）。RORの高さ＝相互作用の強さと
    短絡せず、あくまで「一緒に報告されやすいか」のシグナルとして扱う。

    いずれかのセルが0の場合は Haldane-Anscombe 補正(全セルに0.5加算)を施す。
    戻り値: {"ror","ci_low","ci_high","expected","obs_exp"} 算出不能時は None。
    """
    a, b, c = co, solo_a - co, solo_b - co
    d = n_total - solo_a - solo_b + co
    if a <= 0 or b < 0 or c < 0 or d <= 0:
        return None
    if min(a, b, c, d) == 0:  # ゼロセルがあると対数・除算が発散するため補正
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    ci_low = math.exp(math.log(ror) - 1.96 * se)
    ci_high = math.exp(math.log(ror) + 1.96 * se)
    expected = solo_a * solo_b / n_total if n_total else 0  # 独立を仮定した期待併用報告数
    return {
        "ror": round(ror, 2),
        "ci_low": round(ci_low, 2),
        "ci_high": round(ci_high, 2),
        "expected": round(expected, 1),
        "obs_exp": round(co / expected, 2) if expected else None,
    }


def co_report_stats(drug_a_en: str, drug_b_en: str, use_cache: bool = True, polite: bool = True) -> dict:
    """
    2剤の英語名(openFDA検索名)から、併用報告の件数・死亡/重篤割合を集計して返す。

    戻り値:
      {
        "drug_a", "drug_b": 検索に使った英語名,
        "co_reports_total": 併用報告の総数,
        "co_reports_death": うち死亡転帰(reactionoutcome=5)の件数,
        "co_reports_serious": うち重篤(serious=1)の件数,
        "death_ratio": 死亡件数 / 総数,
        "serious_ratio": 重篤件数 / 総数,
        "solo_a_total" / "solo_b_total": 各剤単独の全報告数(比較用ベースライン),
        "top_reactions": [(有害事象名, 件数), ...],  # 併用報告での頻出事象Top5
        "ror": 報告オッズ比, "ror_ci_low"/"ror_ci_high": 95%信頼区間,
        "expected_co": 独立を仮定した期待併用報告数, "obs_exp": 実測/期待比
        （RORが算出不能な場合これらは None）
      }
    """
    safe_a, safe_b = _escape(drug_a_en), _escape(drug_b_en)
    cache_file = CACHE_DIR / f"{safe_a}__{safe_b}.json"
    if use_cache and cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    base_q = f'patient.drug.medicinalproduct:"{safe_a}"+AND+patient.drug.medicinalproduct:"{safe_b}"'

    # 6回のAPI呼び出しを逐次に投げると約7秒かかる。互いに独立なので並列化する。
    # ラウンド1: 併用報告総数・各剤単独報告数（total>0かに依存しない3本）を同時取得。
    # ラウンド2: 死亡/重篤/頻出事象（いずれもtotal>0のときだけ必要）を同時取得。
    # この2ラウンド化で実効ラウンドトリップが6→2に減る。
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        f_total = ex.submit(_total, base_q)
        f_solo_a = ex.submit(_total, f'patient.drug.medicinalproduct:"{safe_a}"')
        f_solo_b = ex.submit(_total, f'patient.drug.medicinalproduct:"{safe_b}"')
        total = f_total.result()
        solo_a = f_solo_a.result()
        solo_b = f_solo_b.result()

    if total:
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
            f_death = ex.submit(_total, f"{base_q}+AND+patient.reaction.reactionoutcome:5")
            f_serious = ex.submit(_total, f"{base_q}+AND+serious:1")
            f_top = ex.submit(_top_reactions, base_q)
            death = f_death.result()
            serious = f_serious.result()
            top = f_top.result()
        ror = _ror(total, solo_a, solo_b, _db_total())
    else:
        death = serious = 0
        top = []
        ror = None

    result = {
        "drug_a": drug_a_en, "drug_b": drug_b_en,
        "co_reports_total": total,
        "co_reports_death": death,
        "co_reports_serious": serious,
        "death_ratio": round(death / total, 4) if total else 0.0,
        "serious_ratio": round(serious / total, 4) if total else 0.0,
        "solo_a_total": solo_a,
        "solo_b_total": solo_b,
        "top_reactions": top,
        "ror": ror["ror"] if ror else None,
        "ror_ci_low": ror["ci_low"] if ror else None,
        "ror_ci_high": ror["ci_high"] if ror else None,
        "expected_co": ror["expected"] if ror else None,
        "obs_exp": ror["obs_exp"] if ror else None,
    }
    cache_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    if polite:
        time.sleep(0.5)  # 連続バッチ実行時のopenFDAサーバ負荷軽減（対話リクエストではpolite=Falseで省く）
    return result


# 添付文書から推定した英語名（例 "Rosuvastatin Calcium"）は塩・水和物が付くことがある。
# FAERSは中核INN名（rosuvastatin）での報告が多いので、末尾の対イオン塩・水和物語を
# 落とした候補も作る。先頭にある語（"sodium bicarbonate"のsodium等）は名前の一部なので
# 落とさない＝末尾一致のみ剥がす。
_EN_SALT_RE = re.compile(
    r"\s+(?:calcium|sodium|potassium|magnesium|zinc|"
    r"hydrochloride|hydrobromide|hydroiodide|sulfate|sulphate|phosphate|nitrate|"
    r"maleate|fumarate|succinate|tartrate|bitartrate|citrate|besilate|besylate|"
    r"mesilate|mesylate|tosilate|tosylate|acetate|benzoate|pamoate|embonate|"
    r"gluconate|lactate|aspartate|edisylate|napadisilate|"
    r"hydrate|hemihydrate|monohydrate|dihydrate|trihydrate|anhydrous)+$",
    re.IGNORECASE,
)


def _guess_candidates(raw: str):
    """添付文書由来の英語名(raw)から、openFDA照合用の候補を優先順に返す。

    末尾の塩・水和物を剥がした中核名(より広くヒット)を先に、剥がす前の正式名を後に置く。
    両方とも小文字化し、重複は除く。"""
    if not raw:
        return []
    full = re.sub(r"\s+", " ", raw).strip().lower()
    core = _EN_SALT_RE.sub("", full).strip()
    out = []
    for c in (core, full):
        if c and c not in out:
            out.append(c)
    return out


def _validate_name(en: str, use_cache: bool = True, polite: bool = True) -> bool:
    """英語名enがopenFDAに実報告として存在するか（単独報告数>0か）を確認する。

    添付文書からの自動推定名は綴り違い・誤抽出・剥がしすぎがあり得るため、実在しない
    名前で“それらしい嘘の統計”を出さないよう、採用前にこの実在確認を必ず通す。
    名前単位でcache/openfda/_val_*.jsonにキャッシュする。"""
    safe = _escape(en)
    cache_file = CACHE_DIR / f"_val_{safe}.json"
    if use_cache and cache_file.exists():
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))["ok"]
        except Exception:
            pass
    cnt = _total(f'patient.drug.medicinalproduct:"{safe}"')
    ok = cnt > 0
    cache_file.write_text(json.dumps({"name": en, "count": cnt, "ok": ok}, ensure_ascii=False),
                          encoding="utf-8")
    if polite:
        time.sleep(0.3)
    return ok


def _resolve_english(query: str, matched_name: str, guess: str = None,
                     use_cache: bool = True, polite: bool = True):
    """検索クエリ→添付文書上の正式名→正式名の中核成分名(severity._name_core)の順に
    drug_name_map.json と照合する。

    drug_name_map.json は中核成分名（例: 「アムロジピン」）で登録する運用のため、
    クエリが先発品名（例: 「ノルバスク」）や塩・水和物付きの正式名（例:
    「アムロジピンベシル酸塩」）であっても、正式名を中核名へ正規化することで
    対応表のエントリに到達できるようにする。

    手動対応表で見つからない場合のフォールバックとして、添付文書から推定した英語名
    (guess)を openFDA で実在確認(_validate_name)してから採用する。これにより手動
    メンテに頼らず大半の薬剤をカバーしつつ、実在しない推定名で誤った統計を出すのを防ぐ
    （確認に通らなければ従来通り None＝fda_name_missing にフォールバック）。
    """
    for candidate in (query, matched_name, severity._name_core(matched_name)):
        en = to_english(candidate) if candidate else None
        if en:
            return en  # 検証済みの手動対応表を最優先
    for cand in _guess_candidates(guess):
        if _validate_name(cand, use_cache=use_cache, polite=polite):
            return cand
    return None


def lookup_pair(query_a: str, matched_a: str, query_b: str, matched_b: str,
                guess_a: str = None, guess_b: str = None,
                use_cache: bool = True, polite: bool = True):
    """日本語の薬剤名2つ(検索クエリと添付文書上の正式名)から
    openFDA併用報告統計を取得する高レベル関数。

    guess_a/guess_b には添付文書から推定した英語名(pmda_lookupのenglish_name_guess)を
    渡せる。手動対応表に無い薬でも、推定名をopenFDAで実在確認できれば統計を取得する。
    いずれの方法でも英語名が確定できない場合は None を返す。"""
    en_a = _resolve_english(query_a, matched_a, guess_a, use_cache=use_cache, polite=polite)
    en_b = _resolve_english(query_b, matched_b, guess_b, use_cache=use_cache, polite=polite)
    if not en_a or not en_b:
        return None
    return co_report_stats(en_a, en_b, use_cache=use_cache, polite=polite)


if __name__ == "__main__":
    import sys
    args = sys.argv[1:]
    a, b = (args[0], args[1]) if len(args) >= 2 else ("amlodipine", "diltiazem")
    stats = co_report_stats(a, b)
    print(f"\n=== {stats['drug_a']} + {stats['drug_b']} 併用報告 (openFDA/FAERS) ===")
    print(f"  併用報告総数        : {stats['co_reports_total']:,}")
    print(f"  うち死亡転帰        : {stats['co_reports_death']:,}  ({stats['death_ratio']:.1%})")
    print(f"  うち重篤(serious=1) : {stats['co_reports_serious']:,}  ({stats['serious_ratio']:.1%})")
    print(f"  単剤報告数(参考)    : {stats['drug_a']}={stats['solo_a_total']:,} / {stats['drug_b']}={stats['solo_b_total']:,}")
    print("  頻出有害事象Top5    :")
    for term, cnt in stats["top_reactions"]:
        print(f"    - {term}: {cnt:,}")
