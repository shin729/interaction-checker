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
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import severity

BASE = "https://api.fda.gov/drug/event.json"
NAME_MAP_FILE = Path(__file__).parent / "drug_name_map.json"
CACHE_DIR = Path(__file__).parent / "cache" / "openfda"
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


def co_report_stats(drug_a_en: str, drug_b_en: str, use_cache: bool = True) -> dict:
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
        "top_reactions": [(有害事象名, 件数), ...]  # 併用報告での頻出事象Top5
      }
    """
    safe_a, safe_b = _escape(drug_a_en), _escape(drug_b_en)
    cache_file = CACHE_DIR / f"{safe_a}__{safe_b}.json"
    if use_cache and cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    base_q = f'patient.drug.medicinalproduct:"{safe_a}"+AND+patient.drug.medicinalproduct:"{safe_b}"'

    total = _total(base_q)
    death = _total(f"{base_q}+AND+patient.reaction.reactionoutcome:5") if total else 0
    serious = _total(f"{base_q}+AND+serious:1") if total else 0
    solo_a = _total(f'patient.drug.medicinalproduct:"{safe_a}"')
    solo_b = _total(f'patient.drug.medicinalproduct:"{safe_b}"')
    top = _top_reactions(base_q) if total else []

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
    }
    cache_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    time.sleep(0.5)  # openFDAサーバへの負荷軽減
    return result


def _resolve_english(query: str, matched_name: str):
    """検索クエリ→添付文書上の正式名→正式名の中核成分名(severity._name_core)の順に
    drug_name_map.json と照合する。

    drug_name_map.json は中核成分名（例: 「アムロジピン」）で登録する運用のため、
    クエリが先発品名（例: 「ノルバスク」）や塩・水和物付きの正式名（例:
    「アムロジピンベシル酸塩」）であっても、正式名を中核名へ正規化することで
    対応表のエントリに到達できるようにする（配合剤の複合名は中核化されないため
    対象外＝従来通り fda_name_missing として扱われる）。
    """
    for candidate in (query, matched_name, severity._name_core(matched_name)):
        en = to_english(candidate) if candidate else None
        if en:
            return en
    return None


def lookup_pair(query_a: str, matched_a: str, query_b: str, matched_b: str, use_cache: bool = True):
    """日本語の薬剤名2つ(検索クエリと添付文書上の正式名)から
    openFDA併用報告統計を取得する高レベル関数。
    drug_name_map.json に対応する英語名が見つからない場合は None を返す。"""
    en_a = _resolve_english(query_a, matched_a)
    en_b = _resolve_english(query_b, matched_b)
    if not en_a or not en_b:
        return None
    return co_report_stats(en_a, en_b, use_cache=use_cache)


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
