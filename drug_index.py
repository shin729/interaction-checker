# -*- coding: utf-8 -*-
"""
ローカル薬剤名インデックス（サジェスト高速化・打ち間違い許容・英語名入力対応）。

drug_name_map.json（日本語一般名→英語INN, 273件）と drug_classes.json
（薬効系統メンバー）から、ローカルに薬剤名の索引を一度だけ構築する。
これにより、よく使う薬についてはPMDAへ毎回問い合わせずに即座に入力候補
（サジェスト）を返せる。日本語のあいまい一致（打ち間違い許容）と、
英語名入力（amlodipine → アムロジピン）にも対応する。

設計上の使い分け（app.py側で組み合わせる）:
  - サジェスト: まず local_suggest()（即時・耐タイポ・英語対応）。ローカル辞書に
    無い薬だけ pmda_lookup.suggest() へフォールバック（網羅性を維持）。PMDAは
    1回約2秒かかるため、一次候補をローカルにするだけで体感が大きく変わる。
  - 本検索の英語名解決: resolve_input()。英語名が辞書に**完全一致**したときだけ
    日本語一般名へ変換する。あいまい一致（打ち間違い補正）は本検索には適用しない
    ――医療ツールで検索対象の薬剤名を黙って別の薬に書き換えるのは危険なため、
    打ち間違い許容はあくまで「候補提示（サジェスト）」に留め、実際に照合する
    薬剤名は利用者が選んだ／入力したものを尊重する。
"""
import difflib
import json
import re
from pathlib import Path

_DIR = Path(__file__).parent
_MAP_FILE = _DIR / "drug_name_map.json"
_CLASSES_FILE = _DIR / "drug_classes.json"

# 英語名入力の判定（ラテン文字・数字・空白・ハイフン・ピリオドのみで構成される語）
_ENGLISH_RE = re.compile(r"[A-Za-z0-9 .\-]+")

_index = None  # 遅延構築（プロセス内で一度だけ）


def _build():
    global _index
    try:
        jp_to_en = json.loads(_MAP_FILE.read_text(encoding="utf-8"))
    except Exception:
        jp_to_en = {}

    jp_names = set(jp_to_en.keys())
    # 薬効系統表のメンバー（drug_name_map に無い一般名も候補に含める）
    try:
        classes = json.loads(_CLASSES_FILE.read_text(encoding="utf-8")).get("classes", [])
        for cls in classes:
            for m in cls.get("members", []):
                if m:
                    jp_names.add(m)
    except Exception:
        pass

    # 英語名（小文字）→ 日本語一般名。同じ英語に複数の日本語が紐づくことは稀だが、
    # 念のため最初に出たものを採用する。
    en_to_jp = {}
    for jp, en in jp_to_en.items():
        if en:
            en_to_jp.setdefault(en.lower(), jp)

    _index = {
        "jp_names": sorted(jp_names),
        "en_to_jp": en_to_jp,
        "en_names": sorted(en_to_jp.keys()),
    }
    return _index


def _idx():
    return _index if _index is not None else _build()


def _is_english(q: str) -> bool:
    return bool(q) and bool(_ENGLISH_RE.fullmatch(q))


def _dedupe(seq, limit):
    out, seen = [], set()
    for x in seq:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
            if len(out) >= limit:
                break
    return out


def _suggest_japanese(q: str, idx: dict, limit: int) -> list:
    names = idx["jp_names"]
    qlen = len(q)

    prefix = sorted((n for n in names if n.startswith(q)), key=lambda n: (len(n), n))
    substr = sorted((n for n in names if q in n and not n.startswith(q)), key=lambda n: (len(n), n))

    # 打ち間違い許容(1): 入力長ぶんの先頭部分が近い名前（例 アモロジ→アムロジピン）
    scored = []
    for n in names:
        if n.startswith(q) or q in n:
            continue
        ratio = difflib.SequenceMatcher(None, q, n[:qlen]).ratio()
        if ratio >= 0.7:
            scored.append((ratio, len(n), n))
    scored.sort(key=lambda t: (-t[0], t[1]))
    fuzzy_prefix = [n for _, _, n in scored]

    # 打ち間違い許容(2): 名前全体としての近似一致（例 クラリスマイシン→クラリスロマイシン）
    fuzzy_full = difflib.get_close_matches(q, names, n=limit, cutoff=0.6)

    return _dedupe(prefix + substr + fuzzy_prefix + fuzzy_full, limit)


def _suggest_english(q: str, idx: dict, limit: int) -> list:
    en_names = idx["en_names"]  # すべて小文字
    en_to_jp = idx["en_to_jp"]

    prefix = sorted((e for e in en_names if e.startswith(q)), key=lambda e: (len(e), e))
    substr = sorted((e for e in en_names if q in e and not e.startswith(q)), key=lambda e: (len(e), e))
    # 英語は同じ語幹（-sartan, -mycin 等）で多数が緩く一致しやすいので、しきい値を高めにして無関係な候補を抑える
    fuzzy = difflib.get_close_matches(q, en_names, n=limit, cutoff=0.7)

    ordered_en, seen = [], set()
    for e in prefix + substr + fuzzy:
        if e not in seen:
            seen.add(e)
            ordered_en.append(e)

    # 候補は日本語一般名で返す（本検索はPMDAを日本語名で引くため、選んでそのまま検索できる）
    return _dedupe((en_to_jp[e] for e in ordered_en if e in en_to_jp), limit)


def local_suggest(query: str, limit: int = 8) -> list:
    """ローカル辞書から入力候補（日本語一般名のリスト）を即時に返す。

    日本語入力は前方一致→部分一致→打ち間違い許容の順、英語入力は英語名で照合して
    対応する日本語一般名へ変換して返す。ローカルに該当が無ければ空リスト
    （呼び出し側でPMDA検索にフォールバックする想定）。
    """
    q = (query or "").strip()
    if len(q) < 2:
        return []
    idx = _idx()
    if _is_english(q):
        return _suggest_english(q.lower(), idx, limit)
    return _suggest_japanese(q, idx, limit)


def resolve_input(query: str) -> str:
    """本検索向けの入力正規化。英語名が辞書に**完全一致**したときだけ日本語一般名へ
    変換する（例: "amlodipine" → "アムロジピン"）。

    あいまい一致はここでは行わない（検索対象の薬剤を黙って書き換えない）。日本語入力・
    未知の英語名はそのまま返す（後段のPMDA検索／グレースフルなエラーに委ねる）。
    """
    q = (query or "").strip()
    if not q or not _is_english(q):
        return q
    return _idx()["en_to_jp"].get(q.lower(), q)


if __name__ == "__main__":
    import sys
    tests = sys.argv[1:] or ["アムロ", "amlo", "クラリスマイシン", "アモロジ", "valsartan", "あ"]
    for t in tests:
        print(f"{t!r:>22} -> local_suggest: {local_suggest(t)}")
    print("--- resolve_input ---")
    for t in ("amlodipine", "AMLODIPINE", "warfarin", "アムロジピン", "unknownXYZ"):
        print(f"{t!r:>14} -> {resolve_input(t)!r}")
