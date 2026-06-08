# -*- coding: utf-8 -*-
"""
2剤の添付文書情報(pmda_lookup)とopenFDA併用報告統計(openfda_lookup)から、
相互作用の強さを「強・中・弱・記載なし」で判定するモジュール。

判定方針（添付文書の記載を一次根拠、openFDAの実報告を補助根拠とする）:
  強     = どちらかの添付文書の「併用禁忌」に相手の薬剤名が記載されている
  中     = どちらかの添付文書の「併用注意」に相手の薬剤名が記載されている
  弱     = 添付文書には直接の記載がないが、openFDA(FAERS)の併用報告で
           死亡転帰が一定数報告されている、または「DRUG INTERACTION」が
           頻出有害事象として観察される（添付文書に未記載の実世界シグナル）
  記載なし = 添付文書・openFDAいずれにも有意な記載/シグナルが見当たらない
           （※「相互作用が無いことの確認」ではなく「記載が見つからなかった」
           という意味。添付文書は報告・検討された組み合わせしか載らないため、
           「記載なし」は「安全」を意味しない。badge色をあえて緑にせず中立色に
           しているのもこのため）

※ 添付文書は「報告・検討された組み合わせ」しか列挙しないため、
  名称の直接記載が無くても実際には相互作用が起こり得る
  （ワルファリン+クラリスロマイシンで実際に確認: 添付文書側に
  クラリスロマイシン/マクロライドの記載は無いが、openFDAでは
  「DRUG INTERACTION」「INR増加」が上位に出現）。
  → openFDAを弱シグナルの補助根拠に使う設計はこの実例に基づく。
※ openFDAによる「弱」判定の閾値(死亡20件など)は暫定値。
  実際の薬剤ペアでの検証結果をもとに調整が必要
  （project_interaction_checker.md 参照）。
"""
import re

_SALT_SUFFIX_RE = re.compile(
    r"([\s　]*(?:ピボキシル|アキセチル|プロキセチル|シレキセチル|メドキソミル|サキセチル)|"
    r"塩酸塩|臭化水素酸塩|硫酸塩|リン酸エステル|リン酸塩|クエン酸塩|マレイン酸塩|フマル酸塩|"
    r"コハク酸塩|ベシル酸塩|メシル酸塩|トシル酸塩|酢酸塩|安息香酸塩|パモ酸塩|ヨウ化物|"
    r"ナトリウム|カリウム|カルシウム|水和物|無水物)+$"
)

_OPENFDA_DEATH_THRESHOLD = 20
_OPENFDA_SIGNAL_TERMS = ("DRUG INTERACTION", "TOXICITY TO VARIOUS AGENTS")

_COMBO_SUFFIX_RE = re.compile(r"配合剤?$")


def _name_core(name: str) -> str:
    """「アムロジピンベシル酸塩」->「アムロジピン」「セフカペン ピボキシル塩酸塩水和物」
    ->「セフカペン」のように、塩・水和物・エステル(プロドラッグ)修飾語の表記を
    落とした中核名を返す"""
    return _SALT_SUFFIX_RE.sub("", name).strip() or name


def _split_combo(name: str) -> list:
    """「テルミサルタン・アムロジピンベシル酸塩配合剤」のような配合剤の複合名を
    成分名のリストに分解する（「・」を含まない通常の薬剤名は[name]のまま返す）。

    先発品（ブランド名）で配合剤を検索すると、PMDAの検索結果は個々の成分名を
    「・」で連結した複合名で返ってくるが、相手剤の添付文書は通常テルミサルタン／
    アムロジピンのように成分名単位で言及するため、複合名のままでは照合できない。
    """
    if "・" not in name:
        return [name]
    parts = [p.strip() for p in name.split("・") if p.strip()]
    if parts:
        parts[-1] = _COMBO_SUFFIX_RE.sub("", parts[-1]).strip() or parts[-1]
    return parts


def _expand_names(query: str, matched_name: str) -> tuple:
    """検索クエリ・マッチ名・（配合剤なら）分解した成分名をまとめ、重複排除して返す"""
    seen, result = set(), []
    for n in (query, matched_name, *_split_combo(matched_name)):
        if n and n not in seen:
            seen.add(n)
            result.append(n)
    return tuple(result)


def _text_mentions(text, *names) -> bool:
    if not text:
        return False
    candidates = set()
    for n in names:
        if not n:
            continue
        candidates.add(n)
        candidates.add(_name_core(n))
    return any(len(c) >= 2 and c in text for c in candidates)


def _openfda_signal(fda_stats):
    if not fda_stats or not fda_stats.get("co_reports_total"):
        return None
    death = fda_stats.get("co_reports_death", 0)
    terms = {t for t, _ in fda_stats.get("top_reactions", [])}
    hit_terms = [t for t in _OPENFDA_SIGNAL_TERMS if t in terms]
    if death >= _OPENFDA_DEATH_THRESHOLD or hit_terms:
        return death, hit_terms
    return None


def classify(pmda_a: dict, pmda_b: dict, query_a: str, query_b: str, fda_stats: dict = None) -> dict:
    """
    戻り値: {"level": "強"|"中"|"弱"|"記載なし", "reason": str}
    """
    name_a = pmda_a.get("matched_name") or query_a
    name_b = pmda_b.get("matched_name") or query_b
    names_a = _expand_names(query_a, name_a)
    names_b = _expand_names(query_b, name_b)

    contra_a = pmda_a.get("contraindicated_combinations")
    contra_b = pmda_b.get("contraindicated_combinations")
    if _text_mentions(contra_a, *names_b):
        return {"level": "強", "reason": "添付文書の「併用禁忌」に相手剤の記載あり",
                "mechanism_text": contra_a, "mechanism_names": names_b}
    if _text_mentions(contra_b, *names_a):
        return {"level": "強", "reason": "添付文書の「併用禁忌」に相手剤の記載あり",
                "mechanism_text": contra_b, "mechanism_names": names_a}

    caution_a = pmda_a.get("caution_combinations")
    caution_b = pmda_b.get("caution_combinations")
    if _text_mentions(caution_a, *names_b):
        return {"level": "中", "reason": "添付文書の「併用注意」に相手剤の記載あり",
                "mechanism_text": caution_a, "mechanism_names": names_b}
    if _text_mentions(caution_b, *names_a):
        return {"level": "中", "reason": "添付文書の「併用注意」に相手剤の記載あり",
                "mechanism_text": caution_b, "mechanism_names": names_a}

    signal = _openfda_signal(fda_stats)
    if signal:
        death, hit_terms = signal
        parts = []
        if death >= _OPENFDA_DEATH_THRESHOLD:
            parts.append(f"死亡転帰{death}件")
        if hit_terms:
            parts.append("頻出事象に" + "・".join(hit_terms))
        return {
            "level": "弱",
            "reason": "添付文書に直接記載はないが、openFDA(FAERS)併用報告で"
                      + "・".join(parts) + "（添付文書未記載の実世界シグナル）",
        }

    return {"level": "記載なし", "reason": "添付文書・openFDA併用報告のいずれにも有意な記載/シグナルが見当たらない"
                                       "（注: 相互作用が無いことの確認ではなく、報告・記載が見つからなかったという意味です）"}


if __name__ == "__main__":
    import sys
    import pmda_lookup
    import openfda_lookup

    pairs = [tuple(sys.argv[1:3])] if len(sys.argv) >= 3 else [
        ("アムロジピン", "ジルチアゼム"),
        ("ワルファリン", "クラリスロマイシン"),
    ]
    for qa, qb in pairs:
        pa, pb = pmda_lookup.lookup(qa), pmda_lookup.lookup(qb)
        fda = openfda_lookup.lookup_pair(qa, qb)
        result = classify(pa, pb, qa, qb, fda)
        print(f"\n=== {qa} × {qb} ===")
        print(f"  添付文書: {pa['matched_name']} / {pb['matched_name']}")
        print(f"  判定: 【{result['level']}】")
        print(f"  根拠: {result['reason']}")
        if fda:
            print(f"  openFDA併用報告: 総数{fda['co_reports_total']:,} / 死亡{fda['co_reports_death']:,}"
                  f"({fda['death_ratio']:.1%}) / 重篤{fda['co_reports_serious']:,}({fda['serious_ratio']:.1%})")
