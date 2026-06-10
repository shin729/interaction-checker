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

import pk_numbers

_SALT_SUFFIX_RE = re.compile(
    r"([\s　]*(?:ピボキシル|アキセチル|プロキセチル|シレキセチル|メドキソミル|サキセチル)|"
    r"塩酸塩|臭化水素酸塩|硫酸塩|リン酸エステル|リン酸塩|クエン酸塩|マレイン酸塩|フマル酸塩|"
    r"コハク酸塩|ベシル酸塩|メシル酸塩|トシル酸塩|酢酸塩|安息香酸塩|パモ酸塩|ヨウ化物|"
    r"ナトリウム|カリウム|カルシウム|水和物|無水物)+$"
)

# 「弱」シグナルは、併用報告の頻出有害事象に相互作用そのものを示す語が
# 出現するかで判定する。
#
# 判定閾値の変遷:
#   (1) 旧: 死亡20件以上の絶対数 → 頻用薬ほど無条件に超える
#   (2) ROR導入: CI下限>1 → 併用処方されるペアをほぼ全部拾ってしまう（特異度が低い）
#   (3) 現行: 相互作用語の出現 → これに変更
# validation_set.json での検証で、RORは真の相互作用(ワルファリン×クラリス ROR3.9等)と
# 単なる併用処方の偽陽性(レボフロキサシン×アムロジピン ROR4.8等)が値域で完全に重なり
# 区別できないのに対し、「DRUG INTERACTION」語は両者を完全分離できることが判明したため、
# 弱判定のトリガーを語の出現に変えた（RORは表示用の参考値として残す）。
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


def _name_candidates(*names) -> set:
    """薬剤名群を、その中核名(塩・水和物等を落とした名)も含めた照合候補集合にする"""
    candidates = set()
    for n in names:
        if not n:
            continue
        candidates.add(n)
        candidates.add(_name_core(n))
    return candidates


def _text_mentions(text, *names) -> bool:
    if not text:
        return False
    return any(len(c) >= 2 and c in text for c in _name_candidates(*names))


def _pk_interaction_with(pk_text, *names) -> bool:
    """16.7薬物相互作用テキストに、相手剤名に紐づく「定量的なPK変化」があるか。

    16.7は相互作用がある薬剤だけでなく「影響を及ぼさなかった」薬剤も列挙する
    （例: シタグリプチンの16.7に『以下の薬物の薬物動態に明らかな影響を及ぼさ
    なかった…メトホルミン…』）。そのため相手剤名が出るだけでは相互作用ありとは
    言えない。相手剤名を含む文に AUC/Cmax 等の数値変化が抽出できる場合のみ真とし、
    「影響なし」記載を相互作用ありと誤判定しないようにする。"""
    if not pk_text:
        return False
    return bool(pk_numbers.extract_all(pk_text, partner_names=_name_candidates(*names)))


def _openfda_signal(fda_stats):
    """openFDA併用報告から弱シグナルの有無を判定する。

    成立条件: 頻出有害事象(Top5)に「DRUG INTERACTION」等の相互作用語が出現する
    ＝報告者が相互作用そのものを有害事象として挙げている。RORは値域が真偽で
    重なり区別できないため判定には使わず、参考値として戻り値に含めるのみ。

    戻り値: シグナルの根拠を説明する辞書（無ければ None）
    """
    if not fda_stats or not fda_stats.get("co_reports_total"):
        return None
    terms = {t for t, _ in fda_stats.get("top_reactions", [])}
    hit_terms = [t for t in _OPENFDA_SIGNAL_TERMS if t in terms]
    if not hit_terms:
        return None
    return {"ror": fda_stats.get("ror"), "ci_low": fda_stats.get("ror_ci_low"),
            "ci_high": fda_stats.get("ror_ci_high"), "death": fda_stats.get("co_reports_death", 0),
            "total": fda_stats.get("co_reports_total", 0), "hit_terms": hit_terms}


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

    # 「16.7 薬物相互作用」セクションも判定根拠に使う。ここは併用禁忌/併用注意の
    # 根拠となるPKデータが載る節で、相手剤名があれば薬物動態学的相互作用が文書化
    # されている＝併用注意相当(中)とみなす。特に旧書式の多段組PDFでは、相手剤が
    # 10.2併用注意から参照されていてもPyMuPDFの抽出順が乱れて併用注意ブロックに
    # 入らないことがあり(例: スボレキサント×ジルチアゼム)、16.7は連続ブロックで
    # 取れるためここを見ないと「PK欄に数値が出ているのに判定は記載なし」という
    # 不整合が起きる。具体的な変化量(AUC何倍等)はPK欄に別途数値表示される。
    pk_a = pmda_a.get("pk_interactions")
    pk_b = pmda_b.get("pk_interactions")
    if _pk_interaction_with(pk_a, *names_b):
        return {"level": "中", "reason": "添付文書の「薬物相互作用(16.7)」に相手剤との定量的な薬物動態変化の記載あり",
                "mechanism_text": pk_a, "mechanism_names": names_b}
    if _pk_interaction_with(pk_b, *names_a):
        return {"level": "中", "reason": "添付文書の「薬物相互作用(16.7)」に相手剤との定量的な薬物動態変化の記載あり",
                "mechanism_text": pk_b, "mechanism_names": names_a}

    signal = _openfda_signal(fda_stats)
    if signal:
        parts = ["頻出有害事象に" + "・".join(signal["hit_terms"]) + "が出現"]
        if signal["ror"] is not None:
            parts.append(f"併用報告ROR {signal['ror']:.1f}（参考値）")
        if signal["death"]:
            parts.append(f"うち死亡転帰{signal['death']}件")
        return {
            "level": "弱",
            "reason": "添付文書に直接記載はないが、openFDA(FAERS)併用報告で"
                      + "、".join(parts)
                      + "（添付文書未記載の実世界の相互作用シグナル）",
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
        fda = openfda_lookup.lookup_pair(
            qa, pa.get("matched_name") or qa, qb, pb.get("matched_name") or qb)
        result = classify(pa, pb, qa, qb, fda)
        print(f"\n=== {qa} × {qb} ===")
        print(f"  添付文書: {pa['matched_name']} / {pb['matched_name']}")
        print(f"  判定: 【{result['level']}】")
        print(f"  根拠: {result['reason']}")
        if fda:
            print(f"  openFDA併用報告: 総数{fda['co_reports_total']:,} / 死亡{fda['co_reports_death']:,}"
                  f"({fda['death_ratio']:.1%}) / 重篤{fda['co_reports_serious']:,}({fda['serious_ratio']:.1%})")
