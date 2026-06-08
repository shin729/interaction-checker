# -*- coding: utf-8 -*-
"""
添付文書の併用禁忌・併用注意の記載から、相互作用の根底にある薬物動態学的な
機序（CYP3A4阻害、グルクロン酸抱合、P-糖蛋白関与、ビタミンK代謝拮抗等）を
キーワードで拾い出し、判定根拠が「どのロジックで照合されたか」を示す
補助情報として提示するモジュール。

PDFの表組みテキストは列単位に抽出されるため、本来は同じ行（同じ薬剤）に
属する「機序・危険因子」列の説明と「薬剤名等」列の薬剤名が、抽出後の
テキストでは数百文字離れて並ぶことがある（実例: クラリスロマイシンの
併用注意欄で「本剤のCYP3Aに対する阻害作用により、左記薬剤の代謝が
阻害される」という説明が、薬剤名「ワルファリンカリウム」の約540文字
手前に位置していた）。一方でワルファリンの併用禁忌欄のように、薬剤名の
直後（数十〜200文字程度）に機序文が来るケースもある。

テーブル構造を厳密に再構成するのは現実的ではないため、相手薬剤名の
出現位置を起点に前方・後方に window を取り、その範囲内に出現する
機序キーワードを拾うヒューリスティックを採る。
"""
import re

_MECHANISM_PATTERNS = [
    ("CYP3A4", re.compile(r"CYP\s*3\s*A")),
    ("CYP2C9", re.compile(r"CYP\s*2\s*C\s*9")),
    ("CYP2C19", re.compile(r"CYP\s*2\s*C\s*19")),
    ("CYP2D6", re.compile(r"CYP\s*2\s*D\s*6")),
    ("CYP1A2", re.compile(r"CYP\s*1\s*A\s*2")),
    ("CYP2B6", re.compile(r"CYP\s*2\s*B\s*6")),
    ("CYP(分子種不明)", re.compile(r"CYP(?!\s*[0-9])")),
    ("肝代謝酵素阻害", re.compile(r"肝薬物代謝酵素|代謝酵素[をが]?阻害")),
    ("P-糖蛋白(P-gp)", re.compile(r"P-?\s*糖[蛋た]ん?白|P-?gp")),
    ("グルクロン酸抱合", re.compile(r"グルクロン酸抱合")),
    ("硫酸抱合", re.compile(r"硫酸抱合")),
    ("血漿蛋白結合競合", re.compile(r"血漿蛋白(結合|からの遊離)")),
    ("トランスポーター(OATP等)", re.compile(r"トランスポーター|OATP")),
    ("血小板凝集への影響", re.compile(r"血小板凝集")),
    ("QT延長", re.compile(r"QT\s*(間隔)?延長")),
    ("中枢神経抑制", re.compile(r"中枢神経(系)?(の)?(抑制|興奮)")),
    ("出血傾向の増強", re.compile(r"出血(傾向|のリスク|の危険)")),
    ("ビタミンK代謝拮抗", re.compile(r"ビタミンK")),
    ("セロトニン症候群", re.compile(r"セロトニン症候群")),
    ("腎排泄・尿細管への影響", re.compile(r"尿細管|腎(での)?排泄")),
    ("消化管吸収への影響", re.compile(r"(消化管|腸管)(から)?の?吸収")),
]

_WINDOW_BEFORE = 650
_WINDOW_AFTER = 250
_MAX_TAGS = 4


def extract_near(text, *names, window_before=_WINDOW_BEFORE, window_after=_WINDOW_AFTER):
    """
    text中で names のいずれかが最初に出現する位置を起点に、その前後の window
    から機序キーワードを拾い出す。重複ラベルを除き、出現順を保って最大
    _MAX_TAGS 件まで返す（CYPの分子種が特定できた場合は「CYP(分子種不明)」を
    重ねて出さない）。
    """
    if not text:
        return []
    collapsed = re.sub(r"\s+", "", text)
    pos = -1
    for n in names:
        if not n or len(n) < 2:
            continue
        idx = collapsed.find(n)
        if idx != -1 and (pos == -1 or idx < pos):
            pos = idx
    if pos == -1:
        return []

    window = collapsed[max(0, pos - window_before):min(len(collapsed), pos + window_after)]

    found = []
    cyp_specific_hit = False
    for label, pattern in _MECHANISM_PATTERNS:
        if not pattern.search(window):
            continue
        if label == "CYP(分子種不明)" and cyp_specific_hit:
            continue
        if label.startswith("CYP") and label != "CYP(分子種不明)":
            cyp_specific_hit = True
        found.append(label)
        if len(found) >= _MAX_TAGS:
            break
    return found


if __name__ == "__main__":
    import pmda_lookup

    cases = [
        ("ワルファリン", "クラリスロマイシン", "caution_combinations"),
        ("ワルファリン", "メナテトレノン", "contraindicated_combinations"),
    ]
    for self_name, other_name, block in cases:
        info = pmda_lookup.lookup(self_name)
        text = info.get(block)
        tags = extract_near(text, other_name)
        print(f"\n=== {self_name}（{block}）× {other_name} ===")
        print(f"  機序タグ: {tags}")
