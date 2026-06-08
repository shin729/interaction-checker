# -*- coding: utf-8 -*-
"""
PMDA添付文書PDFから抽出したテキストを画面表示用に整形するユーティリティ。

PyMuPDFのテキスト抽出はPDFの段組み・セル幅で改行されるため、
そのまま<pre>等で表示すると「横幅10文字程度の細い列」が
延々と続くように見えて読みにくい（実際の表組みの縦書き見出しが
1文字ずつ改行される問題も含む）。
本モジュールは句点・括弧・空行を手がかりに行を連結し、
自然な文章の流れに「reflow」し直す。
"""
import re

_FLUSH_SUFFIXES = ("。", "、")
_BRACKET_CLOSE = "）)】］"


def reflow(text: str) -> str:
    """改行で分断された行を、文の区切り（句点等）まで連結して読みやすくする"""
    if not text:
        return text
    lines = [ln.strip() for ln in text.splitlines()]
    paragraphs = []
    buf = ""
    for ln in lines:
        if not ln:
            if buf:
                paragraphs.append(buf)
                buf = ""
            continue
        buf += ln
        if buf.endswith(_FLUSH_SUFFIXES) or (buf and buf[-1] in _BRACKET_CLOSE):
            paragraphs.append(buf)
            buf = ""
    if buf:
        paragraphs.append(buf)
    return "\n".join(paragraphs)


def collapse(text: str) -> str:
    """文中の改行・余分な空白を取り除き、1つの連続したテキストにする（短い引用向け）"""
    if not text:
        return text
    return re.sub(r"\s+", "", text).strip()
