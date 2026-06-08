# -*- coding: utf-8 -*-
"""
verify_pmda_access.py で取得した添付文書PDFのテキストから、
「相互作用」章（併用禁忌・併用注意を含む）を切り出すパーサ。

検証で判明した通り、添付文書には新旧2つの書式が混在している:
  - 新書式（2024年新記載要領）: 番号付き見出し
        10. 相互作用 / 10.1 併用禁忌（併用しないこと） /
        10.2 併用注意（併用に注意すること） / 16.7 薬物相互作用
    例: ワルファリン、クラリスロマイシン
  - 旧書式（番号なし）: 「相互作用」が単独見出しとして現れ、
        直後に「併用禁忌（併用しないこと）」または
        「併用注意（併用に注意すること）」の表が続く。
        次の大見出し（副作用／薬物動態／適用上の注意 等）まで。
    例: ロキソプロフェン（まだ新書式に改訂されていない）
両方を検出し、フォールバックする。
"""
import json
import re
from pathlib import Path

SAMPLE_DIR = Path(__file__).parent / "pmda_sample"

NEW_FORMAT = {
    "interactions": r"10\s*[\.．]\s*相\s*互\s*作\s*用",
    "contraindicated_combo": r"10\s*[\.．]\s*1\s*併\s*用\s*禁\s*忌",
    "caution_combo": r"10\s*[\.．]\s*2\s*併\s*用\s*注\s*意",
    "end": r"11\s*[\.．]\s*副\s*作\s*用",
}

OLD_FORMAT = {
    "interactions": r"\n[ \t　]*相\s*互\s*作\s*用[ \t　]*\n",
    # 見出し直前の改行は「相互作用」見出し側に消費されている場合があるため
    # 先頭の \n を必須にせず独立検索する
    "contraindicated_combo": r"併\s*用\s*禁\s*忌\s*（\s*併\s*用\s*し\s*な\s*い\s*こ\s*と\s*）",
    "caution_combo": r"併\s*用\s*注\s*意\s*（\s*併\s*用\s*に\s*注\s*意\s*す\s*る\s*こ\s*と\s*）",
    "end": r"\n[ \t　]*(?:副作用|薬物動態|臨床成績|薬効薬理|適用上の注意|取扱い上の注意)[ \t　]*\n",
}

PK_INTERACTIONS = r"16\s*[\.．]\s*7\s*[\.．]?\s*薬\s*物\s*相\s*互\s*作\s*用"


def _search(pattern, text, start=0):
    m = re.search(pattern, text[start:])
    return (start + m.start(), start + m.end()) if m else None


def _extract(text, fmt):
    """
    注意: PyMuPDFのテキスト抽出は、添付文書の表（セル位置で組まれている）の
    読み取り順を完全には保証しない。実際にワルファリンのPDFでは
    「10.2 併用注意」の表が「10. 相互作用」の見出しより前の文字位置に
    現れる（テーブルレイアウト由来の抽出順の乱れ）。
    そのため "見出しの出現順" を前提にせず、各小見出しを独立に
    全文検索し、それらの位置から外接区間を逆算する。
    """
    interactions = _search(fmt["interactions"], text)
    if interactions is None:
        return None

    contraindicated = _search(fmt["contraindicated_combo"], text)
    caution = _search(fmt["caution_combo"], text)

    anchors = [interactions[0]]
    if contraindicated:
        anchors.append(contraindicated[0])
    if caution:
        anchors.append(caution[0])
    span_start, span_last = min(anchors), max(anchors)

    end = _search(fmt["end"], text, span_last)
    sec_end = end[0] if end else span_last + 4000

    contraindicated_block = None
    caution_block = None
    if contraindicated:
        c_end = caution[0] if (caution and caution[0] > contraindicated[0]) else \
            (sec_end if contraindicated[0] <= sec_end else contraindicated[0] + 2000)
        contraindicated_block = text[contraindicated[0]:c_end].strip()
    if caution:
        c_end = sec_end if caution[0] <= sec_end else caution[0] + 2500
        caution_block = text[caution[0]:c_end].strip()

    return {
        "interactions_section": text[span_start:sec_end].strip(),
        "contraindicated_combinations": contraindicated_block,
        "caution_combinations": caution_block,
    }


def parse(text: str) -> dict:
    result = _extract(text, NEW_FORMAT)
    fmt_name = "new"
    if result is None:
        result = _extract(text, OLD_FORMAT)
        fmt_name = "old"
    if result is None:
        result = {
            "interactions_section": None,
            "contraindicated_combinations": None,
            "caution_combinations": None,
        }
        fmt_name = "unknown"

    pk = _search(PK_INTERACTIONS, text)
    pk_block = None
    if pk:
        end = _search(NEW_FORMAT["end"], text, pk[1]) or _search(OLD_FORMAT["end"], text, pk[1])
        pk_end = end[0] if end else pk[0] + 4000
        pk_block = text[pk[0]:pk_end].strip()

    result["format"] = fmt_name
    result["pk_interactions"] = pk_block
    return result


if __name__ == "__main__":
    for txt_path in sorted(SAMPLE_DIR.glob("*_tenpu.txt")):
        text = txt_path.read_text(encoding="utf-8")
        parsed = parse(text)
        out_path = txt_path.with_name(txt_path.stem.replace("_tenpu", "_parsed") + ".json")
        out_path.write_text(json.dumps(parsed, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"\n=== {txt_path.stem.replace('_tenpu','')}  [format={parsed['format']}] -> {out_path.name} ===")
        for k, v in parsed.items():
            if k == "format":
                continue
            if v:
                preview = re.sub(r"\s+", " ", v)[:160]
                print(f"  [{k}] ({len(v)} chars): {preview} ...")
            else:
                print(f"  [{k}] (見つからず)")
