# -*- coding: utf-8 -*-
"""
添付文書の「16.7 薬物相互作用」「併用注意」等のテキストから、
血中濃度・AUC・Cmax・クリアランス等の定量的な変化（倍率・%）に
言及している文を抽出するモジュール。

このツールの核となる目的: 併用時に薬物動態パラメータが
「具体的にどれくらい変化するか」は添付文書ではほとんど言及されず、
書籍等でも横断的に検索するのが難しい。それを拾い上げて数値で
見やすく提示することで、添付文書・書籍では分からない部分を補う。

実例（クラリスロマイシン 16.7.2 テオフィリン）:
  「テオフィリンの血清中濃度はCmaxで1.26倍、AUCで1.19倍上昇し、
    クリアランスは16.4%減少した」
  → [{"metric":"Cmax","value":1.26,"unit":"倍","direction":"上昇"}, ...]
"""
import re

_METRIC_RE = r"(Cmax|AUC(?:[0-9\-–〜~]*)?|Tmax|t\s*1\s*/\s*2|半減期|血中濃度|血漿中濃度|血清中濃度|クリアランス|CL)"
_NUM_RE = r"[0-9０-９]+(?:[.．][0-9０-９]+)?"
# 「1.26倍」のような単一値だけでなく「X〜Y倍」のような範囲表記も拾う
# （添付文書には単一値が無く範囲のみ記載のケースがあるため、ユーザー要望で対応）
_RANGE_RE = rf"({_NUM_RE})(?:\s*[～〜~\-–]\s*({_NUM_RE}))?"
_CONNECT_RE = r"[はがでにと、・約]{0,4}"
_UNIT_DIR_RE = r"(倍|[%％])(?:程度|に|まで)?(上昇|増加|増強|低下|減少|減弱|延長|短縮)?"

_PAIR_RE = re.compile(rf"{_METRIC_RE}{_CONNECT_RE}{_RANGE_RE}{_UNIT_DIR_RE}")

_ZEN2HAN = str.maketrans("０１２３４５６７８９．", "0123456789.")


def _to_float(s: str) -> float:
    return float(s.translate(_ZEN2HAN))


def _fmt_num(v: float) -> str:
    """10.0 -> '10', 1.26 -> '1.26' のように余分な小数点以下のゼロを除く"""
    return f"{v:g}"


def _normalize_metric(m: str) -> str:
    m = re.sub(r"\s+", "", m)
    if m.startswith("AUC"):
        return "AUC"
    if m.startswith("t1") or "半減期" in m:
        return "半減期(t1/2)"
    return m


def extract(text: str):
    """
    テキストから定量的なPK変化の言及を文単位で抽出する。

    戻り値: [{"sentence": 該当文, "changes": [{"metric","value","unit","direction"}, ...]}]
    1文に複数の数値言及（Cmax/AUC/クリアランス等）があればまとめて返す。
    PDFのテキスト抽出で生じる文中の改行・余分な空白は除去して判定する。
    """
    if not text:
        return []
    results = []
    for raw_sentence in re.split(r"(?<=。)", text):
        sentence = re.sub(r"[ \t　]*\n[ \t　]*", "", raw_sentence)
        sentence = re.sub(r"\s+", " ", sentence).strip()
        if not sentence:
            continue
        changes = []
        for m in _PAIR_RE.finditer(sentence.replace(" ", "")):
            metric, value, value_max, unit, direction = m.groups()
            unit = unit.translate(str.maketrans("％", "%"))
            v = _to_float(value)
            v_max = _to_float(value_max) if value_max else None
            changes.append({
                "metric": _normalize_metric(metric),
                "value": v,
                "value_max": v_max,
                "unit": unit,
                "direction": direction,
                # 単一値「1.26倍」/ 範囲「X〜Y倍」のどちらでも表示できるラベル
                "value_label": f"{_fmt_num(v)}〜{_fmt_num(v_max)}{unit}" if v_max is not None else f"{_fmt_num(v)}{unit}",
            })
        # 「Cmaxで1.26倍、AUCで1.19倍上昇し」のように方向(上昇/低下等)が
        # 列挙の最後にしか書かれない場合、後続の値から方向を継承する
        for i in range(len(changes) - 2, -1, -1):
            if changes[i]["direction"] is None and changes[i + 1]["direction"]:
                changes[i]["direction"] = changes[i + 1]["direction"]
        for c in changes:
            c["direction"] = c["direction"] or "変化"
        if changes:
            results.append({"sentence": sentence, "changes": changes})
    return results


def extract_all(*texts):
    """複数のテキスト（pk_interactions / caution_combinations 等）から重複なく抽出する"""
    seen = set()
    merged = []
    for text in texts:
        for item in extract(text):
            if item["sentence"] in seen:
                continue
            seen.add(item["sentence"])
            merged.append(item)
    return merged


if __name__ == "__main__":
    import sys
    import pmda_lookup

    names = sys.argv[1:] or ["クラリスロマイシン", "アムロジピン", "ワルファリン", "ジルチアゼム", "ロキソプロフェン"]
    for name in names:
        info = pmda_lookup.lookup(name)
        items = extract_all(info.get("pk_interactions"), info.get("caution_combinations"),
                            info.get("contraindicated_combinations"))
        print(f"\n=== {name}（{info['matched_name']}）数値PKデータ: {len(items)}件 ===")
        for item in items:
            tags = " / ".join(f"{c['metric']} {c['value_label']}{c['direction']}" for c in item["changes"])
            print(f"  [{tags}]")
            print(f"    {item['sentence'][:150]}")
