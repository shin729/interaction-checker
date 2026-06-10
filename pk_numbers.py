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

# AUCは「AUC0-∞」「AUC0-t」「AUCτ」等の下付き表記を伴う（特にインタビューフォーム）。
# 下付き部分に ∞・t・τ・inf も許容する。
_METRIC_RE = r"(Cmax|AUC(?:[0-9\-–〜~tτ∞]|inf|∞)*|Tmax|t\s*1\s*/\s*2|半減期|血中濃度|血漿中濃度|血清中濃度|クリアランス|CL)"
_NUM_RE = r"[0-9０-９]+(?:[.．][0-9０-９]+)?"
# 「1.26倍」のような単一値だけでなく「X〜Y倍」のような範囲表記も拾う
# （添付文書には単一値が無く範囲のみ記載のケースがあるため、ユーザー要望で対応）
_RANGE_RE = rf"({_NUM_RE})(?:\s*[～〜~\-–]\s*({_NUM_RE}))?"
_CONNECT_RE = r"[はがでにと、・約]{0,4}"
# ペア形式の接続部は「は、それぞれ」「は各々」等が入るため、それぞれ/各々も許容する
_PAIRED_CONNECT_RE = r"(?:それぞれ|各々|[はがでにと、・約]){0,8}"
_UNIT_DIR_RE = r"(倍|[%％])(?:程度|に|まで)?(上昇|増加|増強|低下|減少|減弱|延長|短縮)?"

_PAIR_RE = re.compile(rf"{_METRIC_RE}{_CONNECT_RE}{_RANGE_RE}{_UNIT_DIR_RE}")

# 「Cmax及びAUCは22%及び105%増加」「Cmax及びAUCは、それぞれ22%及び105%上昇」のように、
# 指標を2つ並べてから数値を2つまとめて書く形式（旧書式の添付文書・IFに多い）。
# metric1↔num1、metric2↔num2 を順番で対応づける。方向は両指標に共通で末尾に1度だけ。
_PAIRED_RE = re.compile(
    rf"{_METRIC_RE}(?:及び|並びに|、){_METRIC_RE}{_PAIRED_CONNECT_RE}"
    rf"({_NUM_RE})(倍|[%％])(?:及び|並びに|、)({_NUM_RE})(倍|[%％])?"
    rf"(?:程度|に|まで)?(上昇|増加|増強|低下|減少|減弱|延長|短縮)?"
)

# 範囲表記の上限/下限比がこの値以上なら「広すぎて現場判断に使いにくい」低信頼とみなす
# （例: 1〜10倍=比10は弾く。1.5〜2.2倍=比約1.5は通常表示）。ユーザー要望に基づく。
_WIDE_RANGE_RATIO = 3.0

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


def _make_change(metric, value, unit, direction, value_max=None) -> dict:
    unit = unit.translate(str.maketrans("％", "%"))
    v = _to_float(value)
    v_max = _to_float(value_max) if value_max else None
    # 範囲が広すぎる（上限/下限比が_WIDE_RANGE_RATIO以上）と現場判断に使えないため低信頼フラグ
    wide = v_max is not None and v > 0 and (v_max / v) >= _WIDE_RANGE_RATIO
    return {
        "metric": _normalize_metric(metric),
        "value": v,
        "value_max": v_max,
        "unit": unit,
        "direction": direction,
        "wide": wide,
        # 単一値「1.26倍」/ 範囲「X〜Y倍」のどちらでも表示できるラベル
        "value_label": f"{_fmt_num(v)}〜{_fmt_num(v_max)}{unit}" if v_max is not None else f"{_fmt_num(v)}{unit}",
    }


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
        ns = sentence.replace(" ", "")
        collected = []  # (出現位置, change) ― 最後に位置順に並べ替えて文中の語順を保つ
        consumed = []   # _PAIRED_RE が消費した区間。単一パターンの二重マッチを防ぐ
        # ① 「Cmax及びAUCは22%及び105%増加」形式（指標2つ＋数値2つ）を先に拾う
        for m in _PAIRED_RE.finditer(ns):
            m1, m2, n1, u1, n2, u2, direction = m.groups()
            u2 = u2 or u1  # 2つ目に単位が無ければ1つ目を流用（「22%及び105%」等）
            collected.append((m.start(), _make_change(m1, n1, u1, direction)))
            collected.append((m.start() + 1, _make_change(m2, n2, u2, direction)))
            consumed.append(m.span())
        # ② 「Cmaxで1.26倍」形式（指標→数値が隣接）。①が消費した区間は除外する
        for m in _PAIR_RE.finditer(ns):
            if any(s <= m.start() < e for s, e in consumed):
                continue
            metric, value, value_max, unit, direction = m.groups()
            collected.append((m.start(), _make_change(metric, value, unit, direction, value_max)))
        collected.sort(key=lambda x: x[0])
        changes = [c for _, c in collected]
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


def _mentions(sentence: str, names) -> bool:
    """文(sentence)に names のいずれかが含まれるか。PDF抽出由来の空白は無視して照合する。"""
    s = sentence.replace(" ", "").replace("　", "")
    return any(n and len(n) >= 2 and n.replace(" ", "") in s for n in names)


def extract_all(*texts, partner_names=None):
    """複数のテキスト（pk_interactions / caution_combinations 等）から重複なく抽出する。

    partner_names を渡すと、相手剤名（その中核名・配合成分名を含む）を含む文だけに絞る。
    添付文書の相互作用欄は薬剤ごとに記載が並ぶため、フィルタしないと「いま調べている
    ペアとは無関係な第三の薬剤（例: アムロジピン×ジルチアゼムを調べているのに、
    アムロジピン添文中のシンバスタチンとの相互作用）」の数値まで拾ってしまう。
    相手剤名が同じ文に出現しない数値は、ペアの相互作用とは言い切れないため除外する
    （取りこぼすより、無関係な数値を相互作用として誤提示しない方を優先する設計）。
    """
    seen = set()
    merged = []
    for text in texts:
        for item in extract(text):
            if partner_names and not _mentions(item["sentence"], partner_names):
                continue
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
