# -*- coding: utf-8 -*-
"""
CYP/トランスポーターの役割（cyp_roles.json）から、薬剤ペアの相互作用を
機序的に「予測」するモジュール。

添付文書・IF・openFDAは「報告・検討された組み合わせ」しか拾えないため、
記載が無い＝安全とは限らない。一方、薬AがCYP3A4の感受性基質で薬Bが
CYP3A4の強い阻害薬なら、両者の添付文書が相手に言及していなくても
「AのAUCが大きく上昇しうる」と機序から予測できる（Lexicomp等の臨床
DDIチェッカーと同じ考え方）。

これはあくまで予測であり、文書化された根拠とは区別して提示する。
"""
import json
from pathlib import Path

import severity

_FILE = Path(__file__).parent / "cyp_roles.json"
_cache = None

# 阻害薬の強さ×基質の感受性から、予測される阻害の程度（FDA区分の語彙に対応）
_INHIBIT_LEVEL = {
    ("強", True): "強い", ("強", False): "中等度",
    ("中", True): "中等度", ("中", False): "弱い",
    ("弱", True): "弱い", ("弱", False): "弱い",
}
_INDUCE_LEVEL = {
    ("強", True): "強い", ("強", False): "中等度",
    ("中", True): "中等度", ("中", False): "弱い",
    ("弱", True): "弱い", ("弱", False): "弱い",
}

# 強さ→副詞形（「強く阻害」「中程度に阻害」「弱く阻害」）
_ADVERB = {"強": "強く", "中": "中程度に", "弱": "弱く"}


def _roles_table():
    global _cache
    if _cache is None:
        _cache = json.loads(_FILE.read_text(encoding="utf-8")).get("roles", {})
    return _cache


def _roles_for(name):
    """薬剤名（matched_name等）からCYP役割を引く。中核名で前方一致照合。"""
    if not name:
        return None, None
    core = severity._name_core(name)
    for key, roles in _roles_table().items():
        if name == key or core == key or name.startswith(key) or core.startswith(key):
            return key, roles
    return None, None


def _pairwise(victim, vroles, perp, proles):
    """victim が基質、perp が阻害/誘導薬となる相互作用を予測する。"""
    out = []
    substrate = vroles.get("substrate", {})
    inhibitor = proles.get("inhibitor", {})
    inducer = proles.get("inducer", {})
    for enz, sens in substrate.items():
        sensitive = (sens == "感受性")
        if enz in inhibitor:
            level = _INHIBIT_LEVEL[(inhibitor[enz], sensitive)]
            out.append({
                "victim": victim, "perpetrator": perp, "enzyme": enz, "kind": "阻害",
                "level": level,
                "effect": f"{victim}の血中濃度(AUC)が上昇し、作用・副作用が強まる可能性",
                "basis": f"{perp}は{enz}を{_ADVERB[inhibitor[enz]]}阻害、{victim}は{enz}の"
                         f"{'感受性基質' if sensitive else '基質'}",
            })
        if enz in inducer:
            level = _INDUCE_LEVEL[(inducer[enz], sensitive)]
            out.append({
                "victim": victim, "perpetrator": perp, "enzyme": enz, "kind": "誘導",
                "level": level,
                "effect": f"{victim}の血中濃度(AUC)が低下し、効果が減弱する可能性",
                "basis": f"{perp}は{enz}を{_ADVERB[inducer[enz]]}誘導、{victim}は{enz}の"
                         f"{'感受性基質' if sensitive else '基質'}",
            })
    return out


def predict(name_a, name_b):
    """2剤の薬剤名から機序的な相互作用予測のリストを返す（無ければ空）。

    戻り値: [{victim, perpetrator, enzyme, kind(阻害/誘導), level(強い/中等度/弱い),
              effect, basis}, ...]
    """
    ka, ra = _roles_for(name_a)
    kb, rb = _roles_for(name_b)
    if not ra or not rb:
        return []
    preds = _pairwise(ka, ra, kb, rb) + _pairwise(kb, rb, ka, ra)
    # 強い→中等度→弱いの順に並べ、見やすくする
    order = {"強い": 0, "中等度": 1, "弱い": 2}
    preds.sort(key=lambda p: order.get(p["level"], 3))
    return preds


if __name__ == "__main__":
    import sys
    pairs = [tuple(sys.argv[1:3])] if len(sys.argv) >= 3 else [
        ("トリアゾラム", "イトラコナゾール"),
        ("シンバスタチン", "クラリスロマイシン"),
        ("ジゴキシン", "ベラパミル"),
        ("ワルファリン", "フルコナゾール"),
        ("スボレキサント", "リファンピシン"),
    ]
    for a, b in pairs:
        print(f"\n=== {a} × {b} ===")
        for p in predict(a, b):
            print(f"  [{p['level']}{p['kind']}] {p['effect']}")
            print(f"     根拠: {p['basis']}")
        if not predict(a, b):
            print("  予測なし（役割表に未登録 or 同一酵素の基質×阻害/誘導の関係なし）")
