# -*- coding: utf-8 -*-
"""
薬力学的(PD)な相加リスクを検出するモジュール。

CYPやP-糖蛋白の薬物動態(PK)相互作用（interaction_predict.py）とは別の軸で、
『同じ方向の作用を持つ薬を併用すると相加的に危険』というタイプの相互作用を拾う。
例: QT延長薬×QT延長薬、セロトニン作動薬同士、抗凝固薬＋NSAID＝出血、
    オピオイド＋ベンゾジアゼピン＝呼吸抑制、RAA系＋K保持性利尿薬＝高カリウム。

添付文書は個々のペアを網羅しないため、この種の相加リスクは記載が漏れやすい。
pd_risks.json の薬剤→リスクフラグ表を使い、両剤が同じフラグを持てば警告する。
これはあくまでスクリーニングであり、実際の危険度は用量・投与経路・患者要因で異なる。
判定(verdict)には反映せず、機序予測と同じく『参考』として提示する。
"""
import json
from pathlib import Path

import severity

_FILE = Path(__file__).parent / "pd_risks.json"
_cache = None  # {"risks": {...}, "index": {薬剤中核名: set(flag)}}


def _data():
    global _cache
    if _cache is None:
        raw = json.loads(_FILE.read_text(encoding="utf-8"))
        risks = raw.get("risks", {})
        index = {}
        for flag, info in risks.items():
            for drug in info.get("drugs", []):
                index.setdefault(drug, set()).add(flag)
        _cache = {"risks": risks, "index": index}
    return _cache


def _flags_for(name):
    """薬剤名（matched_name等）からPDリスクフラグの集合を引く。中核名で前方一致照合。"""
    if not name:
        return set()
    data = _data()
    core = severity._name_core(name)
    for key, flags in data["index"].items():
        if name == key or core == key or name.startswith(key) or core.startswith(key):
            return flags
    return set()


def predict(name_a, name_b):
    """2剤が共有するPDリスクフラグを返す（無ければ空）。

    戻り値: [{"flag", "label", "concern"}, ...]（pd_risks.jsonの並び順で安定化）
    """
    shared = _flags_for(name_a) & _flags_for(name_b)
    if not shared:
        return []
    risks = _data()["risks"]
    order = list(risks.keys())
    out = []
    for flag in sorted(shared, key=order.index):
        info = risks[flag]
        out.append({"flag": flag, "label": info["label"], "concern": info["concern"]})
    return out


def group_shared(names):
    """複数薬のリストから、2剤以上が共有するPDリスクフラグ群を返す（多剤マトリクス用）。

    ポリファーマシーでは、同じ相加リスク（例: QT延長）を持つ薬が処方内に複数あると
    リスクが積み上がる。処方全体を薬剤名リストで受け取り、フラグごとに該当薬をまとめ、
    2剤以上該当するフラグだけを返す。

    戻り値: [{"flag", "label", "concern", "members":[薬剤名...]}, ...]（pd_risks.jsonの順）
    """
    risks = _data()["risks"]
    members = {}
    for name in names:
        for flag in _flags_for(name):
            members.setdefault(flag, []).append(name)
    out = []
    for flag in risks:
        ms = members.get(flag, [])
        if len(ms) >= 2:
            out.append({"flag": flag, "label": risks[flag]["label"],
                        "concern": risks[flag]["concern"], "members": ms})
    return out


if __name__ == "__main__":
    import sys
    pairs = [tuple(sys.argv[1:3])] if len(sys.argv) >= 3 else [
        ("クラリスロマイシン", "ハロペリドール"),   # QT延長×QT延長
        ("セルトラリン", "トラマドール"),           # セロトニン症候群
        ("ワルファリン", "ロキソプロフェン"),       # 出血
        ("オキシコドン", "アルプラゾラム"),         # 中枢・呼吸抑制
        ("エナラプリル", "スピロノラクトン"),       # 高カリウム
        ("アムロジピン", "アトルバスタチン"),       # 相加リスクなし
    ]
    for a, b in pairs:
        ws = predict(a, b)
        print(f"\n=== {a} × {b} ===")
        if not ws:
            print("  相加的なPDリスクの共有フラグなし")
        for w in ws:
            print(f"  [{w['flag']}] {w['label']}")
            print(f"     {w['concern']}")
