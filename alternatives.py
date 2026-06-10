# -*- coding: utf-8 -*-
"""
薬剤名から「同じ薬効系統の他の薬」を提示するモジュール。

ある薬剤の組み合わせに具体的な数値データが無いとき、同系統の別の薬
（例: ロキソプロフェン→ジクロフェナク・セレコキシブ等）なら数値が
見つかることがあるため、調べ直しの候補として提示する。

系統の定義は drug_classes.json（手動キュレーション）。薬剤名は塩・水和物等の
修飾を落とした中核名（severity._name_core）で照合する。
"""
import json
from pathlib import Path

import severity

_FILE = Path(__file__).parent / "drug_classes.json"
_cache = None


def _classes():
    global _cache
    if _cache is None:
        _cache = json.loads(_FILE.read_text(encoding="utf-8")).get("classes", [])
    return _cache


def _matches(member: str, name: str, core: str) -> bool:
    """系統メンバー名 member が薬剤名 name（中核名 core）に該当するか。
    『ロキソプロフェン』が『ロキソプロフェンナトリウム水和物』に前方一致する等。"""
    return (
        name == member or core == member
        or name.startswith(member) or core.startswith(member)
    )


def find_alternatives(name: str, exclude=()):
    """name が属する薬効系統の、他のメンバーを返す。

    戻り値: {"label": 系統名, "matched": 該当メンバー, "members": [他メンバー...]}
            または該当系統が無ければ None。
    exclude に渡した薬剤名（相手剤など）と中核名が一致するものは候補から除く。
    """
    if not name:
        return None
    core = severity._name_core(name)
    exclude_cores = {severity._name_core(e) for e in exclude if e}
    for cls in _classes():
        for member in cls["members"]:
            if _matches(member, name, core):
                others = [
                    m for m in cls["members"]
                    if m != member and m not in exclude_cores
                ]
                if not others:
                    return None
                return {"label": cls["label"], "matched": member, "members": others}
    return None


if __name__ == "__main__":
    import sys
    for n in (sys.argv[1:] or ["ロキソプロフェンナトリウム水和物", "アムロジピンベシル酸塩", "スボレキサント"]):
        alt = find_alternatives(n)
        print(f"\n{n} ->")
        if alt:
            print(f"  系統: {alt['label']}（該当: {alt['matched']}）")
            print(f"  同系統: {' / '.join(alt['members'])}")
        else:
            print("  系統不明（drug_classes.json 未登録）")
