# -*- coding: utf-8 -*-
"""
多剤マトリクス（checker.check_matrix / matrix_openfda_signals）のオフライン単体テスト。

check_matrix は添付文書取得(pmda_lookup)と判定(severity)に依存するため、ネットワークを
避けてロジック（正規化・重複排除・found/not_found分割・三角ペア生成・needs_openfda判定・
上限/下限のエラー）だけを固定できるよう、依存をフェイクへ差し替えて検証する。
差し替えは各テスト内で必ず元へ戻す。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import checker

_failures = []


def check(label, got, expected):
    ok = got == expected
    print(f"  {'[OK]  ' if ok else '[FAIL]'} {label}: 期待[{expected!r}] 実際[{got!r}]")
    if not ok:
        _failures.append(label)


def check_true(label, cond):
    print(f"  {'[OK]  ' if cond else '[FAIL]'} {label}")
    if not cond:
        _failures.append(label)


class _patch:
    """checker のモジュール属性を一時的に差し替えるコンテキストマネージャ。"""
    def __init__(self, **kw):
        self.kw = kw
        self.old = {}

    def __enter__(self):
        for k, v in self.kw.items():
            self.old[k] = getattr(checker, k)
            setattr(checker, k, v)
        return self

    def __exit__(self, *a):
        for k, v in self.old.items():
            setattr(checker, k, v)


def _fake_pmda(found_names):
    """found_names に含まれる名前だけ found=True を返す偽の pmda_lookup 代替。"""
    class _FakePmda:
        @staticmethod
        def lookup(q, use_cache=True, polite=True):
            return {"found": q in found_names, "matched_name": q,
                    "contraindicated_combinations": None, "caution_combinations": None,
                    "pk_interactions": None, "english_name_guess": None}
    return _FakePmda


def _fake_severity(level_for):
    """ペア(frozenset)→レベル を引く偽の severity 代替（classifyのみ使用）。"""
    class _FakeSev:
        @staticmethod
        def classify(pa, pb, qa, qb, fda):
            lvl = level_for.get(frozenset((qa, qb)), "記載なし")
            return {"level": lvl, "reason": f"test:{lvl}"}
    return _FakeSev


def test_errors_offline():
    print("\n[check_matrix] 入力エラー（ネットワーク不要）")
    check("空入力→2つ以上", checker.check_matrix([])["error"], "薬剤名を2つ以上入力してください。")
    check("1剤→2つ以上", checker.check_matrix(["アムロジピン"])["error"], "薬剤名を2つ以上入力してください。")
    over = checker.check_matrix([f"薬{i}" for i in range(11)])["error"]
    check_true("11剤→上限エラー", over is not None and "10剤" in over)


def test_assembly():
    print("\n[check_matrix] 組み立て（依存をフェイク化）")
    names = ["A", "A", "B", "C", "UNK"]  # A重複、UNKは未収載
    levels = {frozenset(("A", "B")): "強", frozenset(("A", "C")): "中"}  # 他は記載なし
    with _patch(pmda_lookup=_fake_pmda({"A", "B", "C"}), severity=_fake_severity(levels)):
        d = checker.check_matrix(names)
    check("エラーなし", d["error"], None)
    check("重複排除して3剤", [x["input"] for x in d["drugs"]], ["A", "B", "C"])
    check("未収載はnot_found", d["not_found"], ["UNK"])
    check("ペア数=C(3,2)=3", len(d["cells"]), 3)
    lv = {frozenset((c["a"], c["b"])): c["level"] for c in d["cells"]}
    check("A×B=強", lv[frozenset(("A", "B"))], "強")
    check("A×C=中", lv[frozenset(("A", "C"))], "中")
    check("B×C=記載なし", lv[frozenset(("B", "C"))], "記載なし")
    need = {frozenset((c["a"], c["b"])) for c in d["cells"] if c["needs_openfda"]}
    check("openFDA遅延対象は記載なしのB×Cのみ", need, {frozenset(("B", "C"))})
    # 三角インデックスが a_idx < b_idx を満たす
    check_true("a_idx<b_idxを満たす", all(c["a_idx"] < c["b_idx"] for c in d["cells"]))


def test_openfda_signals_shape():
    print("\n[matrix_openfda_signals] 並列取得の形状（pair_openfda_signalをフェイク化）")
    def fake_signal(a, b):
        # A×B は弱シグナルあり、それ以外はなし
        if frozenset((a, b)) == frozenset(("A", "B")):
            return {"level": "弱", "ror": 3.1, "available": True}
        return {"level": "記載なし", "ror": None, "available": True}
    with _patch(pair_openfda_signal=fake_signal):
        pairs = [{"a": "A", "b": "B", "a_idx": 0, "b_idx": 1},
                 {"a": "B", "b": "C", "a_idx": 1, "b_idx": 2}]
        res = checker.matrix_openfda_signals(pairs)
    check("件数はペア数と一致", len(res), 2)
    byidx = {(r["a_idx"], r["b_idx"]): r for r in res}
    check("A×B(0,1)は弱", byidx[(0, 1)]["level"], "弱")
    check("A×B ROR保持", byidx[(0, 1)]["ror"], 3.1)
    check("B×C(1,2)は記載なし", byidx[(1, 2)]["level"], "記載なし")
    check("空入力→空リスト", checker.matrix_openfda_signals([]), [])


def main():
    print("=" * 70)
    print("多剤マトリクスの単体テスト（オフライン・依存フェイク化）")
    print("=" * 70)
    for fn in (test_errors_offline, test_assembly, test_openfda_signals_shape):
        fn()
    print("\n" + "-" * 70)
    if _failures:
        print(f"[FAIL] {len(_failures)}件の失敗: " + "、".join(_failures))
        return 1
    print("[OK] 全テスト通過")
    return 0


if __name__ == "__main__":
    sys.exit(main())
