# -*- coding: utf-8 -*-
"""
相互作用チェッカーの精度回帰テスト。

validation_set.json の各ペアを実際のパイプライン（pmda_lookup + openfda_lookup +
severity.classify + pk_numbers）に通し、判定レベルが期待値と一致するか検証する。

結果は3区分:
  ✅ 一致     : 判定が期待値どおり
  ⚠️ 既知ギャップ: 期待値と不一致だが known_gap=true（現状の限界として把握済み）
  ❌ 回帰      : 期待値と不一致かつ known_gap ではない＝デグレード（要修正）

回帰が1件でもあれば exit code 1 を返す（CI/コミット前チェックに使える）。
expect_pk=true のペアでPK数値が抽出できない場合も警告する。

使い方:
  python validate.py            # 全ペアを検証（キャッシュ利用）
  python validate.py --no-cache # キャッシュを使わず再取得して検証
"""
import json
import sys
from pathlib import Path

import checker
import openfda_lookup
import pmda_lookup
import severity

# Windowsコンソール(cp932)でも日本語・記号が化けないようUTF-8出力にする
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

SET_FILE = Path(__file__).parent / "validation_set.json"


def _evaluate(pair, use_cache=True):
    qa, qb = pair["a"], pair["b"]
    pa = pmda_lookup.lookup(qa, use_cache=use_cache)
    pb = pmda_lookup.lookup(qb, use_cache=use_cache)
    if not pa["found"] or not pb["found"]:
        return {"level": "取得失敗", "pk": 0, "ror": None,
                "note": f"添付文書未取得 a={pa['found']} b={pb['found']}"}
    fda = openfda_lookup.lookup_pair(
        qa, pa.get("matched_name") or qa, qb, pb.get("matched_name") or qb, use_cache=use_cache)
    verdict = severity.classify(pa, pb, qa, qb, fda)
    pk = checker._build_pk_changes(pa, pb, qa, qb)
    return {"level": verdict["level"], "pk": len(pk),
            "ror": fda.get("ror") if fda else None, "note": ""}


def _unit_checks():
    """合成データによる判定ロジックの単体ガード。実ペアでは添付文書層が
    ほぼ全ての相互作用を拾うため「弱」の正例が得にくい。弱判定が死んでいない
    こと、およびROR過剰発火が再発しないことを合成fda_statsで固定する。"""
    empty = {"contraindicated_combinations": None, "caution_combinations": None,
             "pk_interactions": None, "matched_name": "X"}
    emptyY = dict(empty, matched_name="Y")

    def fda(ror, terms):
        return {"co_reports_total": 500, "co_reports_death": 0, "ror": ror,
                "ror_ci_low": ror * 0.9, "ror_ci_high": ror * 1.1,
                "top_reactions": [(t, 100) for t in terms]}

    cases = [
        ("相互作用語あり→弱", fda(2.0, ["DRUG INTERACTION", "NAUSEA"]), "弱"),
        ("相互作用語なし高ROR→記載なし(過剰発火ガード)", fda(80.0, ["NAUSEA", "FALL"]), "記載なし"),
        ("openFDAデータ無し→記載なし", None, "記載なし"),
    ]
    print("\n[単体ガード]")
    ok = True
    for label, f, exp in cases:
        got = severity.classify(empty, emptyY, "X", "Y", f)["level"]
        mark = "[OK]  " if got == exp else "[FAIL]"
        if got != exp:
            ok = False
        print(f"  {mark} {label}: 期待[{exp}] 実際[{got}]")
    return ok


def main():
    use_cache = "--no-cache" not in sys.argv
    data = json.loads(SET_FILE.read_text(encoding="utf-8"))
    pairs = data["pairs"]

    match = gaps = regressions = pk_warn = 0
    print(f"\n{'判定':<6}{'期待':<6}{'結果':<14}{'PK':<4}{'ROR':<7} ペア")
    print("-" * 78)
    regression_rows, gap_rows = [], []
    for pair in pairs:
        got = _evaluate(pair, use_cache=use_cache)
        exp = pair["expected_level"]
        is_gap = pair.get("known_gap", False)
        ok = got["level"] == exp
        name = f"{pair['a']} × {pair['b']}"

        if ok:
            mark, status = "[OK]  ", "一致"
            match += 1
        elif is_gap:
            mark, status = "[GAP] ", "既知ギャップ"
            gaps += 1
            gap_rows.append((name, exp, got))
        else:
            mark, status = "[FAIL]", "回帰(要確認)"
            regressions += 1
            regression_rows.append((name, exp, got))

        ror = f"{got['ror']:.1f}" if got["ror"] is not None else "-"
        print(f"{mark} {got['level']:<5}{exp:<6}{status:<14}{got['pk']:<4}{ror:<7} {name}")

        if pair.get("expect_pk") and got["pk"] == 0 and got["level"] != "取得失敗":
            pk_warn += 1
            print(f"         (!) PK数値が期待されるが0件: {name}")

    units_ok = _unit_checks()

    total = len(pairs)
    print("-" * 78)
    print(f"一致 {match}/{total} ｜ 既知ギャップ {gaps} ｜ 回帰 {regressions} ｜ PK欠落警告 {pk_warn}"
          f" ｜ 単体ガード {'OK' if units_ok else 'FAIL'}")

    if regression_rows:
        print("\n[FAIL] 回帰（known_gapに無い不一致＝デグレードの可能性）:")
        for name, exp, got in regression_rows:
            print(f"   - {name}: 期待[{exp}] → 実際[{got['level']}] {got['note']}")
    if gap_rows:
        print("\n[GAP] 既知ギャップ（把握済みの限界。多くはROR弱判定の特異度問題）:")
        for name, exp, got in gap_rows:
            ror = f"ROR{got['ror']:.1f}" if got["ror"] is not None else ""
            print(f"   - {name}: 期待[{exp}] → 実際[{got['level']}] {ror}")

    return 1 if (regressions or not units_ok) else 0


if __name__ == "__main__":
    sys.exit(main())
