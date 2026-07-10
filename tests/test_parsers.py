# -*- coding: utf-8 -*-
"""
パース系純粋関数のオフライン単体テスト（ネットワーク不要・stdlibのみ）。

このツールの心臓部＝「添付文書テキストからAUC等の数値変化を抽出し、FDA区分へ翻訳する」
(pk_numbers)、「併用報告からRORを計算する」(openfda_lookup._ror)、「新旧書式の添付文書から
相互作用の章を切り出す」(parse_tenpu) の3系統を、合成入力で固定する。いずれも決定論的な
純粋関数で、正規表現・境界値・数式のリグレッションをネットワークなしで即検出する。

  python tests/test_parsers.py   # 失敗が1件でもあれば exit 1
"""
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import openfda_lookup as fda
import parse_tenpu as pt
import pk_numbers as pk

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


def _cat(metric, value, unit, direction, vmax=None):
    """_make_change→_fda_auc_category の合成ヘルパ（fdaキーはextract側で付くため直接呼ぶ）"""
    c = pk._make_change(metric, str(value), unit, direction, str(vmax) if vmax is not None else None)
    return pk._fda_auc_category(c)


# ---------- pk_numbers ----------

def test_extract_multi_metric():
    print("\n[pk_numbers.extract] 1文に複数指標（Cmax/AUC/クリアランス）")
    # 実例（クラリスロマイシン 16.7.2 テオフィリン）
    s = "テオフィリンの血清中濃度はCmaxで1.26倍、AUCで1.19倍上昇し、クリアランスは16.4%減少した。"
    items = pk.extract(s)
    check("1文にまとまる", len(items), 1)
    got = [(c["metric"], c["value_label"], c["direction"]) for c in items[0]["changes"]]
    check("3指標を語順どおり抽出",
          got, [("Cmax", "1.26倍", "上昇"), ("AUC", "1.19倍", "上昇"), ("クリアランス", "16.4%", "減少")])
    # AUC 1.19倍は1.25倍未満なのでFDA区分は付かない
    check("AUC1.19倍はFDA区分なし", items[0]["changes"][1]["fda"], None)


def test_extract_paired_form():
    print("\n[pk_numbers.extract] 指標2つ＋数値2つの併記形式")
    s = "本剤併用によりCmax及びAUCは、それぞれ22%及び105%上昇した。"
    items = pk.extract(s)
    got = [(c["metric"], c["value_label"], c["direction"], c["fda"]) for c in items[0]["changes"]]
    check("Cmax22%↑・AUC105%↑（AUCは中等度阻害相当）",
          got, [("Cmax", "22%", "上昇", None),
                ("AUC", "105%", "上昇", "中等度阻害相当(AUC2〜5倍)")])


def test_fda_auc_category():
    print("\n[pk_numbers._fda_auc_category] AUC倍率→FDA相互作用区分（境界値）")
    check("5倍→強い阻害", _cat("AUC", 5, "倍", "上昇"), "強い阻害相当(AUC5倍以上)")
    check("3倍→中等度阻害", _cat("AUC", 3, "倍", "上昇"), "中等度阻害相当(AUC2〜5倍)")
    check("1.5倍→弱い阻害", _cat("AUC", 1.5, "倍", "上昇"), "弱い阻害相当(AUC1.25〜2倍)")
    check("1.19倍→区分なし(1.25未満)", _cat("AUC", 1.19, "倍", "上昇"), None)
    check("105%上昇→中等度(fold2.05)", _cat("AUC", 105, "%", "上昇"), "中等度阻害相当(AUC2〜5倍)")
    check("90%減少→強い誘導", _cat("AUC", 90, "%", "減少"), "強い誘導相当(AUC80%以上減少)")
    check("60%減少→中等度誘導", _cat("AUC", 60, "%", "減少"), "中等度誘導相当(AUC50〜80%減少)")
    check("30%減少→弱い誘導", _cat("AUC", 30, "%", "減少"), "弱い誘導相当(AUC20〜50%減少)")
    check("0.04倍低下→強い誘導(96%減)", _cat("AUC", 0.04, "倍", "低下"), "強い誘導相当(AUC80%以上減少)")
    check("非AUC(Cmax)→None", _cat("Cmax", 5, "倍", "上昇"), None)


def test_normalize_metric():
    print("\n[pk_numbers._normalize_metric] 指標の表記ゆれ正規化")
    check("AUC下付き(0-∞)→AUC", pk._normalize_metric("AUC0-∞"), "AUC")
    check("AUCτ→AUC", pk._normalize_metric("AUCτ"), "AUC")
    check("t1/2→半減期", pk._normalize_metric("t 1 / 2"), "半減期(t1/2)")
    check("血清中濃度はそのまま", pk._normalize_metric("血清中濃度"), "血清中濃度")


def test_make_change_wide_range():
    print("\n[pk_numbers._make_change] 範囲が広すぎる場合の低信頼フラグ")
    check("1〜10倍は広すぎ(比10)→wide", pk._make_change("AUC", "1", "倍", "上昇", "10")["wide"], True)
    check("1.5〜2.2倍は通常(比1.47)→非wide", pk._make_change("AUC", "1.5", "倍", "上昇", "2.2")["wide"], False)


def test_mentions_and_filter():
    print("\n[pk_numbers._mentions / extract_all] 相手剤フィルタ")
    check("空白を無視して一致", pk._mentions("本 剤は アムロジピン と併用", {"アムロジピン"}), True)
    check("無関係名は不一致", pk._mentions("本剤はワルファリンと併用", {"アムロジピン"}), False)
    text = "クラリスロマイシン併用でAUCが3倍上昇した。リファンピシン併用でAUCが5倍低下した。"
    check("フィルタなしは両文", len(pk.extract_all(text)), 2)
    filt = pk.extract_all(text, partner_names={"クラリスロマイシン"})
    check("相手剤フィルタで1文に絞る", len(filt), 1)
    check("残るのはクラリスの3倍上昇",
          [(c["metric"], c["value_label"], c["direction"]) for c in filt[0]["changes"]],
          [("AUC", "3倍", "上昇")])


# ---------- openfda_lookup ----------

def test_ror_normal():
    print("\n[openfda_lookup._ror] 報告オッズ比の計算")
    r = fda._ror(100, 1000, 1000, 1_000_000)
    check("ROR値", r["ror"], 123.22)
    check("期待併用報告数", r["expected"], 1.0)
    check("実測/期待比", r["obs_exp"], 100.0)
    check_true("CIが値をまたぐ(low<ror<high)", r["ci_low"] < r["ror"] < r["ci_high"])


def test_ror_zerocell_and_invalid():
    print("\n[openfda_lookup._ror] ゼロセル補正・算出不能")
    z = fda._ror(50, 50, 1000, 1_000_000)  # b=solo_a-co=0 → Haldane補正
    check_true("ゼロセルでもNoneでなく算出される", z is not None and "ror" in z)
    check("co=0→算出不能(None)", fda._ror(0, 1000, 1000, 1_000_000), None)
    check("co>solo_a→b<0でNone", fda._ror(100, 50, 1000, 1_000_000), None)
    check("d<=0→None", fda._ror(100, 1000, 1000, 1500), None)


def test_guess_candidates():
    print("\n[openfda_lookup._guess_candidates] 英語名の末尾塩・水和物を剥がす")
    check("末尾Calciumを剥がした中核名を先頭に",
          fda._guess_candidates("Rosuvastatin Calcium"), ["rosuvastatin", "rosuvastatin calcium"])
    check("先頭sodiumは名前の一部→剥がさない",
          fda._guess_candidates("sodium bicarbonate"), ["sodium bicarbonate"])
    check("塩なしはそのまま1件", fda._guess_candidates("amlodipine"), ["amlodipine"])
    check("空文字→空リスト", fda._guess_candidates(""), [])


def test_escape():
    print("\n[openfda_lookup._escape] クエリ用URLエンコード")
    check("スペースを%20に", fda._escape("cefcapene pivoxil"), "cefcapene%20pivoxil")


# ---------- parse_tenpu ----------

def test_parse_new_format():
    print("\n[parse_tenpu.parse] 新書式（番号付き見出し）")
    text = (
        "9. 高齢者への投与\n本剤は慎重に投与する。\n"
        "10. 相互作用\n本剤はCYP3A4で代謝される。\n"
        "10.1 併用禁忌（併用しないこと）\nイトラコナゾール（併用により作用が増強）\n"
        "10.2 併用注意（併用に注意すること）\nクラリスロマイシン（AUC上昇）\n"
        "11. 副作用\n主な副作用は頭痛。\n"
        "16.7 薬物相互作用\nクラリスロマイシン併用でAUCが2.5倍上昇した。\n"
        "17. 臨床成績\n略。\n"
    )
    r = pt.parse(text)
    check("新書式と判定", r["format"], "new")
    check_true("禁忌ブロックにイトラコナゾール", "イトラコナゾール" in (r["contraindicated_combinations"] or ""))
    check_true("注意ブロックにクラリスロマイシン", "クラリスロマイシン" in (r["caution_combinations"] or ""))
    check_true("16.7にAUC2.5倍", "2.5倍" in (r["pk_interactions"] or ""))
    check_true("次章(副作用)は注意ブロックに混入しない", "主な副作用" not in (r["caution_combinations"] or ""))


def test_parse_old_format():
    print("\n[parse_tenpu.parse] 旧書式（番号なし見出し）へフォールバック")
    text = (
        "効能・効果\n略。\n"
        "相互作用\n本剤は腎排泄。\n"
        "併用注意（併用に注意すること）\nリチウム（血中濃度上昇）\n"
        "副作用\n主な副作用は浮腫。\n"
    )
    r = pt.parse(text)
    check("旧書式と判定", r["format"], "old")
    check_true("注意ブロックにリチウム", "リチウム" in (r["caution_combinations"] or ""))


def test_parse_unknown_format():
    print("\n[parse_tenpu.parse] 相互作用章が無い→unknownで全None")
    r = pt.parse("この文書には相互作用の章がありません。副作用のみ記載。")
    check("unknownと判定", r["format"], "unknown")
    check("禁忌None", r["contraindicated_combinations"], None)
    check("注意None", r["caution_combinations"], None)
    check("16.7None", r["pk_interactions"], None)


def test_parse_real_sample():
    print("\n[parse_tenpu.parse] 実サンプル添付文書での回帰ガード（存在すれば）")
    sample = _ROOT / "pmda_sample" / "クラリスロマイシン_tenpu.txt"
    if not sample.exists():
        print("  [SKIP] サンプル未配置のためスキップ")
        return
    r = pt.parse(sample.read_text(encoding="utf-8"))
    check_true("新旧いずれかの書式として認識", r["format"] in ("new", "old"))
    check_true("相互作用セクションが空でない", bool(r["interactions_section"]))


def main():
    print("=" * 70)
    print("パース系純粋関数の単体テスト（オフライン）")
    print("=" * 70)
    for fn in (
        test_extract_multi_metric, test_extract_paired_form, test_fda_auc_category,
        test_normalize_metric, test_make_change_wide_range, test_mentions_and_filter,
        test_ror_normal, test_ror_zerocell_and_invalid, test_guess_candidates, test_escape,
        test_parse_new_format, test_parse_old_format, test_parse_unknown_format,
        test_parse_real_sample,
    ):
        fn()
    print("\n" + "-" * 70)
    if _failures:
        print(f"[FAIL] {len(_failures)}件の失敗: " + "、".join(_failures))
        return 1
    print("[OK] 全テスト通過")
    return 0


if __name__ == "__main__":
    sys.exit(main())
