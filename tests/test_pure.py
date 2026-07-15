# -*- coding: utf-8 -*-
"""
純粋関数のオフライン単体テスト（ネットワーク・キャッシュ不要、stdlibのみ）。

validate.py は実ペアをパイプライン全体に通す統合テストで、ネットワーク/キャッシュに
依存し実行が重い。一方こちらは、名寄せ・入力正規化・程度集約・機序予測といった
決定論的な純粋関数を、合成入力で高速(<1秒)に固定する回帰ネット。パース/判定ロジックを
リファクタしても即座に壊れを検出できるようにするのが目的。

  python tests/test_pure.py      # 失敗が1件でもあれば exit 1

pytest等の追加依存は入れず、requirements.txtの軽量方針を保つ。
"""
import sys
from pathlib import Path

# プロジェクトルート（このファイルの親の親）をインポートパスに追加
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

try:
    sys.stdout.reconfigure(encoding="utf-8")  # Windowsコンソール(cp932)対策
except Exception:
    pass

import checker
import drug_index
import interaction_predict
import severity

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


def test_name_core():
    print("\n[severity._name_core] 塩・水和物・エステルを落とした中核名")
    check("ベシル酸塩を除去", severity._name_core("アムロジピンベシル酸塩"), "アムロジピン")
    check("エステル+塩+水和物を除去", severity._name_core("セフカペン ピボキシル塩酸塩水和物"), "セフカペン")
    check("修飾語なしはそのまま", severity._name_core("ワルファリン"), "ワルファリン")


def test_split_combo():
    print("\n[severity._split_combo] 配合剤の複合名を成分へ分解")
    check("配合剤を成分に分解",
          severity._split_combo("テルミサルタン・アムロジピンベシル酸塩配合剤"),
          ["テルミサルタン", "アムロジピンベシル酸塩"])
    check("単剤はそのまま1要素", severity._split_combo("アムロジピン"), ["アムロジピン"])


def test_expand_names():
    print("\n[severity._expand_names] クエリ+マッチ名+成分名を重複排除")
    got = severity._expand_names("アムロジピン", "アムロジピンベシル酸塩")
    check_true("tupleで返る", isinstance(got, tuple))
    check("重複排除された集合", set(got), {"アムロジピン", "アムロジピンベシル酸塩"})


def test_text_mentions():
    print("\n[severity._text_mentions] 相手剤名（中核名含む）の言及照合")
    check("塩付き名でも中核名で一致",
          severity._text_mentions("本剤はアムロジピンとの併用で…", "アムロジピンベシル酸塩"), True)
    check("無関係な薬剤名は不一致",
          severity._text_mentions("本剤はワルファリンとの併用で…", "アムロジピン"), False)


def test_resolve_input():
    print("\n[drug_index.resolve_input] 英語名の完全一致のみ日本語一般名へ正規化")
    check("英語名→日本語一般名", drug_index.resolve_input("amlodipine"), "アムロジピン")
    check("日本語はそのまま", drug_index.resolve_input("アムロジピン"), "アムロジピン")
    check("未知の語は書き換えない（医療安全）", drug_index.resolve_input("zzunknownzz"), "zzunknownzz")


def test_local_suggest():
    print("\n[drug_index.local_suggest] ローカル辞書サジェスト")
    got = drug_index.local_suggest("アムロ")
    check_true("前方一致でアムロジピンを含む", "アムロジピン" in got)
    check("2文字未満は空", drug_index.local_suggest("ア"), [])


def test_predict():
    print("\n[interaction_predict.predict] CYP役割からの機序予測")
    preds = interaction_predict.predict("トリアゾラム", "イトラコナゾール")
    check_true("予測が1件以上返る", len(preds) >= 1)
    if preds:
        p = preds[0]
        check("victimは基質側(トリアゾラム)", p["victim"], "トリアゾラム")
        check("強いCYP3A4阻害を予測", (p["level"], p["kind"], p["enzyme"]), ("強い", "阻害", "CYP3A4"))
    check("役割表に無いペアは空", interaction_predict.predict("zz1", "zz2"), [])

    # 拡充データの重要ペア（コルヒチン+マクロライドは潜在的に致死的な古典的相互作用）
    cp = {(p["level"], p["enzyme"]) for p in interaction_predict.predict("コルヒチン", "クラリスロマイシン")}
    check_true("コルヒチン×クラリス: 強いCYP3A4阻害を予測", ("強い", "CYP3A4") in cp)
    check_true("コルヒチン×クラリス: P-gp阻害も予測", any(e == "P-gp" for _, e in cp))
    tk = interaction_predict.predict("チカグレロル", "イトラコナゾール")
    check_true("チカグレロル×イトラコナゾール: 強いCYP3A4阻害",
               any(p["level"] == "強い" and p["enzyme"] == "CYP3A4" for p in tk))
    dg = interaction_predict.predict("ジゴキシン", "ジルチアゼム")
    check_true("ジゴキシン×ジルチアゼム: P-gp阻害を予測（ジルチアゼムP-gp追加）",
               any(p["enzyme"] == "P-gp" for p in dg))


def test_magnitude():
    print("\n[checker._magnitude] AUC変化から最大の程度区分を1つ選ぶ（表示専用）")
    mid = [{"source": "S", "changes": [
        {"metric": "AUC", "fda": "中等度阻害相当(AUC2〜5倍)", "value_label": "3倍", "direction": "上昇"}]}]
    strong = [{"source": "S", "changes": [
        {"metric": "AUC", "fda": "弱い阻害相当(AUC1.25〜2倍)", "value_label": "1.5倍", "direction": "上昇"},
        {"metric": "AUC", "fda": "強い阻害相当(AUC5倍以上)", "value_label": "6倍", "direction": "上昇"}]}]
    noauc = [{"source": "S", "changes": [
        {"metric": "Cmax", "fda": "中等度阻害相当(AUC2〜5倍)", "value_label": "3倍", "direction": "上昇"}]}]

    m = checker._magnitude(mid)
    check("中等度阻害を拾う", (m["tier"], m["kind"]), ("中等度", "阻害"))
    check("複数あれば最強(強い)を選ぶ", checker._magnitude(strong)["tier"], "強い")
    check("AUC以外の指標は対象外→None", checker._magnitude(noauc), None)
    check("変化なし→None", checker._magnitude([]), None)
    check("Noneを渡してもNone", checker._magnitude(None), None)


def test_cyp_roles_integrity():
    print("\n[cyp_roles.json / drug_name_map.json] 知識テーブルの整合性")
    import json
    root = Path(__file__).resolve().parent.parent
    roles = json.loads((root / "cyp_roles.json").read_text(encoding="utf-8"))["roles"]
    nmap = json.loads((root / "drug_name_map.json").read_text(encoding="utf-8"))
    check_true("cyp_rolesは50薬以上", len(roles) >= 50)
    # プレフィックス衝突: あるキーが別キーの接頭辞だと _roles_for(前方一致)で誤シャドウしうる
    keys = list(roles.keys())
    collisions = [(a, b) for a in keys for b in keys if a != b and b.startswith(a)]
    check("キーのプレフィックス衝突なし", collisions, [])
    # 予測できる薬はopenFDA照会もできるべき（食品のグレープフルーツのみ例外）
    missing = [k for k in roles if k not in nmap and k != "グレープフルーツ"]
    check("cyp_roles薬は全てname_mapにある(食品除く)", missing, [])
    # 値の妥当性: substrateは感受性/基質、inhibitor/inducerは強/中/弱のみ
    bad = []
    for name, r in roles.items():
        for v in (r.get("substrate") or {}).values():
            if v not in ("感受性", "基質"):
                bad.append((name, "substrate", v))
        for role in ("inhibitor", "inducer"):
            for v in (r.get(role) or {}).values():
                if v not in ("強", "中", "弱"):
                    bad.append((name, role, v))
    check("役割値は既定の語彙のみ", bad, [])


def main():
    print("=" * 70)
    print("純粋関数の単体テスト（オフライン）")
    print("=" * 70)
    for fn in (test_name_core, test_split_combo, test_expand_names, test_text_mentions,
               test_resolve_input, test_local_suggest, test_predict, test_magnitude,
               test_cyp_roles_integrity):
        fn()
    print("\n" + "-" * 70)
    if _failures:
        print(f"[FAIL] {len(_failures)}件の失敗: " + "、".join(_failures))
        return 1
    print("[OK] 全テスト通過")
    return 0


if __name__ == "__main__":
    sys.exit(main())
