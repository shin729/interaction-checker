# -*- coding: utf-8 -*-
"""
薬剤間相互作用チェッカー Webアプリ（Flask）。

このモジュールはFlaskのルーティング（HTTP入出力）だけを担う薄いアダプタ。
「2剤名→結果」の中核ロジックは checker.py（Flask非依存）にあり、判定の詳細は
severity.py / pk_numbers.py / openfda_lookup.py の各docstringを参照。

  - 添付文書(PMDA)とopenFDA(FAERS)併用報告から相互作用の強さ「強・中・弱・記載なし」を判定
  - 根拠を優先順位付きで表示
    ①openFDA死亡報告 → ②openFDA重篤報告 → ③添付文書併用禁忌 → ④添付文書併用注意

CYP分類表（cyp_roles.json）の拡充は今後の課題（project_interaction_checker.md参照）。
"""
from flask import Flask, jsonify, render_template, request

import checker
import drug_index
import pmda_lookup

app = Flask(__name__)


@app.route("/suggest")
def suggest():
    q = request.args.get("q", "").strip()
    if len(q) < 3:
        return jsonify([])
    # まずローカル辞書（即時・打ち間違い許容・英語名入力対応）。よく使う薬はPMDA往復なしで返す。
    # ローカルに無い薬（配合剤ブランド名・辞書未収載の薬）だけPMDA前方一致検索へフォールバック。
    names = drug_index.local_suggest(q)
    if not names:
        try:
            names = pmda_lookup.suggest(q, polite=False)
        except Exception:
            names = []
    return jsonify(names)


@app.route("/", methods=["GET", "POST"])
def index():
    result = None
    error = None
    query_a = query_b = ""

    # POST（フォーム送信）に加え、GET ?a=&b= でも実行できる（マトリクスのセルから
    # 特定ペアの詳細へ直接リンクするため。共有リンクにも使える）。
    if request.method == "POST":
        raw_a = request.form.get("drug_a", "")
        raw_b = request.form.get("drug_b", "")
        run = True
    else:
        raw_a = request.args.get("a", "")
        raw_b = request.args.get("b", "")
        run = bool(raw_a and raw_b)

    if run:
        outcome = checker.check_interaction(raw_a, raw_b)
        result = outcome["result"]
        error = outcome["error"]
        query_a = outcome["query_a"]
        query_b = outcome["query_b"]

    return render_template(
        "index.html",
        result=result,
        error=error,
        query_a=query_a,
        query_b=query_b,
    )


_LEVEL_ORDER = {"強": 0, "中": 1, "弱": 2}


@app.route("/matrix", methods=["GET", "POST"])
def matrix():
    data = None
    if request.method == "POST":
        data = checker.check_matrix(request.form.getlist("drugs"))
        if data and not data.get("error"):
            # テンプレート描画用の補助データ（プレゼン層なのでルート側で組み立てる）:
            #  grid_map   : (行,列)→セル の索引（三角グリッド描画用）
            #  notable    : 強/中の注意ペア一覧（重要度順のサマリ）
            #  need_pairs : openFDA遅延取得の対象（記載なしセル）をJSへ渡す
            data["grid_map"] = {f'{c["a_idx"]}-{c["b_idx"]}': c for c in data["cells"]}
            data["notable"] = sorted(
                (c for c in data["cells"] if c["level"] in _LEVEL_ORDER),
                key=lambda c: _LEVEL_ORDER[c["level"]])
            data["need_pairs"] = [
                {"a": c["a"], "b": c["b"], "a_idx": c["a_idx"], "b_idx": c["b_idx"]}
                for c in data["cells"] if c["needs_openfda"]]
    return render_template("matrix.html", data=data)


@app.route("/matrix/openfda", methods=["POST"])
def matrix_openfda():
    """マトリクスの記載なしセルについて、openFDA弱シグナルを遅延取得する（JSONで返す）。"""
    pairs = request.get_json(silent=True) or []
    try:
        return jsonify(checker.matrix_openfda_signals(pairs))
    except Exception:
        return jsonify([])


if __name__ == "__main__":
    app.run(debug=True, port=5050)
