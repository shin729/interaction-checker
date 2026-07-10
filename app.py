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

    if request.method == "POST":
        outcome = checker.check_interaction(
            request.form.get("drug_a", ""),
            request.form.get("drug_b", ""),
        )
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


if __name__ == "__main__":
    app.run(debug=True, port=5050)
