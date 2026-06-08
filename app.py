# -*- coding: utf-8 -*-
"""
薬剤間相互作用チェッカー Webアプリ（Flask）

2つの薬剤名を入力すると:
  - 添付文書(PMDA)とopenFDA(FAERS)併用報告から相互作用の強さ「強・中・弱・記載なし」を判定
  - 根拠を優先順位付きで表示
    ①openFDA死亡報告 → ②openFDA重篤報告 → ③添付文書併用禁忌 → ④添付文書併用注意

PubMed/CYP分類表の統合は今後の拡張（project_interaction_checker.md参照）。
"""
from flask import Flask, jsonify, render_template, request

import mechanism
import openfda_lookup
import pk_numbers
import pmda_lookup
import severity
import text_format

app = Flask(__name__)


def _build_pk_changes(pmda_a, pmda_b, query_a, query_b):
    """両剤の添付文書から、AUC/Cmax/クリアランス等の定量的な変化（倍率・%）の言及を集める。

    添付文書では「具体的にどれくらい変化するか」の記載が少なく、書籍等でも
    横断的に探しにくいため、本ツールではこれを最優先で見やすく提示する。
    """
    items = []
    for label_name, info in (
        (pmda_a.get("matched_name") or query_a, pmda_a),
        (pmda_b.get("matched_name") or query_b, pmda_b),
    ):
        for entry in pk_numbers.extract_all(
            info.get("pk_interactions"),
            info.get("caution_combinations"),
            info.get("contraindicated_combinations"),
        ):
            items.append({
                "source": label_name,
                "changes": entry["changes"],
                "sentence": entry["sentence"],
            })
    return items


def _build_evidence(pmda_a, pmda_b, query_a, query_b, fda_stats):
    """根拠を優先順位順に並べたリストを返す（①死亡→②重篤→③併用禁忌→④併用注意）"""
    items = []

    if fda_stats and fda_stats["co_reports_death"] > 0:
        items.append({
            "priority": 1,
            "label": "openFDA(FAERS) 死亡転帰の併用報告",
            "summary": f"{fda_stats['drug_a']}と{fda_stats['drug_b']}の併用報告"
                       f"{fda_stats['co_reports_total']:,}件中、死亡転帰"
                       f"{fda_stats['co_reports_death']:,}件（{fda_stats['death_ratio']:.1%}）",
            "detail": None,
            "numeric": True,
        })

    if fda_stats and fda_stats["co_reports_serious"] > 0:
        items.append({
            "priority": 2,
            "label": "openFDA(FAERS) 重篤(serious)併用報告",
            "summary": f"併用報告{fda_stats['co_reports_total']:,}件中、重篤報告"
                       f"{fda_stats['co_reports_serious']:,}件（{fda_stats['serious_ratio']:.1%}）。"
                       f"頻出有害事象: "
                       + "、".join(f"{t}({c:,}件)" for t, c in fda_stats["top_reactions"][:5]),
            "detail": None,
            "numeric": True,
        })

    for label_name, info, other_query, other_name in (
        (pmda_a.get("matched_name") or query_a, pmda_a, query_b, pmda_b.get("matched_name") or query_b),
        (pmda_b.get("matched_name") or query_b, pmda_b, query_a, pmda_a.get("matched_name") or query_a),
    ):
        other_names = severity._expand_names(other_query, other_name)
        contraindicated = info.get("contraindicated_combinations")
        if contraindicated and severity._text_mentions(contraindicated, *other_names):
            items.append({
                "priority": 3,
                "label": f"添付文書（{label_name}）の「併用禁忌」",
                "summary": f"{other_name}が併用禁忌として記載されています。",
                "detail": text_format.reflow(contraindicated),
                "numeric": False,
            })

        caution = info.get("caution_combinations")
        if caution and severity._text_mentions(caution, *other_names):
            items.append({
                "priority": 4,
                "label": f"添付文書（{label_name}）の「併用注意」",
                "summary": f"{other_name}が併用注意として記載されています。",
                "detail": text_format.reflow(caution),
                "numeric": False,
            })

    items.sort(key=lambda x: x["priority"])
    return items


@app.route("/suggest")
def suggest():
    q = request.args.get("q", "").strip()
    if len(q) < 3:
        return jsonify([])
    try:
        names = pmda_lookup.suggest(q)
    except Exception:
        names = []
    return jsonify(names)


@app.route("/", methods=["GET", "POST"])
def index():
    result = None
    error = None
    query_a = query_b = ""

    if request.method == "POST":
        query_a = request.form.get("drug_a", "").strip()
        query_b = request.form.get("drug_b", "").strip()

        if not query_a or not query_b:
            error = "薬剤名を2つとも入力してください。"
        else:
            try:
                pmda_a = pmda_lookup.lookup(query_a)
                pmda_b = pmda_lookup.lookup(query_b)
            except Exception as e:
                error = f"添付文書の検索中にエラーが発生しました: {e}"
                pmda_a = pmda_b = None

            if pmda_a is not None:
                if not pmda_a["found"]:
                    error = f"「{query_a}」の添付文書が見つかりませんでした。薬剤名（一般名）を確認してください。"
                elif not pmda_b["found"]:
                    error = f"「{query_b}」の添付文書が見つかりませんでした。薬剤名（一般名）を確認してください。"
                else:
                    try:
                        fda_stats = openfda_lookup.lookup_pair(
                            query_a, pmda_a.get("matched_name") or query_a,
                            query_b, pmda_b.get("matched_name") or query_b,
                        )
                    except Exception:
                        fda_stats = None

                    verdict = severity.classify(pmda_a, pmda_b, query_a, query_b, fda_stats)
                    mech_text = verdict.pop("mechanism_text", None)
                    mech_names = verdict.pop("mechanism_names", ())
                    verdict["mechanisms"] = (
                        mechanism.extract_near(mech_text, *mech_names) if mech_text else []
                    )
                    evidence = _build_evidence(pmda_a, pmda_b, query_a, query_b, fda_stats)
                    pk_changes = _build_pk_changes(pmda_a, pmda_b, query_a, query_b)

                    result = {
                        "name_a": pmda_a["matched_name"],
                        "name_b": pmda_b["matched_name"],
                        "verdict": verdict,
                        "evidence": evidence,
                        "pk_changes": pk_changes,
                        "fda_stats": fda_stats,
                        "fda_name_missing": fda_stats is None,
                    }

    return render_template(
        "index.html",
        result=result,
        error=error,
        query_a=query_a,
        query_b=query_b,
    )


if __name__ == "__main__":
    app.run(debug=True, port=5050)
