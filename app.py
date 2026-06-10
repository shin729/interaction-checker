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

import alternatives
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

    各剤の添文に並ぶ「他の薬剤との相互作用」を無差別に拾うと、いま調べている
    ペアと無関係な第三の薬剤の数値まで出てしまうため、相手剤名（中核名・配合成分名を
    含む）を含む文だけに絞り込む。

    添付文書に加え、インタビューフォーム(IF)から抽出済みの数値(if_pk_items)も
    相手剤フィルタして併合する。IFは添付文書より詳細な相互作用試験データを載せるため
    数値の網羅が増える。同じ数値（指標+値+方向が一致）は添付文書を優先して重複排除する。
    """
    name_a = pmda_a.get("matched_name") or query_a
    name_b = pmda_b.get("matched_name") or query_b
    items = []
    for label_name, info, other_query, other_name in (
        (name_a, pmda_a, query_b, name_b),
        (name_b, pmda_b, query_a, name_a),
    ):
        partner_names = set()
        for n in severity._expand_names(other_query, other_name):
            partner_names.add(n)
            partner_names.add(severity._name_core(n))

        seen_sigs = set()

        def _norm_dir(d):
            # 添付文書「増加/減少」とIF「上昇/低下」は同義。重複排除のため粗い向きに正規化する
            if d in ("上昇", "増加", "増強", "延長"):
                return "up"
            if d in ("低下", "減少", "減弱", "短縮"):
                return "down"
            return "flat"

        def _sig(changes):
            return frozenset((c["metric"], c["value_label"], _norm_dir(c["direction"])) for c in changes)

        # ① 添付文書（16.7・併用注意・併用禁忌）
        for entry in pk_numbers.extract_all(
            info.get("pk_interactions"),
            info.get("caution_combinations"),
            info.get("contraindicated_combinations"),
            partner_names=partner_names,
        ):
            seen_sigs.add(_sig(entry["changes"]))
            items.append({
                "source": label_name, "source_type": "添付文書",
                "changes": entry["changes"], "sentence": entry["sentence"],
            })

        # ② インタビューフォーム（相手剤フィルタ＋添付文書と重複する数値は除外）
        for entry in (info.get("if_pk_items") or []):
            if not pk_numbers._mentions(entry["sentence"], partner_names):
                continue
            sig = _sig(entry["changes"])
            if sig in seen_sigs:
                continue
            seen_sigs.add(sig)
            items.append({
                "source": label_name, "source_type": "インタビューフォーム",
                "changes": entry["changes"], "sentence": entry["sentence"],
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

    if fda_stats and fda_stats.get("ror") is not None:
        lo, hi = fda_stats["ror_ci_low"], fda_stats["ror_ci_high"]
        if lo > 1:
            interp = "各剤単独の報告頻度から期待される水準より有意に多く併用報告されています（信頼区間の下限>1）。"
        else:
            interp = "信頼区間の下限が1を下回るため、統計的に有意な偏りとは言えません。"
        items.append({
            "priority": 2,
            "label": "openFDA(FAERS) 併用報告の不均衡（ROR・報告オッズ比）",
            "summary": f"併用報告{fda_stats['co_reports_total']:,}件。"
                       f"独立を仮定した期待併用報告数{fda_stats['expected_co']:,.1f}件に対し、"
                       f"ROR {fda_stats['ror']:.1f}（95%CI {lo:.1f}–{hi:.1f}）。{interp} "
                       f"※RORは相互作用だけでなく「併用処方の多さ」も反映するため、"
                       f"高ROR＝相互作用が強い、と短絡できない点にご注意ください。",
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
        mentioned_in_block = False
        if caution and severity._text_mentions(caution, *other_names):
            mentioned_in_block = True
            items.append({
                "priority": 4,
                "label": f"添付文書（{label_name}）の「併用注意」",
                "summary": f"{other_name}が併用注意として記載されています。",
                "detail": text_format.reflow(caution),
                "numeric": False,
            })

        # 16.7 薬物相互作用セクションの記載。併用禁忌/注意ブロックに既に相手剤が
        # 出ている場合は重複になるので追加しない（旧書式PDFで抽出順が乱れ、注意
        # ブロックに入らず16.7だけに相手剤が出るケースを根拠欄に反映するのが主目的）。
        pk_text = info.get("pk_interactions")
        if not mentioned_in_block and severity._pk_interaction_with(pk_text, *other_names):
            items.append({
                "priority": 4,
                "label": f"添付文書（{label_name}）の「薬物相互作用(16.7)」",
                "summary": f"{other_name}との薬物動態学的相互作用が記載されています"
                           f"（具体的な変化量は上の「薬物動態への影響」を参照）。",
                "detail": text_format.reflow(pk_text),
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

                    # 同じ薬効系統の別の薬（調べ直し候補）。相手剤は除外する。
                    alt_a = alternatives.find_alternatives(pmda_a["matched_name"], exclude=(pmda_b["matched_name"],))
                    alt_b = alternatives.find_alternatives(pmda_b["matched_name"], exclude=(pmda_a["matched_name"],))

                    result = {
                        "name_a": pmda_a["matched_name"],
                        "name_b": pmda_b["matched_name"],
                        "verdict": verdict,
                        "evidence": evidence,
                        "pk_changes": pk_changes,
                        "fda_stats": fda_stats,
                        "fda_name_missing": fda_stats is None,
                        "alt_a": alt_a,
                        "alt_b": alt_b,
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
