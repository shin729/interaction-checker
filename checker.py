# -*- coding: utf-8 -*-
"""
相互作用チェックの中核ロジック（Flask非依存のサービス層）。

Webルート(app.py)から「2剤の入力→結果dict」の組み立てを分離したモジュール。
Flaskのrequest/responseに一切依存しないため、以下が可能になる:
  - 単体テスト（Flaskを起動せずに検証できる。validate.pyもここを直接叩く）
  - 再利用（将来のJSON API・CLIから同じロジックを呼べる）

方針・各パートの詳細な「なぜ」は severity.py / pk_numbers.py / openfda_lookup.py の
docstringを参照。ここは主に「取得済みデータを優先順位付きで組み立てる」層。
"""
import concurrent.futures

import alternatives
import drug_index
import interaction_predict
import mechanism
import openfda_lookup
import pk_numbers
import pmda_lookup
import severity
import text_format


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


_MAGNITUDE_RANK = {"強い": 3, "中等度": 2, "弱い": 1}


def _magnitude(pk_changes):
    """pk_changes（AUC等の数値変化）から、FDA程度区分が最も強いものを1つ選んで返す。

    判定バッジ（強/中/弱/記載なし）は「どこに記載されたか（禁忌欄/注意欄）」を写すだけで
    「どれくらいの程度か（AUC何倍か）」を反映していない。これを補うため、取得済みの
    AUC変化量（pk_numbers が各変化に付与した `fda` 区分）を集約し、相互作用の大きさの
    目安として判定の横に併記する。判定ロジックには一切影響させない、表示専用の指標。
    取得済みデータの集計のみ＝追加のネットワーク取得なし＝速度影響なし。

    戻り値: {"tier","kind","value_label","direction","source"} または None（AUC数値なし）
    """
    best = None
    for item in pk_changes or []:
        for c in item.get("changes", []):
            fda = c.get("fda")
            if not fda or c.get("metric") != "AUC":
                continue
            tier = next((t for t in _MAGNITUDE_RANK if fda.startswith(t)), None)
            if tier is None:
                continue
            rank = _MAGNITUDE_RANK[tier]
            if best is None or rank > best["rank"]:
                best = {
                    "tier": tier,                                   # 強い / 中等度 / 弱い
                    "rank": rank,
                    "kind": "阻害" if "阻害" in fda else ("誘導" if "誘導" in fda else ""),
                    "value_label": c.get("value_label"),            # 例: 5.1倍
                    "direction": c.get("direction"),                # 例: 上昇
                    "source": item.get("source"),
                }
    return best


def check_interaction(raw_a, raw_b):
    """2剤の生入力から結果を組み立てる中核関数（Flask非依存）。

    戻り値: {"result": dict|None, "error": str|None, "query_a": str, "query_b": str}
      - result   : 相互作用の結果dict（成功時のみ。テンプレートがそのまま描画する）
      - error    : ユーザー向けエラーメッセージ（失敗時のみ）
      - query_a/b: 英語名を日本語一般名へ正規化した後の検索クエリ（フォーム再表示用）
    """
    # 英語名で入力された場合は日本語一般名へ正規化する（amlodipine→アムロジピン）。
    # PMDA検索は日本語名で引くため。完全一致のみ変換し、それ以外はそのまま。
    query_a = drug_index.resolve_input(raw_a or "")
    query_b = drug_index.resolve_input(raw_b or "")

    if not query_a or not query_b:
        return {"result": None, "error": "薬剤名を2つとも入力してください。",
                "query_a": query_a, "query_b": query_b}

    try:
        # 2剤の添付文書取得は互いに独立。逐次だと1剤約6秒×2かかるため並列化する。
        # 対話リクエストでは polite=False にして取得後の待機(sleep)を省く。
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
            fut_a = ex.submit(pmda_lookup.lookup, query_a, use_cache=True, polite=False)
            fut_b = ex.submit(pmda_lookup.lookup, query_b, use_cache=True, polite=False)
            pmda_a = fut_a.result()
            pmda_b = fut_b.result()
    except Exception as e:
        return {"result": None, "error": f"添付文書の検索中にエラーが発生しました: {e}",
                "query_a": query_a, "query_b": query_b}

    if not pmda_a["found"]:
        return {"result": None,
                "error": f"「{query_a}」の添付文書が見つかりませんでした。薬剤名（一般名）を確認してください。",
                "query_a": query_a, "query_b": query_b}
    if not pmda_b["found"]:
        return {"result": None,
                "error": f"「{query_b}」の添付文書が見つかりませんでした。薬剤名（一般名）を確認してください。",
                "query_a": query_a, "query_b": query_b}

    try:
        fda_stats = openfda_lookup.lookup_pair(
            query_a, pmda_a.get("matched_name") or query_a,
            query_b, pmda_b.get("matched_name") or query_b,
            guess_a=pmda_a.get("english_name_guess"),
            guess_b=pmda_b.get("english_name_guess"),
            polite=False,
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

    # CYP/P-糖蛋白の役割からの機序ベース予測（ローカル表のみ＝ネットワーク不要）。
    # 添付文書・openFDAに記載が無くても機序的に注意すべきペアを補う。判定(verdict)には
    # 一切影響させず、文書化された根拠とは別枠の「参考予測」として表示する。
    predictions = interaction_predict.predict(pmda_a["matched_name"], pmda_b["matched_name"])

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
        "predictions": predictions,
        "magnitude": _magnitude(pk_changes),
    }
    return {"result": result, "error": None, "query_a": query_a, "query_b": query_b}


# ---- 多剤併用マトリクス（添付文書ベースを即時表示、openFDAは遅延取得）----

# 一度にチェックできる薬剤数の上限。ペア数は N(N-1)/2 で増えるため、openFDA遅延取得の
# 最悪時間を現実的に抑える目的で上限を設ける（10剤＝45ペア）。
MATRIX_MAX_DRUGS = 10


def check_matrix(raw_names):
    """複数薬剤の総当たり相互作用マトリクスを組み立てる。

    ペアの重い部分（openFDA併用報告の通信＝1ペア約2.5秒）は含めず、添付文書ベースの判定
    （強/中/記載なし）だけを即座に返す。openFDAでしか出ない「弱」は、記載なしセルに対して
    matrix_openfda_signals() で遅延取得する（needs_openfda=True のセルが対象）。

    戻り値: {
      "drugs":     [{"input","name","found"}...],   # 添付文書が見つかった薬（重複排除後）
      "not_found": [str...],                          # 見つからなかった入力
      "cells":     [{"a_idx","b_idx","a","b","a_name","b_name","level","reason","needs_openfda"}...],
      "n":         見つかった薬の数,
      "error":     str|None,
    }
    """
    # 英語名正規化・空除去・重複排除（入力順を保つ）
    names, seen = [], set()
    for raw in raw_names or []:
        q = drug_index.resolve_input((raw or "").strip())
        if q and q not in seen:
            seen.add(q)
            names.append(q)

    if len(names) < 2:
        return {"drugs": [], "not_found": [], "cells": [], "n": 0,
                "error": "薬剤名を2つ以上入力してください。"}
    if len(names) > MATRIX_MAX_DRUGS:
        return {"drugs": [], "not_found": [], "cells": [], "n": 0,
                "error": f"一度にチェックできるのは{MATRIX_MAX_DRUGS}剤までです"
                         f"（{len(names)}剤が入力されました）。数を減らしてください。"}

    # 添付文書は薬ごとに1回だけ取得（ペアごとではない）。互いに独立なので並列化する。
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(len(names), 8)) as ex:
        infos = list(ex.map(
            lambda q: pmda_lookup.lookup(q, use_cache=True, polite=False), names))

    found, found_infos, not_found = [], [], []
    for q, info in zip(names, infos):
        if info.get("found"):
            found.append({"input": q, "name": info.get("matched_name") or q, "found": True})
            found_infos.append(info)
        else:
            not_found.append(q)

    cells = []
    for a in range(len(found)):
        for b in range(a + 1, len(found)):
            qa, qb = found[a]["input"], found[b]["input"]
            # openFDAは遅延取得のため、ここでは fda_stats=None で添付文書のみから判定する。
            verdict = severity.classify(found_infos[a], found_infos[b], qa, qb, None)
            level = verdict["level"]
            cells.append({
                "a_idx": a, "b_idx": b, "a": qa, "b": qb,
                "a_name": found[a]["name"], "b_name": found[b]["name"],
                "level": level, "reason": verdict["reason"],
                # openFDAで「弱」に上がりうるのは添付文書に記載が無い（記載なし）ペアだけ。
                # 強/中は既に文書化済みなので遅延取得の対象にしない（無駄な通信を省く）。
                "needs_openfda": level == "記載なし",
            })

    return {"drugs": found, "not_found": not_found, "cells": cells,
            "n": len(found), "error": None}


def pair_openfda_signal(raw_a, raw_b):
    """1ペアのopenFDA弱シグナルを取得する（マトリクスの記載なしセルを遅延で埋める用）。

    添付文書ベースで既に強/中と判定済みのペアには呼ばない前提。openFDA併用報告に
    相互作用シグナル（DRUG INTERACTION等）があれば「弱」、無ければ「記載なし」を返す。

    戻り値: {"level": "弱"|"記載なし", "ror": float|None, "available": bool}
    """
    query_a = drug_index.resolve_input((raw_a or "").strip())
    query_b = drug_index.resolve_input((raw_b or "").strip())
    pa = pmda_lookup.lookup(query_a, use_cache=True, polite=False)
    pb = pmda_lookup.lookup(query_b, use_cache=True, polite=False)
    if not pa.get("found") or not pb.get("found"):
        return {"level": "記載なし", "ror": None, "available": False}
    try:
        fda_stats = openfda_lookup.lookup_pair(
            query_a, pa.get("matched_name") or query_a,
            query_b, pb.get("matched_name") or query_b,
            guess_a=pa.get("english_name_guess"), guess_b=pb.get("english_name_guess"),
            polite=False)
    except Exception:
        fda_stats = None
    signal = severity._openfda_signal(fda_stats)
    if signal:
        return {"level": "弱", "ror": signal.get("ror"), "available": True}
    return {"level": "記載なし",
            "ror": fda_stats.get("ror") if fda_stats else None,
            "available": fda_stats is not None}


def matrix_openfda_signals(pairs):
    """マトリクスの複数ペアについてopenFDA弱シグナルをまとめて取得する。

    ペアごとのopenFDA通信は互いに独立なので、サーバ側でスレッド並列化する（gunicornの
    ワーカー数に依存せず1リクエスト内で並列取得できる）。openFDAはネットワークI/O待ちが
    主体でGILが解放されるため、スレッド並列がそのまま効く。

    引数 pairs: [{"a","b","a_idx","b_idx"}...]（check_matrix の needs_openfda セル）
    戻り値:     [{"a_idx","b_idx","level","ror","available"}...]
    """
    pairs = (pairs or [])[:MATRIX_MAX_DRUGS * (MATRIX_MAX_DRUGS - 1) // 2]  # 上限ペア数で安全弁
    if not pairs:
        return []

    def _one(p):
        sig = pair_openfda_signal(p.get("a", ""), p.get("b", ""))
        return {"a_idx": p.get("a_idx"), "b_idx": p.get("b_idx"), **sig}

    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        return list(ex.map(_one, pairs))
