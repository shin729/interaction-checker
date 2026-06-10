# -*- coding: utf-8 -*-
"""
PMDA「医療用医薬品 添付文書等情報検索」から薬剤名を渡して
添付文書PDFを取得し、相互作用関連セクションを構造化して返すモジュール。

verify_pmda_access.py で確立した検索フロー（フォーム全項目POST方式）と
parse_tenpu.py の章節抽出ロジックを、アプリから呼び出せる関数にまとめたもの。
ローカルにJSONキャッシュを持ち、同じ薬剤名への再検索を避ける。
"""
import json
import re
import time
from pathlib import Path

import fitz  # PyMuPDF
import requests
from bs4 import BeautifulSoup

import pk_numbers
from parse_tenpu import parse as parse_interactions

BASE = "https://www.pmda.go.jp"
SEARCH_URL = f"{BASE}/PmdaSearch/iyakuSearch/"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36",
}

CACHE_DIR = Path(__file__).parent / "cache"
CACHE_DIR.mkdir(exist_ok=True)
SUGGEST_CACHE_DIR = CACHE_DIR / "suggest"
SUGGEST_CACHE_DIR.mkdir(exist_ok=True)

_UNSAFE_FILENAME_RE = re.compile(r'[\\/:*?"<>|]')


def _cache_path(name: str) -> Path:
    safe = _UNSAFE_FILENAME_RE.sub("_", name)
    return CACHE_DIR / f"{safe}.json"


def _suggest_cache_path(prefix: str) -> Path:
    safe = _UNSAFE_FILENAME_RE.sub("_", prefix)
    return SUGGEST_CACHE_DIR / f"{safe}.json"


def _submit_search(session: requests.Session, name: str) -> str:
    r = session.get(SEARCH_URL, headers=HEADERS, timeout=30)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "lxml")
    form = soup.find("form", {"id": "iyakuSearchForm"})

    overrides = {
        "nameWord": name,
        "iyakuHowtoNameSearchRadioValue": "1",
        "howtoMatchRadioValue": "2",
    }
    seen, payload = set(), []
    for inp in form.find_all(["input", "select", "textarea"]):
        name_attr = inp.get("name")
        if not name_attr:
            continue
        tag = inp.name
        if tag == "input":
            itype = (inp.get("type") or "text").lower()
            if itype in ("image", "button"):
                continue
            if itype in ("checkbox", "radio"):
                if not inp.has_attr("checked"):
                    continue
                value = inp.get("value", "on")
            else:
                value = inp.get("value", "")
        elif tag == "select":
            opt = inp.find("option", selected=True) or inp.find("option")
            value = opt.get("value", "") if opt is not None else ""
        else:
            value = inp.text or ""

        if name_attr in overrides and name_attr not in seen:
            value = overrides[name_attr]
            seen.add(name_attr)
        payload.append((name_attr, value))
    for k, v in overrides.items():
        if k not in seen:
            payload.append((k, v))
    payload += [("btnA.x", "15"), ("btnA.y", "10")]

    r2 = session.post(SEARCH_URL, data=payload, headers={**HEADERS, "Referer": SEARCH_URL}, timeout=30)
    r2.raise_for_status()
    return r2.text


def _general_list_links(html: str):
    soup = BeautifulSoup(html, "lxml")
    seen, links = set(), []
    for a in soup.select('a[href*="GeneralList"]'):
        href = a.get("href", "")
        text = a.get_text(strip=True)
        if not href or not text:
            continue
        full = requests.compat.urljoin(SEARCH_URL, href)
        if full not in seen:
            seen.add(full)
            links.append((text, full))
    return links


def _detail_pdf_links(session: requests.Session, detail_url: str):
    """詳細ページから添付文書PDFとインタビューフォーム(IF)PDFのURLを返す。

    添付文書は URL に "resultdatasetpdf"、IFは "/interview/" を含むのが目印。
    IFは添付文書より詳細な薬物相互作用試験データ（AUC/Cmaxの数値）を載せるため、
    「程度の具体性」を補う第2のソースとして使う。"""
    r = session.get(detail_url, headers={**HEADERS, "Referer": SEARCH_URL}, timeout=30)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "lxml")
    tenpu, interview = None, None
    for a in soup.find_all("a", href=True):
        full = requests.compat.urljoin(detail_url, a["href"])
        low = full.lower()
        if tenpu is None and "resultdatasetpdf" in low:
            tenpu = full
        if interview is None and "/interview/" in low:
            interview = full
    return tenpu, interview


def _pdf_text(session: requests.Session, pdf_url: str, referer: str) -> str:
    r = session.get(pdf_url, headers={**HEADERS, "Referer": referer}, timeout=60)
    r.raise_for_status()
    doc = fitz.open(stream=r.content, filetype="pdf")
    text = "\n".join(page.get_text() for page in doc)
    doc.close()
    return text


def lookup(drug_name: str, use_cache: bool = True) -> dict:
    """
    薬剤名(一般名/販売名の前方一致)から添付文書を検索し、相互作用情報を返す。

    戻り値:
      {
        "query": 検索した薬剤名,
        "matched_name": 検索でヒットした薬剤の正式名,
        "detail_url": PMDA詳細ページURL,
        "pdf_url": 添付文書PDF URL,
        "format": "new" | "old" | "unknown",
        "interactions_section": str | None,
        "contraindicated_combinations": str | None,
        "caution_combinations": str | None,
        "pk_interactions": str | None,
        "if_pdf_url": str | None,           # インタビューフォームPDFのURL
        "if_pk_items": [ {"sentence", "changes"}, ... ],  # IFから抽出したPK数値文
        "found": bool,
      }
    """
    cache_file = _cache_path(drug_name)
    if use_cache and cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    result = {
        "query": drug_name, "matched_name": None, "detail_url": None, "pdf_url": None,
        "format": "unknown", "interactions_section": None,
        "contraindicated_combinations": None, "caution_combinations": None,
        "pk_interactions": None, "if_pdf_url": None, "if_pk_items": [], "found": False,
    }

    with requests.Session() as session:
        html = _submit_search(session, drug_name)
        links = _general_list_links(html)
        if not links:
            cache_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
            return result

        matched_name, detail_url = links[0]
        pdf_url, if_pdf_url = _detail_pdf_links(session, detail_url)
        result.update(matched_name=matched_name, detail_url=detail_url, pdf_url=pdf_url,
                      if_pdf_url=if_pdf_url, found=True)
        if pdf_url:
            text = _pdf_text(session, pdf_url, detail_url)
            parsed = parse_interactions(text)
            result.update(parsed)
        # インタビューフォームから定量的PK数値を補完抽出する（任意・失敗時はスキップ）。
        # 相手剤に依存しない全PK文を保存し、照合（相手剤フィルタ）はapp側で行う。
        if if_pdf_url:
            try:
                if_text = _pdf_text(session, if_pdf_url, detail_url)
                result["if_pk_items"] = pk_numbers.extract(if_text)
            except Exception:
                result["if_pk_items"] = []

    cache_file.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    time.sleep(1)  # PMDAサーバへの負荷軽減（連続アクセス時のマナー）
    return result


def suggest(prefix: str, limit: int = 8, use_cache: bool = True) -> list:
    """
    入力途中の薬剤名（3文字以上）から、PMDAの前方一致検索で候補名（一般名/販売名）を
    返す。入力補助（オートコンプリート）用。

    PMDAは公式の候補API/サジェスト機能を提供していない。部分一致モード
    （howtoMatchRadioValue=1）も試したが、無関係な薬剤名が紛れ込み結果順序も
    最適でなかった（例:「ロキソ」→「アンブロキソール塩酸塩」が先頭に出る等）ため、
    `lookup()`と同じ前方一致モード（howtoMatchRadioValue=2）を流用する方が
    タイプ中の候補として明らかに適切だった。1回の検索に約2秒かかるため、
    プレフィックス文字列単位でローカルキャッシュし（cache/suggest/）、
    同じ語の再検索を避ける。
    """
    prefix = prefix.strip()
    if len(prefix) < 3:
        return []

    cache_file = _suggest_cache_path(prefix)
    if use_cache and cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8"))

    with requests.Session() as session:
        html = _submit_search(session, prefix)
        links = _general_list_links(html)

    seen, names = set(), []
    for text, _ in links:
        if text not in seen:
            seen.add(text)
            names.append(text)
        if len(names) >= limit:
            break

    cache_file.write_text(json.dumps(names, ensure_ascii=False, indent=2), encoding="utf-8")
    time.sleep(1)  # PMDAサーバへの負荷軽減（連続アクセス時のマナー）
    return names


if __name__ == "__main__":
    import sys
    for name in (sys.argv[1:] or ["アムロジピン"]):
        info = lookup(name)
        print(f"\n=== {name} -> matched: {info['matched_name']}  format={info['format']} ===")
        print("  detail:", info["detail_url"])
        print("  pdf   :", info["pdf_url"])
        for key in ("interactions_section", "contraindicated_combinations", "caution_combinations"):
            v = info[key]
            print(f"  [{key}]", f"{len(v)} chars" if v else "なし")
