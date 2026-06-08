# -*- coding: utf-8 -*-
"""
PMDA「医療用医薬品 添付文書等情報検索」から薬剤名で検索し、
添付文書（PDF/XML/SGML）を取得できるか検証するスクリプト。

検証目的:
  1. 薬剤名 -> 検索結果一覧 -> 詳細ページ -> 添付文書ファイルへの導線が
     スクリプトから機械的にたどれるか
  2. 添付文書本文に「相互作用」「併用禁忌」「CYP」等の記載が
     抽出できるか（インタラクションチェッカーの素材になり得るか）
"""
import re
import sys
import time
from pathlib import Path

import requests
from bs4 import BeautifulSoup

BASE = "https://www.pmda.go.jp"
SEARCH_URL = f"{BASE}/PmdaSearch/iyakuSearch/"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36",
}

OUT_DIR = Path(__file__).parent / "pmda_sample"
OUT_DIR.mkdir(exist_ok=True)


def search_drug(session: requests.Session, name: str):
    """薬剤名で検索し、検索結果ページのHTMLを返す

    検証で判明した点:
      - 検索フォームは Struts/Spring 系で、隠しフィールド(`_xxx`)や
        チェックボックスの初期状態を含む「フォーム全項目」を送らないと
        サーバ側のバリデーションを通過せず、フォーム自体が再表示される。
      - 検索ボタンは <input type="image"> なので `btnA.x` / `btnA.y` で送る。
      => 一番確実なのは「フォームHTMLをそのまま読み取り、必要な項目だけ
         上書きして送り返す」方式。
    """
    r = session.get(SEARCH_URL, headers=HEADERS, timeout=30)
    r.raise_for_status()
    soup = BeautifulSoup(r.text, "lxml")
    form = soup.find("form", {"id": "iyakuSearchForm"})

    overrides = {
        "nameWord": name,
        "iyakuHowtoNameSearchRadioValue": "1",  # 1:一般名+販売名 / 2:一般名 / 3:販売名
        "howtoMatchRadioValue": "2",            # 1:部分一致 / 2:前方一致
    }
    seen = set()
    payload = []
    for inp in form.find_all(["input", "select", "textarea"]):
        name_attr = inp.get("name")
        if not name_attr:
            continue
        tag = inp.name
        if tag == "input":
            itype = (inp.get("type") or "text").lower()
            if itype == "image" or itype == "button":
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

    # 画像ボタン（検索実行）のクリック座標を模す
    payload.append(("btnA.x", "15"))
    payload.append(("btnA.y", "10"))

    r2 = session.post(SEARCH_URL, data=payload, headers={**HEADERS, "Referer": SEARCH_URL}, timeout=30)
    r2.raise_for_status()
    return r2.text


def extract_general_list_links(html: str):
    soup = BeautifulSoup(html, "lxml")
    links = []
    for a in soup.select('a[href*="GeneralList"]'):
        href = a.get("href", "")
        text = a.get_text(strip=True)
        if href and text:
            links.append((text, requests.compat.urljoin(BASE + "/PmdaSearch/iyakuSearch/", href)))
    # 重複除去
    seen = set()
    uniq = []
    for t, u in links:
        if u not in seen:
            seen.add(u)
            uniq.append((t, u))
    return uniq


def extract_doc_links(detail_html: str, detail_url: str):
    """詳細ページから添付文書PDF / XML / SGML へのリンクを抽出"""
    soup = BeautifulSoup(detail_html, "lxml")
    pdf, xml, sgml, other = [], [], [], []
    for a in soup.find_all("a", href=True):
        href = a["href"]
        full = requests.compat.urljoin(detail_url, href)
        low = full.lower()
        if "resultdatasetpdf" in low or low.endswith(".pdf"):
            pdf.append(full)
        elif low.endswith(".xml") or "xml" in low:
            xml.append(full)
        elif low.endswith(".sgml") or "sgml" in low:
            sgml.append(full)
    return pdf, xml, sgml


def analyze_pdf_text(pdf_bytes: bytes):
    """PDFからテキストを抽出し、相互作用関連キーワードの有無を調べる"""
    import fitz  # PyMuPDF
    doc = fitz.open(stream=pdf_bytes, filetype="pdf")
    text = "\n".join(page.get_text() for page in doc)
    doc.close()

    keywords = ["相互作用", "併用禁忌", "併用注意", "CYP", "P-糖蛋白", "P糖蛋白", "薬物動態"]
    found = {kw: len(re.findall(re.escape(kw), text)) for kw in keywords}
    return text, found


def run(drug_name: str):
    print(f"\n=== 検索: {drug_name} ===")
    with requests.Session() as session:
        html = search_drug(session, drug_name)
        links = extract_general_list_links(html)
        print(f"検索結果リンク数: {len(links)}")
        for t, u in links[:10]:
            print(f"  - {t}  ->  {u}")

        if not links:
            print("  (検索結果が0件、またはPOSTがサーバ側検索を発火していない可能性)")
            (OUT_DIR / f"{drug_name}_search_result.html").write_text(html, encoding="utf-8")
            print(f"  検索結果ページを保存しました: {drug_name}_search_result.html (中身を目視確認用)")
            return

        # 最初の1件の詳細ページを開く
        title, detail_url = links[0]
        print(f"\n詳細ページを開きます: {title}\n  {detail_url}")
        with requests.Session() as s2:
            s2.get(SEARCH_URL, headers=HEADERS, timeout=30)
            dr = s2.get(detail_url, headers={**HEADERS, "Referer": SEARCH_URL}, timeout=30)
            dr.raise_for_status()
            pdf_links, xml_links, sgml_links = extract_doc_links(dr.text, detail_url)

            print(f"  PDFリンク: {len(pdf_links)}件")
            for u in pdf_links[:5]:
                print(f"    - {u}")
            print(f"  XMLリンク: {len(xml_links)}件  / SGMLリンク: {len(sgml_links)}件")

            if pdf_links:
                pdf_url = pdf_links[0]
                print(f"\n  添付文書PDFを取得して本文解析します: {pdf_url}")
                pr = s2.get(pdf_url, headers={**HEADERS, "Referer": detail_url}, timeout=60)
                pr.raise_for_status()
                fname = OUT_DIR / f"{drug_name}_tenpu.pdf"
                fname.write_bytes(pr.content)
                print(f"    -> 保存: {fname.name} ({len(pr.content):,} bytes)")

                text, found = analyze_pdf_text(pr.content)
                (OUT_DIR / f"{drug_name}_tenpu.txt").write_text(text, encoding="utf-8")
                print("    抽出テキスト中のキーワード出現回数:")
                for kw, cnt in found.items():
                    print(f"      {kw}: {cnt}")

                # 「相互作用」セクション付近を抜粋表示
                m = re.search(r".{0,30}相互作用.{0,400}", text, re.S)
                if m:
                    print("\n    --- 「相互作用」記載の抜粋 ---")
                    print("    " + m.group(0).replace("\n", " ")[:400])
            else:
                print("  PDFリンクが見つかりませんでした。詳細ページHTMLを保存します。")
                (OUT_DIR / f"{drug_name}_detail.html").write_text(dr.text, encoding="utf-8")


if __name__ == "__main__":
    targets = sys.argv[1:] or ["ロキソプロフェン", "ワルファリン", "クラリスロマイシン"]
    for name in targets:
        try:
            run(name)
        except Exception as e:
            print(f"  [ERROR] {name}: {e!r}")
        time.sleep(2)
