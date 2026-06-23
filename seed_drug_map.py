# -*- coding: utf-8 -*-
"""
drug_name_map.json を主要薬で一括拡充する一回限りのseedスクリプト。

各候補(日本語一般名→英語INN名)を openFDA で実在確認(_validate_name: 単独報告>0)し、
ヒットしたものだけを既存マップに併合する。綴り間違い・openFDA未収載は0件として
自動的に弾かれる（誤った英語名で“それらしい嘘の統計”を出さないための安全網）。

実行: PYTHONIOENCODING=utf-8 python seed_drug_map.py
"""
import json
from pathlib import Path

import openfda_lookup

MAP_FILE = Path(__file__).parent / "drug_name_map.json"

# 主要な外来・入院薬を中核成分名(=添付文書一般名のカタカナ)→英語INNで列挙。
# 既存マップに在るものは自動スキップ。確信のあるINN綴りのみ（検証でも弾かれるが二重の保険）。
CANDIDATES = {
    # 循環器・降圧
    "ニフェジピン": "nifedipine", "ベニジピン": "benidipine", "アゼルニジピン": "azelnidipine",
    "シルニジピン": "cilnidipine", "ニカルジピン": "nicardipine", "ベラパミル": "verapamil",
    "アテノロール": "atenolol", "メトプロロール": "metoprolol", "プロプラノロール": "propranolol",
    "ネビボロール": "nebivolol", "イルベサルタン": "irbesartan", "オルメサルタン": "olmesartan",
    "アジルサルタン": "azilsartan", "リシノプリル": "lisinopril", "ペリンドプリル": "perindopril",
    "イミダプリル": "imidapril", "テモカプリル": "temocapril", "アゾセミド": "azosemide",
    "トラセミド": "torasemide", "トリクロルメチアジド": "trichlormethiazide",
    "ヒドロクロロチアジド": "hydrochlorothiazide", "インダパミド": "indapamide",
    "エプレレノン": "eplerenone", "アミオダロン": "amiodarone", "ベプリジル": "bepridil",
    "ピルシカイニド": "pilsicainide", "フレカイニド": "flecainide", "ジソピラミド": "disopyramide",
    "メキシレチン": "mexiletine", "ニトログリセリン": "nitroglycerin", "ニコランジル": "nicorandil",
    "イバブラジン": "ivabradine", "サクビトリルバルサルタン": "sacubitril valsartan",
    # 脂質
    "エゼチミブ": "ezetimibe", "フェノフィブラート": "fenofibrate", "ベザフィブラート": "bezafibrate",
    "ペマフィブラート": "pemafibrate", "ピタバスタチン": "pitavastatin", "フルバスタチン": "fluvastatin",
    "エボロクマブ": "evolocumab",
    # 抗凝固・抗血小板
    "プラスグレル": "prasugrel", "チカグレロル": "ticagrelor", "シロスタゾール": "cilostazol",
    "ジピリダモール": "dipyridamole", "ヘパリン": "heparin", "エノキサパリン": "enoxaparin",
    "フォンダパリヌクス": "fondaparinux", "チクロピジン": "ticlopidine", "サルポグレラート": "sarpogrelate",
    # 抗菌・抗真菌・抗ウイルス・抗結核
    "ミノサイクリン": "minocycline", "ドキシサイクリン": "doxycycline", "バンコマイシン": "vancomycin",
    "リネゾリド": "linezolid", "メトロニダゾール": "metronidazole", "シプロフロキサシン": "ciprofloxacin",
    "モキシフロキサシン": "moxifloxacin", "トスフロキサシン": "tosufloxacin", "ガレノキサシン": "garenoxacin",
    "セファレキシン": "cefalexin", "セフジニル": "cefdinir", "セフトリアキソン": "ceftriaxone",
    "セフェピム": "cefepime", "メロペネム": "meropenem", "クリンダマイシン": "clindamycin",
    "スルファメトキサゾール": "sulfamethoxazole", "ボリコナゾール": "voriconazole",
    "ミカファンギン": "micafungin", "テルビナフィン": "terbinafine", "アシクロビル": "aciclovir",
    "バラシクロビル": "valaciclovir", "オセルタミビル": "oseltamivir", "アマンタジン": "amantadine",
    "リファンピシン": "rifampicin", "イソニアジド": "isoniazid", "エタンブトール": "ethambutol",
    # 消化器
    "ラベプラゾール": "rabeprazole", "ボノプラザン": "vonoprazan", "モサプリド": "mosapride",
    "イトプリド": "itopride", "レバミピド": "rebamipide", "スクラルファート": "sucralfate",
    "酸化マグネシウム": "magnesium oxide", "ルビプロストン": "lubiprostone", "リナクロチド": "linaclotide",
    "ロペラミド": "loperamide", "メサラジン": "mesalazine", "ウルソデオキシコール酸": "ursodeoxycholic acid",
    # 中枢神経・精神
    "フルボキサミン": "fluvoxamine", "デュロキセチン": "duloxetine", "ベンラファキシン": "venlafaxine",
    "ミルナシプラン": "milnacipran", "アミトリプチリン": "amitriptyline", "イミプラミン": "imipramine",
    "クロミプラミン": "clomipramine", "トラゾドン": "trazodone", "スルピリド": "sulpiride",
    "アリピプラゾール": "aripiprazole", "ブロナンセリン": "blonanserin", "パリペリドン": "paliperidone",
    "ハロペリドール": "haloperidol", "クロルプロマジン": "chlorpromazine", "ラモトリギン": "lamotrigine",
    "バルプロ酸": "valproic acid", "カルバマゼピン": "carbamazepine", "レベチラセタム": "levetiracetam",
    "ラコサミド": "lacosamide", "トピラマート": "topiramate", "フェニトイン": "phenytoin",
    "クロナゼパム": "clonazepam", "ジアゼパム": "diazepam", "エチゾラム": "etizolam",
    "ロラゼパム": "lorazepam", "ブロチゾラム": "brotizolam", "エスゾピクロン": "eszopiclone",
    "ゾピクロン": "zopiclone", "スボレキサント": "suvorexant", "レンボレキサント": "lemborexant",
    "ラメルテオン": "ramelteon", "ドネペジル": "donepezil", "メマンチン": "memantine",
    "ガランタミン": "galantamine", "リバスチグミン": "rivastigmine", "レボドパ": "levodopa",
    "プラミペキソール": "pramipexole", "ロピニロール": "ropinirole", "セレギリン": "selegiline",
    # 糖尿病
    "ボグリボース": "voglibose", "ミグリトール": "miglitol", "アカルボース": "acarbose",
    "グリメピリド": "glimepiride", "グリクラジド": "gliclazide", "ナテグリニド": "nateglinide",
    "ミチグリニド": "mitiglinide", "ピオグリタゾン": "pioglitazone", "ビルダグリプチン": "vildagliptin",
    "アログリプチン": "alogliptin", "リナグリプチン": "linagliptin", "テネリグリプチン": "teneligliptin",
    "エンパグリフロジン": "empagliflozin", "ダパグリフロジン": "dapagliflozin",
    "イプラグリフロジン": "ipragliflozin", "カナグリフロジン": "canagliflozin",
    "ルセオグリフロジン": "luseogliflozin", "トホグリフロジン": "tofogliflozin",
    "デュラグルチド": "dulaglutide", "リラグルチド": "liraglutide", "セマグルチド": "semaglutide",
    "インスリングラルギン": "insulin glargine",
    # 呼吸器・アレルギー
    "モンテルカスト": "montelukast", "プランルカスト": "pranlukast", "フェキソフェナジン": "fexofenadine",
    "ロラタジン": "loratadine", "セチリジン": "cetirizine", "レボセチリジン": "levocetirizine",
    "ビラスチン": "bilastine", "デスロラタジン": "desloratadine", "エピナスチン": "epinastine",
    "オロパタジン": "olopatadine", "ベポタスチン": "bepotastine", "クロルフェニラミン": "chlorpheniramine",
    "デキストロメトルファン": "dextromethorphan", "カルボシステイン": "carbocisteine",
    "アンブロキソール": "ambroxol", "サルブタモール": "salbutamol", "プロカテロール": "procaterol",
    "ホルモテロール": "formoterol", "サルメテロール": "salmeterol", "チオトロピウム": "tiotropium",
    "ブデソニド": "budesonide", "フルチカゾン": "fluticasone",
    # 鎮痛・リウマチ・痛風
    "イブプロフェン": "ibuprofen", "ナプロキセン": "naproxen", "メロキシカム": "meloxicam",
    "エトドラク": "etodolac", "ロルノキシカム": "lornoxicam", "ザルトプロフェン": "zaltoprofen",
    "フルルビプロフェン": "flurbiprofen", "インドメタシン": "indometacin", "ケトプロフェン": "ketoprofen",
    "モルヒネ": "morphine", "オキシコドン": "oxycodone", "フェンタニル": "fentanyl",
    "タペンタドール": "tapentadol", "ブプレノルフィン": "buprenorphine", "コデイン": "codeine",
    "アロプリノール": "allopurinol", "フェブキソスタット": "febuxostat", "コルヒチン": "colchicine",
    "ベンズブロマロン": "benzbromarone", "メトトレキサート": "methotrexate",
    "サラゾスルファピリジン": "sulfasalazine", "レフルノミド": "leflunomide",
    # ステロイド・免疫抑制
    "タクロリムス": "tacrolimus", "プレドニゾロン": "prednisolone", "ベタメタゾン": "betamethasone",
    "デキサメタゾン": "dexamethasone", "メチルプレドニゾロン": "methylprednisolone",
    "シクロスポリン": "ciclosporin", "ミコフェノール酸モフェチル": "mycophenolate mofetil",
    "アザチオプリン": "azathioprine",
    # 内分泌・泌尿器・その他
    "レボチロキシン": "levothyroxine", "チアマゾール": "thiamazole", "タムスロシン": "tamsulosin",
    "シロドシン": "silodosin", "ナフトピジル": "naftopidil", "ソリフェナシン": "solifenacin",
    "ミラベグロン": "mirabegron", "デュタステリド": "dutasteride", "フィナステリド": "finasteride",
    "シルデナフィル": "sildenafil", "タダラフィル": "tadalafil", "アレンドロン酸": "alendronate",
    "デノスマブ": "denosumab", "炭酸リチウム": "lithium carbonate",
}


def main():
    current = json.loads(MAP_FILE.read_text(encoding="utf-8"))
    todo = {jp: en for jp, en in CANDIDATES.items() if jp not in current}
    print(f"既存 {len(current)}件 / 候補 {len(CANDIDATES)}件 / 新規検証対象 {len(todo)}件\n")

    added, dropped = {}, []
    for jp, en in todo.items():
        ok = openfda_lookup._validate_name(en, polite=False)
        if ok:
            added[jp] = en
            print(f"  [OK]   {jp} -> {en}")
        else:
            dropped.append((jp, en))
            print(f"  [DROP] {jp} -> {en}  (openFDA 0件: 綴り/未収載の可能性)")

    merged = {**current, **added}
    MAP_FILE.write_text(json.dumps(merged, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"\n追加 {len(added)}件 / 除外 {len(dropped)}件 / 合計 {len(merged)}件 を drug_name_map.json に書き込み")
    if dropped:
        print("除外（要確認）:", "、".join(f"{jp}({en})" for jp, en in dropped))


if __name__ == "__main__":
    main()
