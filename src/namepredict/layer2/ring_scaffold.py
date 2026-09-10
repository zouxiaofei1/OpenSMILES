"""骨架规格 + 保留 SMILES 模板注册表 + 环解析（2026-08-15 三合一；_TEMPLATES 为唯一事实来源，派生全部 ScaffoldSpec/ScaffoldIdentity）。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol, MolFromSmiles
from rdkit import Chem
from namepredict.layer2.parent_skeleton import ParentSkeleton


# ScaffoldSpec 定义（命名/编号元数据；仅 L2 数据，不做命名组装）

@dataclass(frozen=True)
class NumberingPolicy:
    """骨架编号策略：固定编号路径、物化开关、锚点与可取代位。"""
    standard_path: tuple = ()
    materialize_plan: bool = True
    anchors: tuple[str, ...] = ()
    substitutable: frozenset[str] = frozenset()


@dataclass(frozen=True)
class ScaffoldSpec:
    """骨架规格：命名类/词干/环数/编号策略与 locant 前缀等命名元数据。"""
    id: str
    naming_class: str
    stem_en: str | None
    stem_zh: str | None
    n_rings: int
    ring: str
    retained: bool
    numbering: NumberingPolicy
    sub_rules: object | None = None
    principal_slots: object | None = None
    # 五元杂环 locant 前缀（P-61.2.4 指示氢 / P-25.1 杂原子位次）：
    # 空串不补；否则按 `prefix_nh_conditional` 决定无条件注入（1,3- 二唑）
    # 还是仅当环含未取代 NH 时注入（1H- 吡咯型）。
    locant_prefix: str = ""
    prefix_nh_conditional: bool = False

    @property
    def identity(self) -> ScaffoldIdentity:
        """返回本规格的 ScaffoldIdentity。"""
        return identity_of(self)

ScaffoldId = str


@dataclass(frozen=True)
class ScaffoldIdentity:
    """骨架身份：id、命名类、环数与环类型。"""
    id: ScaffoldId
    naming_class: str
    n_rings: int
    ring: str


def identity_of(spec) -> ScaffoldIdentity:
    """由 spec 构造 ScaffoldIdentity。"""
    return ScaffoldIdentity(spec.id, spec.naming_class, spec.n_rings, spec.ring)


# 保留母体 SMILES 模板注册表（唯一事实来源；原 specs.py + retained_templates.py 合并）
_TEMPLATES: dict[str, dict] = {
    # carbocycles
    "benzene":     {"smiles": "c1ccccc1",             "stem_en": "benzene",    "stem_zh": "苯",   "naming_class": "mono_carbo", "fused": True, "fused_prefix": ("benzo", "苯并")},
    "naphthalene": {"smiles": "c1ccc2ccccc2c1",       "stem_en": "naphthalene","stem_zh": "萘",    "naming_class": "naph_family", "fused": True, "fused_prefix": ("naphtho", "萘并")},
    "anthracene":  {"smiles": "c1ccc2cc3ccccc3cc2c1", "stem_en": "anthracene", "stem_zh": "蒽",    "naming_class": "anthra", "fused": True, "fused_prefix": ("anthra", "蒽并")},
    "phenanthrene":{"smiles": "c1ccc2c(c1)ccc1ccccc12", "stem_en": "phenanthrene","stem_zh": "菲", "naming_class": "phenanthrene", "fused": True, "fused_prefix": ("phenanthro", "菲并")},
    "pyrene":      {"smiles": "c1cc2ccc3cccc4ccc(c1)c2c34", "stem_en": "pyrene",  "stem_zh": "芘", "naming_class": "pyrene", "fused": True},
    # monocyclic heteroarenes
    "furan":       {"smiles": "c1ccoc1",    "stem_en": "furan",       "stem_zh": "呋喃",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("furo", "呋喃并")},
    "thiophene":   {"smiles": "c1ccsc1",    "stem_en": "thiophene",   "stem_zh": "噻吩",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("thieno", "噻吩并")},
    "pyrrole":     {"smiles": "c1cc[nH]c1", "stem_en": "pyrrole",     "stem_zh": "吡咯",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyridine":    {"smiles": "n1ccccc1",   "stem_en": "pyridine",    "stem_zh": "吡啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrido", "吡啶并")},
    "pyridazine":  {"smiles": "c1ccnnc1",   "stem_en": "pyridazine",  "stem_zh": "哒嗪",   "naming_class": "monohetero", "fused": True},
    "pyrimidine":  {"smiles": "c1cncnc1",   "stem_en": "pyrimidine",  "stem_zh": "嘧啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrimido", "嘧啶并")},
    "pyrazine":    {"smiles": "c1cnccn1",   "stem_en": "pyrazine",    "stem_zh": "吡嗪",   "naming_class": "monohetero", "fused": True},
    "imidazole":   {"smiles": "c1cnc[nH]1", "stem_en": "imidazole",   "stem_zh": "咪唑",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("imidazo", "咪唑并"), "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyrazole":    {"smiles": "c1ccn[nH]1", "stem_en": "pyrazole",    "stem_zh": "吡唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "oxazole":     {"smiles": "c1cocn1",    "stem_en": "oxazole",     "stem_zh": "噁唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1,3-"},
    "thiazole":    {"smiles": "c1cscn1",    "stem_en": "thiazole",    "stem_zh": "噻唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1,3-"},
    # 其余保留名杂芳环（P-25.1 表 2.8）：异噁唑/三唑/四唑/三嗪原缺失
    "isoxazole":   {"smiles": "c1ccno1",    "stem_en": "1,2-oxazole",  "stem_zh": "1,2-噁唑", "naming_class": "monohetero", "fused": True,"locant_prefix": "1,2-"},
    "triazole":    {"smiles": "c1nc[nH]n1", "stem_en": "1,2,4-triazole", "stem_zh": "1,2,4-三唑", "naming_class": "monohetero","fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "tetrazole":   {"smiles": "c1nnn[nH]1", "stem_en": "1H-tetrazole", "stem_zh": "1H-四唑", "naming_class": "monohetero","fused": True,},
    "triazine":    {"smiles": "c1ncncn1",   "stem_en": "1,3,5-triazine", "stem_zh": "1,3,5-三嗪", "naming_class": "monohetero", "fused": True,"locant_prefix": "1,3,5-"},
    # saturated monohetero rings（P-22.2.2；radical/取代基须识别为环而非开链）
    "pyrrolidine": {"smiles": "C1CCNC1",  "stem_en": "pyrrolidine", "stem_zh": "吡咯烷", "naming_class": "monohetero", "fused": True},
    "piperidine":  {"smiles": "C1CCNCC1", "stem_en": "piperidine",  "stem_zh": "哌啶",   "naming_class": "monohetero", "fused": True},
    "morpholine":  {"smiles": "C1COCCN1", "stem_en": "morpholine",  "stem_zh": "吗啉",   "naming_class": "monohetero", "fused": True},
    "piperazine":  {"smiles": "C1CNCCN1", "stem_en": "piperazine",  "stem_zh": "哌嗪",   "naming_class": "monohetero", "fused": True},
    "oxolane":     {"smiles": "C1CCOC1",  "stem_en": "oxolane",     "stem_zh": "四氢呋喃", "naming_class": "monohetero", "fused": True},
    "oxane":       {"smiles": "C1CCCOC1", "stem_en": "oxane",       "stem_zh": "氧杂环己烷", "naming_class": "monohetero", "fused": True},
    # 小环与含硫饱和杂环（P-22.2.2；3/4 元环与 S 杂环原缺失，致整块取代基丢弃）
    "oxirane":     {"smiles": "C1CO1",    "stem_en": "oxirane",     "stem_zh": "环氧乙烷", "naming_class": "monohetero","fused": True,},
    "aziridine":   {"smiles": "C1CN1",    "stem_en": "aziridine",   "stem_zh": "氮杂环丙烷", "naming_class": "monohetero","fused": True,},
    "oxetane":     {"smiles": "C1COC1",   "stem_en": "oxetane",     "stem_zh": "氧杂环丁烷", "naming_class": "monohetero","fused": True,},
    "azetidine":   {"smiles": "C1CNC1",   "stem_en": "azetidine",   "stem_zh": "氮杂环丁烷", "naming_class": "monohetero","fused": True,},
    "thiolane":    {"smiles": "C1CCSC1",  "stem_en": "thiolane",    "stem_zh": "四氢噻吩", "naming_class": "monohetero","fused": True,},
    "thiane":      {"smiles": "C1CCSCC1", "stem_en": "thiane",      "stem_zh": "四氢噻喃", "naming_class": "monohetero","fused": True,},
    # 双氧/三氧饱和环（缩醛/缩酮、溶剂类骨架）
    "dioxolane":   {"smiles": "C1COCO1",  "stem_en": "1,3-dioxolane", "stem_zh": "1,3-二氧戊环", "naming_class": "monohetero", "locant_prefix": "1,3-"},
    "dioxane":     {"smiles": "C1COCCO1", "stem_en": "1,4-dioxane",   "stem_zh": "1,4-二氧六环", "naming_class": "monohetero", "locant_prefix": "1,4-"},
    "trioxane":    {"smiles": "C1OCOCO1", "stem_en": "1,3,5-trioxane", "stem_zh": "1,3,5-三氧六环", "naming_class": "monohetero", "locant_prefix": "1,3,5-"},
    # 饱和 5 元双杂环（噁唑烷/咪唑烷/噻唑烷）
    "oxazolidine": {"smiles": "C1NCCO1",  "stem_en": "1,3-oxazolidine", "stem_zh": "1,3-噁唑烷", "naming_class": "monohetero", "locant_prefix": "1,3-"},
    "imidazolidine":{"smiles": "C1NCCN1", "stem_en": "imidazolidine",  "stem_zh": "咪唑烷", "naming_class": "monohetero"},
    "thiazolidine":{"smiles": "C1NCCS1",  "stem_en": "1,3-thiazolidine", "stem_zh": "1,3-噻唑烷", "naming_class": "monohetero", "locant_prefix": "1,3-"},
    # 部分不饱和 5/6 元杂环（P-22.2.2 加氢前缀）
    "dihydrofuran":  {"smiles": "C1C=CCO1",   "stem_en": "2,5-dihydrofuran", "stem_zh": "2,5-二氢呋喃", "naming_class": "monohetero"},
    "dihydropyran":  {"smiles": "C1=COCCC1",  "stem_en": "3,4-dihydro-2H-pyran", "stem_zh": "3,4-二氢-2H-吡喃", "naming_class": "monohetero"},
    "dihydropyrrole":{"smiles": "C1C=CCN1",   "stem_en": "2,5-dihydro-1H-pyrrole", "stem_zh": "2,5-二氢-1H-吡咯", "naming_class": "monohetero", "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "dihydroimidazole":{"smiles": "C1=NCCN1", "stem_en": "4,5-dihydro-1H-imidazole", "stem_zh": "4,5-二氢-1H-咪唑", "naming_class": "monohetero", "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "dihydrothiazole":{"smiles": "C1=NCCS1",  "stem_en": "4,5-dihydro-1,3-thiazole", "stem_zh": "4,5-二氢-1,3-噻唑", "naming_class": "monohetero", "locant_prefix": "1,3-"},
    # fused 5+6
    "indole":         {"smiles": "c1ccc2[nH]ccc2c1", "stem_en": "1H-indole",      "stem_zh": "吲哚",     "naming_class": "fused56", "fused": True, "fused_stem": ("indole", "吲哚"), "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "indazole":       {"smiles": "c1ccc2cn[nH]c2c1", "stem_en": "indazole",       "stem_zh": "吲唑",     "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "benzimidazole":  {"smiles": "c1ccc2[nH]cnc2c1", "stem_en": "benzimidazole",  "stem_zh": "苯并咪唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "benzofuran":     {"smiles": "c1ccc2occc2c1",    "stem_en": "benzofuran",     "stem_zh": "苯并呋喃", "naming_class": "fused56", "fused": True, "locant_prefix": "1-"},
    "benzothiophene": {"smiles": "c1ccc2sccc2c1",    "stem_en": "benzothiophene", "stem_zh": "苯并噻吩", "naming_class": "fused56", "fused": True, "locant_prefix": "1-"},
    "benzothiazole":  {"smiles": "c1ccc2scnc2c1",    "stem_en": "benzothiazole",  "stem_zh": "苯并噻唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-"},
    "benzoxazole":    {"smiles": "c1ccc2ocnc2c1",    "stem_en": "benzoxazole",    "stem_zh": "苯并噁唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-"},
    # 三环 6+5+6（13 原子）：两个苯环融合吡咯，N9 邻位桥头 4a/9a。
    "carbazole":      {"smiles": "c1ccc2[nH]c3ccccc3c2c1", "stem_en": "carbazole", "stem_zh": "咔唑", "naming_class": "carbazole", "fused": True, "locant_prefix": "9H-", "prefix_nh_conditional": True},
    # 三环 6+6+6（14 原子）：中间吡啶/含 S 环两侧苯环融合。
    "acridine":       {"smiles": "c1ccc2nc3ccccc3cc2c1",   "stem_en": "acridine",     "stem_zh": "吖啶",   "naming_class": "acridine", "fused": True},
    "phenothiazine":  {"smiles": "c1ccc2Sc3ccccc3Nc2c1",   "stem_en": "phenothiazine", "stem_zh": "吩噻嗪", "naming_class": "phenothiazine", "fused": True, "locant_prefix": "10H-", "prefix_nh_conditional": True},
    # 双环 5+6（9 原子）：苯环并二氧戊环，O1/C2/O3。
    "benzodioxole":   {"smiles": "c1ccc2OCOc2c1",          "stem_en": "benzodioxole", "stem_zh": "苯并二氧杂环戊烯", "naming_class": "benzodioxole", "fused": True, "locant_prefix": "1,3-"},
    # fused 6+6
    "quinoline":    {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline",    "stem_zh": "喹啉",   "naming_class": "naph_family", "fused": True},
    "isoquinoline": {"smiles": "c1nccc2ccccc21", "stem_en": "isoquinoline", "stem_zh": "异喹啉", "naming_class": "naph_family", "fused": True},
    "quinazoline":  {"smiles": "c1ccc2ncncc2c1", "stem_en": "quinazoline",  "stem_zh": "喹唑啉", "naming_class": "naph_family", "fused": True},
    "quinoxaline":  {"smiles": "c1ccc2nccnc2c1", "stem_en": "quinoxaline",  "stem_zh": "喹喔啉", "naming_class": "naph_family", "fused": True},
    # 保留名（表 2.8）：purine=嘌呤（5+6，特殊编号 1–9，PIN 7H-purine，P-25 表 2.8 第16项）；
    # pteridine=蝶啶（6+6 四 N，naph-family 编号，N 在 1/3/5/8、CH 在 2/4/6/7）。无环外 =O，故可入表。
    "purine":       {"smiles": "c1ncc2[nH]cnc2n1", "stem_en": "7H-purine",     "stem_zh": "嘌呤",   "naming_class": "purine", "fused": True, "fused_stem": ("purine", "嘌呤"), "locant_prefix": "7H-", "prefix_nh_conditional": True},
    "pteridine":    {"smiles": "c1cnc2ncncc2n1",    "stem_en": "pteridine",    "stem_zh": "蝶啶",   "naming_class": "naph_family", "fused": True},
    # NOTE: carbonyl mothers（benzoquinone / anthraquinone / chromenone /
    # ortho_benzoquinone）不入表：模板含环外 =O，匹配集会超出环系统原子集。
}


def component_stem(sid: str) -> tuple[str, str] | None:
    """稠合组分词干 (en, zh)：只有 `fused=True` 的模板可作稠合命名零件，否则 None。

    多数等于母体词干；indole/purine 以 `fused_stem` 覆盖去掉指示氢前缀（1H-/7H-），
    因为稠合前缀取去指示氢的组分名（indolo 而非 1H-indolo）。
    """
    entry = _TEMPLATES.get(sid)
    if not entry or not entry.get("fused"):
        return None
    return entry.get("fused_stem") or (entry["stem_en"], entry["stem_zh"])


def retained_fusion_prefix(sid: str) -> tuple[str, str] | None:
    """附加组分的保留稠合前缀 (en, zh)（P-25.3.2.2.3）；无登记则 None（调用方走通用规则）。"""
    entry = _TEMPLATES.get(sid)
    return entry.get("fused_prefix") if entry else None

# 保留 fused 母体的固定编号标签（P-25.4）：融合桥头用字母 locant（3a/7a、4a/8a）。
FUSED56_LABELS: tuple[str, ...] = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")
# purine（嘌呤，5+6 九原子）保留传统编号：桥头为 C4/C5 得纯数字 4、5，无 a/b 字母位（P-25.3.3 例外）。
PURINE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "5", "6", "7", "8", "9")
# carbazole（13 原子）：N9，苯环 1-4 / 5-8，桥头 4a/8a/9a/9b（P-25.4.1.4）。
CARBAZOLE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a", "9", "9a", "9b")
# acridine（14 原子）：N10，对位 C9 连接取代基；苯环 1-4 / 5-8。
ACRIDINE_LABELS: tuple[str, ...] = ("1","2","3","4","4a","5","6","7","8","8a","9","9a","10","10a")
# phenothiazine（14 原子）：S5、N10；苯环 1-4 / 6-9。
PHENOTHIAZINE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "9", "9a", "10", "10a", "10b")
NAPH_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")
# phenanthrene(14 原子): 端环 1-4 / 5-8, 桥头 4a/4b/8a/8b, 中环 9/10(P-25.4.1 传统编号)。
PHENANTHRENE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "8b", "9", "10")
# pyrene(16 原子): 外周 1-10, 稠合碳 3a/5a/8a/8b/10a/10b(P-25.3.3.3.1 推荐编号)。
PYRENE_LABELS: tuple[str, ...] = ("1", "2", "3", "3a", "4", "5", "5a", "6", "7", "8", "8a", "8b", "9", "10", "10a", "10b")

# 不对称 fused 环（含杂原子）与 1,3-二唑的标准编号：模板原子索引按固定 locant 顺序排列。
# 起点为最优先杂原子，沿环编号绕开融合桥头（桥头只得字母位）。
# 对称碳环（naphthalene/anthracene）与单杂环（pyrrole/pyridine 等单杂原子）不在此列，
# 走 P-14.4 通用枚举；1,3-二唑（咪唑/吡唑/噻唑/噁唑）IUPAC 编号固定（N 必须得 1,3 位、
# 带 H/取代基的 N 为 1），通用枚举会被 principal 最小化翻转，故登记标准编号。
# _STANDARD_ORDERS[spec] = 模板原子按 locant 顺序；_STANDARD_LABELS[spec] = 对应 locant 标签。
_STANDARD_ORDERS: dict[str, tuple[int, ...]] = {
    # naph_family（10 原子）：1,2,3,4,4a,5,6,7,8,8a
    "quinoline":    (4, 5, 6, 7, 8, 9, 0, 1, 2, 3),
    "quinazoline":  (4, 5, 6, 7, 8, 9, 0, 1, 2, 3),
    # purine（嘌呤）9 原子：模板 c1ncc2[nH]cnc2n1 原子按 locant 1–9 顺序（N1、C2、N3、C4、C5、C6、N7、C8、N9；
    # 桥头 C4/C5 得数字位，N7 为指示氢所在）。locant→模板原子：1→1,2→0,3→8,4→7,5→3,6→2,7→4,8→5,9→6。
    "purine":       (1, 0, 8, 7, 3, 2, 4, 5, 6),
    # pteridine（蝶啶）10 原子：模板 c1cnc2ncncc2n1 原子沿外周按 1,2,3,4,4a,5,6,7,8,8a 编号（N 在 1/3/5/8）。
    "pteridine":    (4, 5, 6, 7, 8, 9, 0, 1, 2, 3),
    # fused56（9 原子）：1,2,3,3a,4,5,6,7,7a；杂原子(1)走远离桥头方向，苯环从 3a 起
    "indole":       (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzofuran":   (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzothiophene": (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzothiazole": (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "benzoxazole":  (4, 5, 6, 7, 8, 0, 1, 2, 3),
    "indazole":     (6, 5, 4, 3, 2, 1, 0, 8, 7),
    # phenanthrene(14 原子): 模板原子按标准 locant 序(端环1-4 + 中环9/10)。
    "phenanthrene": (9, 10, 11, 12, 13, 8, 5, 0, 1, 2, 3, 4, 6, 7),
    # pyrene(16 原子): 外周1-10 + 稠合碳 3a/5a/8a/8b/10a/10b(fused_numbering 方向)。
    "pyrene":       (6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 14, 3, 4, 5, 15),
    # carbazole（13 原子）：1,2,3,4,4a,5,6,7,8,8a,9,9a,9b
    "carbazole":    (6, 7, 8, 9, 5, 12, 0, 1, 2, 11, 4, 3, 10),
    # acridine（14 原子）：1,2,3,4,4a,5,6,7,8,8a,9(N 对位),10(N),10a,10b
    "acridine":     (9,8,7,6,5,2,1,0,13,12,11,10,4,3),
    # phenothiazine（14 原子）：1,2,3,4,4a,5(S),6,7,8,9,9a,10(N),10a,10b
    "phenothiazine": (0, 1, 2, 3, 12, 4, 5, 6, 7, 8, 10, 11, 13, 9),
    # benzodioxole（9 原子）：1(O),2(C),3(O),3a,4,5,6,7,7a
    "benzodioxole": (4, 5, 6, 7, 3, 8, 0, 1, 2),
    # monohetero 1,3-二唑（5 原子）：N1(带 H/取代)/C2/N3/C4/C5。
    # 噻唑/噁唑杂原子异元素（S/O vs N）模板匹配无歧义，登记固定编号；
    # 咪唑/吡唑双 N 对称（模板 [nH] 对两个 N 可互换匹配），N1 须动态取
    # 带取代基/H 的 N，走 numbering_engine._narrow_hetero_ring 的 (c) 同元素 N 收窄。
    "thiazole":     (2, 3, 4, 0, 1),
    "oxazole":      (2, 3, 4, 0, 1),
}
_STANDARD_LABELS: dict[str, tuple[str, ...]] = {
    "quinoline": NAPH_LABELS,
    "quinazoline": NAPH_LABELS,
    "pteridine": NAPH_LABELS,
    "purine": PURINE_LABELS,
    "indole": FUSED56_LABELS,
    "benzofuran": FUSED56_LABELS,
    "benzothiophene": FUSED56_LABELS,
    "benzothiazole": FUSED56_LABELS,
    "benzoxazole": FUSED56_LABELS,
    "indazole": FUSED56_LABELS,
    "phenanthrene": PHENANTHRENE_LABELS,
    "pyrene": PYRENE_LABELS,
    "carbazole": CARBAZOLE_LABELS,
    "acridine": ACRIDINE_LABELS,
    "phenothiazine": PHENOTHIAZINE_LABELS,
    "benzodioxole": FUSED56_LABELS,
    "thiazole": ("1", "2", "3", "4", "5"),
    "oxazole": ("1", "2", "3", "4", "5"),
}

# 查询子结构与元素签名，import 时构建一次。
_Q: dict[str, Mol] = {sid: MolFromSmiles(entry["smiles"]) for sid, entry in _TEMPLATES.items()}


def _hydrogenated(mol: Mol) -> Mol | None:
    """返回完全氢化的分子副本（清芳香性 → 全键改单键 → sanitize 补满隐式 H），骨架不可 sanitize 时 None。"""
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    out = rw.GetMol()
    if Chem.SanitizeMol(out, catchErrors=True) != Chem.SanitizeFlags.SANITIZE_NONE:
        return None
    return out


# 完全氢化模板（键级/芳香性抹平的骨架），供 _match_with_map 匹配加氢衍生物（P-25.3.4）。
# 苯例外（P-31.2.3.1）：单环 mancude 碳环的加氢衍生物用 cyclohexene/cyclohexadiene/cyclohexane
# （P-31.1.3 ene 词尾），不写 hydrobenzene；故 benzene 模板不参与氢化骨架匹配，环己烷/环己烯类
# 继续走 carbocycle 单环路径。单环杂环不豁免：无饱和保留名时仍以母体+hydro 表达（P-31.2.3.2）。
_Q_H: dict[str, Mol] = {
    sid: h for sid, q in _Q.items()
    if _TEMPLATES[sid]["naming_class"] != "mono_carbo" and (h := _hydrogenated(q)) is not None
}


_KEKULE_ATOMS: dict[str, frozenset[int]] = {}


def _kekule_double_atoms(sid: str) -> frozenset[int]:
    """模板 Kekulé 双键端点原子集（缓存）：芳香键须化为确定双键，否则稠合单键（萘 4a-8a）会被误当不饱和度。"""
    hit = _KEKULE_ATOMS.get(sid)
    if hit is not None:
        return hit
    q = _Q.get(sid)
    out: set[int] = set()
    if q is not None:
        m = Chem.Mol(q)
        try:
            Chem.Kekulize(m, clearAromaticFlags=True)
            out = {b.GetBeginAtomIdx() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE}
            out |= {b.GetEndAtomIdx() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE}
        except Exception:
            out = set()
    _KEKULE_ATOMS[sid] = frozenset(out)
    return _KEKULE_ATOMS[sid]


def hydrogenated_atoms(mol: Mol, scaffold_id: str, match) -> frozenset[int]:
    """match（模板原子→分子原子）下被加氢的分子原子集：模板某原子承载不饱和双键（Kekulé）而分子中该位已非芳香者（P-31.2.2：hydro 修饰源于双键的饱和）。用原子芳香性而非键级，桥头/带取代基饱和碳均可正确归属。"""
    if not match or mol is None or scaffold_id not in _Q:
        return frozenset()
    out: set[int] = set()
    for qi in _kekule_double_atoms(scaffold_id):
        if qi >= len(match):
            continue
        mi = match[qi]
        if mi < mol.GetNumAtoms() and all(
                b.GetBondType() == Chem.BondType.SINGLE for b in mol.GetAtomWithIdx(mi).GetBonds()):
            out.add(mi)  # 分子中该位已无多重键（原带双键、现饱和）才是加氢位；残留芳香/多重键者未加氢（hybridization 对 NH 会误报 SP2，故查键级）
    return frozenset(out)


def _elem_sig(mol: Mol, atom_ids) -> frozenset:
    """原子集的元素组成签名（(Z, 计数) 冻结集合）。"""
    return frozenset(Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids).items())


_TEMPLATE_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q.items()}


def _spec_from_template(sid: str, entry: dict) -> ScaffoldSpec:
    """由模板条目派生 ScaffoldSpec（n_rings/ring 从 smiles 自动算）。"""
    q = _Q[sid]
    n_rings = len(q.GetRingInfo().AtomRings())
    ring = "carbo" if all(q.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in range(q.GetNumAtoms())) else "hetero"
    std_labels = _STANDARD_LABELS.get(sid, ())
    numbering = NumberingPolicy(standard_path=std_labels)
    return ScaffoldSpec(
        id=sid, naming_class=entry["naming_class"],
        stem_en=entry["stem_en"], stem_zh=entry["stem_zh"],
        n_rings=n_rings, ring=ring, retained=True,
        numbering=numbering,
        locant_prefix=entry.get("locant_prefix", ""),
        prefix_nh_conditional=bool(entry.get("prefix_nh_conditional")),
    )


_ALL_SPECS: tuple[ScaffoldSpec, ...] = tuple(
    _spec_from_template(sid, entry) for sid, entry in _TEMPLATES.items()
)
_BY_ID: dict[str, ScaffoldSpec] = {s.id: s for s in _ALL_SPECS}
_IDENTITIES: dict[str, ScaffoldIdentity] = {s.id: s.identity for s in _ALL_SPECS}


def get_identity(spec_id: str) -> ScaffoldIdentity | None:
    """按 id 查 ScaffoldIdentity（无则 None）。"""
    return _IDENTITIES.get(spec_id)


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    """按 id 查 ScaffoldSpec（无则 None）。"""
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    """返回全部 ScaffoldSpec（由 _TEMPLATES 派生）。"""
    return _ALL_SPECS


def numbering_scaffold_facts(spec_id: str | None, atom_count: int) -> dict | None:
    """物化纯母体 facts；``relative_stereo`` 预留环面约束。"""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.numbering.materialize_plan:
        return None
    labels = spec.numbering.standard_path
    return None if not labels or len(labels) != atom_count else {
        "scaffold_id": spec.id, "labels": labels, "relative_stereo": None,
    }


def fused56_kind_ids() -> frozenset[str]:
    """使用 fused56 位次标签的 kind / scaffold id。"""
    return kind_ids_for("fused56")


def naph_kind_ids() -> frozenset[str]:
    """使用 naph 位次标签的 kind / scaffold id。"""
    return kind_ids_for("naph_family")


def monohetero_kind_ids() -> frozenset[str]:
    """单杂环保留骨架的 kind / scaffold id。"""
    return kind_ids_for("monohetero")


def kind_ids_for(naming_class: str) -> frozenset[str]:
    """给定 naming_class（fused56 / naph_family / monohetero / …）的 id 集合。"""
    return frozenset(s.id for s in _ALL_SPECS if s.naming_class == naming_class)

RetainedEntry = dict


def _entry_for(sid: str) -> RetainedEntry | None:
    """由模板条目派生保留条目（拓扑字段 + 词干）。"""
    entry = _TEMPLATES.get(sid)
    if entry is None:
        return None
    q = _Q[sid]
    hetero_Z = tuple(sorted(
        q.GetAtomWithIdx(i).GetAtomicNum() for i in range(q.GetNumAtoms())
        if q.GetAtomWithIdx(i).GetAtomicNum() != 6
    ))
    n_rings = len(q.GetRingInfo().AtomRings())
    return {
        "kind": sid,
        "n_rings": n_rings,
        "n_atoms": q.GetNumAtoms(),
        "hetero_Z": hetero_Z,
        "topology": "fused" if n_rings > 1 else "mono",
        "aromatic": all(q.GetAtomWithIdx(i).GetIsAromatic() for i in range(q.GetNumAtoms())),
        "en": entry["stem_en"],
        "zh": entry["stem_zh"],
    }


def registry() -> dict[str, RetainedEntry]:
    """返回带词干的保留拓扑注册表（由 _TEMPLATES 派生）。"""
    return {sid: _entry_for(sid) for sid in _TEMPLATES}


def get_entry(scaffold_id: str) -> RetainedEntry | None:
    """按 scaffold_id 查带词干的保留条目。"""
    return _entry_for(scaffold_id)


def match_systems(info: dict) -> list[tuple[str, dict, RetainedEntry]]:
    """返回匹配系统的 (scaffold_id, ring_system, entry)（模板子图同构语义）。"""
    out: list[tuple[str, dict, RetainedEntry]] = []
    for system in info.get("ring_systems") or []:
        sid = match_retained(info, system.get("atom_ids") or ())
        if sid:
            out.append((sid, system, _entry_for(sid)))
    return out


def match_scaffold_ids(info: dict) -> list[str]:
    """返回匹配的 scaffold id 列表。"""
    return [sid for sid, _, _ in match_systems(info)]


def match_retained(info: dict, atom_ids) -> str | None:
    """返回模板精确覆盖 atom_ids 的保留母体 sid，无命中 None；元素签名预过滤后子图同构，同命中取表序第一个（防御性兜底）。"""
    hit = _match_with_map(info, atom_ids)
    return hit[0] if hit else None


def _match_with_map(info: dict, atom_ids) -> tuple[str, tuple[int, ...]] | None:
    """模板精确覆盖 atom_ids 时返回 (sid, match)；match[i] 供 standard_path 把模板原子映射到分子原子。"""
    mol = info["mol"]
  
    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)  
    print(Chem.MolToSmiles(mol))
    for sid, q in _Q.items():
        if _TEMPLATE_ELEM[sid] != elem:
            continue
        for m in mol.GetSubstructMatches(q, uniquify=True):  # print(sid,m)
            if set(m) == atoms:  # print("yes")
                return sid, m
    mol_h = _hydrogenated(mol)  # 精确匹配失败后按完全氢化骨架再比对（加氢衍生物，P-25.3.4）
    if mol_h is not None:
        for sid, qh in _Q_H.items():
            if _TEMPLATE_ELEM[sid] != elem:
                continue
            for m in mol_h.GetSubstructMatches(qh, uniquify=True):
                if set(m) == atoms:
                    return sid, m
    return None


def locant_prefix(spec_id: str | None) -> tuple[str, str, bool]:
    """返回 scaffold 的 locant 前缀 (en, zh, nh_conditional)，无则空；1,3- 二唑无条件注入词干、1H- 吡咯型仅环含未取代 NH 时注入。"""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.locant_prefix:
        return "", "", False
    return spec.locant_prefix, spec.locant_prefix, spec.prefix_nh_conditional


def standard_chain(spec_id: str | None, match: tuple[int, ...] | None) -> list[int] | None:
    """把 fused 模板固定编号映射到分子，返回原子按标准 locant 顺序的列表；无标准顺序/match 长度不符返回 None（走 P-14.4 通用枚举）。"""
    if not match:
        return None
    order = _STANDARD_ORDERS.get(spec_id or "")
    if not order or len(order) != len(match):
        return None
    return [match[t] for t in order]


# 环解析

def _matched_id(info: dict, skeleton: ParentSkeleton) -> str | None:
    """按模板子图同构匹配骨架的 scaffold id。"""
    return match_retained(info, skeleton.atom_ids)


def _generic_carbocycle(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    """无模板命中时的通用环身份兜底：全碳单/多环→carbocycle；非全碳芳香多环→fused_hetero；其余→None（显式失败，避免当开链烷基错名）。"""
    mol = info["mol"]
    all_carbon = all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in skeleton.atom_ids)
    if not all_carbon:
        atoms = set(skeleton.atom_ids)
        # if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atoms):
        n_rings = sum(1 for ring in mol.GetRingInfo().AtomRings() if set(ring) <= atoms)
        if n_rings >= 2:
            return ScaffoldIdentity("fused_hetero", "fused_hetero", n_rings, "hetero")
        return None
    return ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")


def resolve_ring_scaffold(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    """解析骨架的 scaffold 身份（显式/模板匹配/通用兜底：全碳→carbocycle，非全碳多环→fused_hetero，非全碳单环→None）。"""
    direct = get_identity(skeleton.scaffold_id or "")
    if direct:
        return direct
    sid = _matched_id(info, skeleton)
    if sid:
        identity = get_identity(sid)
        if identity:
            return identity
    return _generic_carbocycle(info, skeleton)
