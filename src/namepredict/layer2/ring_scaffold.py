"""骨架规格 + 保留 SMILES 模板注册表 + 环解析（2026-08-15 三合一；_TEMPLATES 为唯一事实来源，派生全部 ScaffoldSpec/ScaffoldIdentity）。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol, MolFromSmarts, MolFromSmiles
from rdkit import Chem
from namepredict.tools import memo
from namepredict.layer2.parent_skeleton import ParentSkeleton
from namepredict.layer1.ring_systems import sssr_rings


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
# xanthene/thioxanthene(14 原子): 外周 1-8, 中央碳 9、O/S 10, 稠合碳 4a/8a/9a/10a(P-25.3.3 传统编号；位次形态同 acridine)。
XANTHENE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a", "9", "9a", "10", "10a")
# cyclopenta[a]phenanthrene(17 原子): 传统甾体编号 1-17 全数字(10/13 为角甲基碳, 非字母桥头)。
STEROID_LABELS: tuple[str, ...] = tuple(str(i) for i in range(1, 18))


# 固定编号登记（P-25.4/P-25.3.3）：模板条目 `standard` = (labels, order) 二元组，
# order 是模板原子按 locant 顺序的下标排列，labels 为对应 locant 标签，两者长度均等于模板原子数。
# 未登记者走 P-14.4 通用枚举（对称碳环 naphthalene/anthracene、单杂环 pyrrole/pyridine 等）。
# order 起点为最优先杂原子，沿环编号绕开融合桥头（桥头只得字母位）；1,3-二唑（咪唑/吡唑/噻唑/噁唑）
# 的编号为 IUPAC 固定（N 须得 1,3 位、带 H/取代基的 N 为 1），通用枚举会被 principal 最小化翻转，故一并登记。
# order 与 labels 同条目、且与 smiles 同表，避免改 SMILES 后编号静默错位（见 _validate_standard_fields）。

# 保留母体 SMILES 模板注册表（唯一事实来源；原 specs.py + retained_templates.py 合并）
_TEMPLATES: dict[str, dict] = {
    # carbocycles
    "benzene":     {"smiles": "c1ccccc1",             "stem_en": "benzene",    "stem_zh": "苯",   "naming_class": "mono_carbo", "fused": True, "fused_prefix": ("benzo", "苯并")},
    "naphthalene": {"smiles": "c1ccc2ccccc2c1",       "stem_en": "naphthalene","stem_zh": "萘",    "naming_class": "naph_family", "fused": True, "fused_prefix": ("naphtho", "萘并")},
    "anthracene":  {"smiles": "c1ccc2cc3ccccc3cc2c1", "stem_en": "anthracene", "stem_zh": "蒽",    "naming_class": "anthra", "fused": True, "fused_prefix": ("anthra", "蒽并")},
    "phenanthrene":{"smiles": "c1ccc2c(c1)ccc1ccccc12", "stem_en": "phenanthrene","stem_zh": "菲", "naming_class": "phenanthrene", "fused": True, "fused_prefix": ("phenanthro", "菲并"), "standard": (PHENANTHRENE_LABELS, (9, 10, 11, 12, 13, 8, 5, 0, 1, 2, 3, 4, 6, 7))},
    "pyrene":      {"smiles": "c1cc2ccc3cccc4ccc(c1)c2c34", "stem_en": "pyrene",  "stem_zh": "芘", "naming_class": "pyrene", "fused": True, "standard": (PYRENE_LABELS, (6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 14, 3, 4, 5, 15))},
    # 表 2.7 第 19 项：茚（PIN 1H-indene，1 位为 CH2 故带指示氢）。5+6 稠合碳环，
    # 位次形态同吲哚/苯并呋喃（1,2,3,3a,4..7,7a），并入 fused56 编号类。
    "indene":      {"smiles": "C1=CCc2ccccc21", "stem_en": "1H-indene", "stem_zh": "1H-茚", "naming_class": "fused56", "fused": True, "fused_stem": ("indene", "茚"), "standard": (FUSED56_LABELS, (2, 1, 0, 8, 7, 6, 5, 4, 3))},
    # 表 2.7 第 8 项：䓛（PIN chrysene；中文按库内约定用「屈」）。6+6+6+6 四环稠烃，
    # 外周 1-6 / 7-12，六桥头 6a,6b,6c,12a,12b,12c（P-25.3.3），暂走优选取向自动编号。
    "chrysene":    {"smiles": "c1ccc2c(c1)ccc1c3ccccc3ccc21", "stem_en": "chrysene", "stem_zh": "屈", "naming_class": "chrysene", "fused": True, "fused_prefix": ("chryseno", "䓛并")},
    # monocyclic heteroarenes
    "furan":       {"smiles": "c1ccoc1",    "stem_en": "furan",       "stem_zh": "呋喃",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("furo", "呋喃并")},
    "thiophene":   {"smiles": "c1ccsc1",    "stem_en": "thiophene",   "stem_zh": "噻吩",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("thieno", "噻吩并")},
    "pyrrole":     {"smiles": "c1cc[nH]c1", "stem_en": "pyrrole",     "stem_zh": "吡咯",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyridine":    {"smiles": "n1ccccc1",   "stem_en": "pyridine",    "stem_zh": "吡啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrido", "吡啶并")},
    "pyridazine":  {"smiles": "c1ccnnc1",   "stem_en": "pyridazine",  "stem_zh": "哒嗪",   "naming_class": "monohetero", "fused": True},
    "pyrimidine":  {"smiles": "c1cncnc1",   "stem_en": "pyrimidine",  "stem_zh": "嘧啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrimido", "嘧啶并")},
    "pyrazine":    {"smiles": "c1cnccn1",   "stem_en": "pyrazine",    "stem_zh": "吡嗪",   "naming_class": "monohetero", "fused": True},
    # 吡喃（P-25.1 表 2.8 保留名，6 元含氧 mancude 母体；原表缺失）：吡喃酮/吡喃并稠环此前落到饱和
    # oxane/oxano 模板，环内 C=C 被静默丢掉（pyran-2-one → oxan-2-one，oxano[3,2-c]pyridine → 缺双键）。
    # 必须排在 oxane 之前：氢化骨架兜底按 _Q_H 表序遍历，先命中者定母体词干。
    "pyran":       {"smiles": "O1C=CC=CC1", "stem_en": "pyran",       "stem_zh": "吡喃",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrano", "吡喃并")},
    "imidazole":   {"smiles": "c1cnc[nH]1", "stem_en": "imidazole",   "stem_zh": "咪唑",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("imidazo", "咪唑并"), "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyrazole":    {"smiles": "c1ccn[nH]1", "stem_en": "pyrazole",    "stem_zh": "吡唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    # 1,3-二唑编号 IUPAC 固定（N 得 1,3 位）；带取代基/H 的 N 走 numbering_engine._narrow_hetero_ring 的 (c) 同元素 N 收窄。
    "oxazole":     {"smiles": "c1cocn1",    "stem_en": "oxazole",     "stem_zh": "噁唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5"), (2, 3, 4, 0, 1))},
    "thiazole":    {"smiles": "c1cscn1",    "stem_en": "thiazole",    "stem_zh": "噻唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5"), (2, 3, 4, 0, 1))},
    # 其余保留名杂芳环（P-25.1 表 2.8）：异噁唑/三唑/四唑/三嗪原缺失
    "isoxazole":   {"smiles": "c1ccno1",    "stem_en": "1,2-oxazole",  "stem_zh": "1,2-噁唑", "naming_class": "monohetero", "fused": True,"locant_prefix": "1,2-"},
    "triazole":    {"smiles": "c1nc[nH]n1", "stem_en": "1,2,4-triazole", "stem_zh": "1,2,4-三唑", "naming_class": "monohetero","fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "tetrazole":   {"smiles": "c1nnn[nH]1", "stem_en": "1H-tetrazole", "stem_zh": "1H-四唑", "naming_class": "monohetero","fused": True,},
    "triazine":    {"smiles": "c1ncncn1",   "stem_en": "1,3,5-triazine", "stem_zh": "1,3,5-三嗪", "naming_class": "monohetero", "fused": True,"locant_prefix": "1,3,5-"},
    # saturated monohetero rings（P-22.2.2；radical/取代基须识别为环而非开链）
    "pyrrolidine": {"smiles": "C1CCNC1",  "stem_en": "pyrrolidine", "stem_zh": "吡咯烷", "naming_class": "monohetero", "fused": False},
    "piperidine":  {"smiles": "C1CCNCC1", "stem_en": "piperidine",  "stem_zh": "哌啶",   "naming_class": "monohetero", "fused": False},
    "morpholine":  {"smiles": "C1COCCN1", "stem_en": "morpholine",  "stem_zh": "吗啉",   "naming_class": "monohetero", "fused": False},
    "piperazine":  {"smiles": "C1CNCCN1", "stem_en": "piperazine",  "stem_zh": "哌嗪",   "naming_class": "monohetero", "fused": False},
    "oxolane":     {"smiles": "C1CCOC1",  "stem_en": "oxolane",     "stem_zh": "四氢呋喃", "naming_class": "monohetero", "fused": False},
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
    # fused56（9 原子）标准编号 1,2,3,3a,4,5,6,7,7a：杂原子(1)走远离桥头方向，苯环从 3a 起。
    "indole":         {"smiles": "c1ccc2[nH]ccc2c1", "stem_en": "1H-indole",      "stem_zh": "吲哚",     "naming_class": "fused56", "fused": True, "fused_stem": ("indole", "吲哚"), "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "indazole":       {"smiles": "c1ccc2cn[nH]c2c1", "stem_en": "indazole",       "stem_zh": "吲唑",     "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (FUSED56_LABELS, (6, 5, 4, 3, 2, 1, 0, 8, 7))},
    "benzimidazole":  {"smiles": "c1ccc2[nH]cnc2c1", "stem_en": "benzimidazole",  "stem_zh": "苯并咪唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "benzofuran":     {"smiles": "c1ccc2occc2c1",    "stem_en": "benzofuran",     "stem_zh": "苯并呋喃", "naming_class": "fused56", "fused": True, "locant_prefix": "1-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzothiophene": {"smiles": "c1ccc2sccc2c1",    "stem_en": "benzothiophene", "stem_zh": "苯并噻吩", "naming_class": "fused56", "fused": True, "locant_prefix": "1-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzothiazole":  {"smiles": "c1ccc2scnc2c1",    "stem_en": "benzothiazole",  "stem_zh": "苯并噻唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzoxazole":    {"smiles": "c1ccc2ocnc2c1",    "stem_en": "benzoxazole",    "stem_zh": "苯并噁唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    # 三环 6+5+6（13 原子）：两个苯环融合吡咯，N9 邻位桥头 4a/9a。
    "carbazole":      {"smiles": "c1ccc2[nH]c3ccccc3c2c1", "stem_en": "carbazole", "stem_zh": "咔唑", "naming_class": "carbazole", "fused": True, "locant_prefix": "9H-", "prefix_nh_conditional": True, "standard": (CARBAZOLE_LABELS, (6, 7, 8, 9, 5, 12, 0, 1, 2, 11, 4, 3, 10))},
    # 三环 6+6+6（14 原子）：中间吡啶/含 S 环两侧苯环融合。
    "acridine":       {"smiles": "c1ccc2nc3ccccc3cc2c1",   "stem_en": "acridine",     "stem_zh": "吖啶",   "naming_class": "acridine", "fused": True, "standard": (ACRIDINE_LABELS, (9, 8, 7, 6, 5, 2, 1, 0, 13, 12, 11, 10, 4, 3))},
    "phenothiazine":  {"smiles": "c1ccc2Sc3ccccc3Nc2c1",   "stem_en": "phenothiazine", "stem_zh": "吩噻嗪", "naming_class": "phenothiazine", "fused": True, "locant_prefix": "10H-", "prefix_nh_conditional": True, "standard": (PHENOTHIAZINE_LABELS, (13, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 5))},
    # 双环 5+6（9 原子）：苯环并二氧戊环，O1/C2/O3。
    "benzodioxole":   {"smiles": "c1ccc2OCOc2c1",          "stem_en": "benzodioxole", "stem_zh": "苯并二氧杂环戊烯", "naming_class": "benzodioxole", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (6, 5, 4, 3, 2, 1, 0, 8, 7))},
    # fused 6+6；naph_family（10 原子）标准编号 1,2,3,4,4a,5,6,7,8,8a。
    "quinoline":    {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline",    "stem_zh": "喹啉",   "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
    "isoquinoline": {"smiles": "c1nccc2ccccc21", "stem_en": "isoquinoline", "stem_zh": "异喹啉", "naming_class": "naph_family", "fused": True},
    "quinazoline":  {"smiles": "c1ccc2ncncc2c1", "stem_en": "quinazoline",  "stem_zh": "喹唑啉", "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
    "quinoxaline":  {"smiles": "c1ccc2nccnc2c1", "stem_en": "quinoxaline",  "stem_zh": "喹喔啉", "naming_class": "naph_family", "fused": True},
    # 表 2.8 第 8 项：噌啉（PIN cinnoline，1,2-二氮杂萘）——N1 邻桥头，N2 次邻；
    # 母体名已固定 N1/N2 位次，故登记 standard（同喹啉 order：N1 起沿环经 4a 绕外周）。
    "cinnoline":    {"smiles": "c1ccc2nnccc2c1", "stem_en": "cinnoline",    "stem_zh": "噌啉", "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},

    # 苯并吡喃（10 原子 6+6）：O 直接连桥头（色烯）或隔一位（异色烯）。保留名 chromene/isochromene
    # 是 P-25.1 表 2.8 的稠合母体，取代「benzo[b]oxane」拼接式（自造体例，两份基准 gold 均 0 见）。
    # 二者连通性不同（O 是否连桥头），氢化骨架匹配路径可区分，故可同表并存。
    # chromene（2H-色烯）：O 连桥头，编号从 O1 起；isochromene（1H-异色烯）：编号自 CH2 起，O 得 2 位。
    "chromene":     {"smiles": "C1=COc2ccccc2C1", "stem_en": "chromene",   "stem_zh": "色烯",  "naming_class": "naph_family", "fused": True, "fused_prefix": ("chromeno", "色烯并"), "standard": (NAPH_LABELS, (2, 1, 0, 9, 8, 7, 6, 5, 4, 3))},
    "isochromene":  {"smiles": "C1=CC2=CC=CC=C2CO1", "stem_en": "isochromene", "stem_zh": "异色烯", "naming_class": "naph_family", "fused": True, "fused_prefix": ("isochromeno", "异色烯并"), "standard": (NAPH_LABELS, (8, 9, 0, 1, 2, 3, 4, 5, 6, 7))},
    # 保留名（表 2.8）：purine=嘌呤（5+6，特殊编号 1–9，PIN 7H-purine，P-25 表 2.8 第16项）；
    # pteridine=蝶啶（6+6 四 N，naph-family 编号，N 在 1/3/5/8、CH 在 2/4/6/7）。无环外 =O，故可入表。
    # purine 9 原子：模板 c1ncc2[nH]cnc2n1 原子按 locant 1–9 序（N1、C2、N3、C4、C5、C6、N7、C8、N9），
    # 桥头 C4/C5 得数字位，N7 为指示氢所在；pteridine 10 原子沿外周按 1,2,3,4,4a,5…8,8a（N 在 1/3/5/8）。
    "purine":       {"smiles": "c1ncc2[nH]cnc2n1", "stem_en": "7H-purine",     "stem_zh": "嘌呤",   "naming_class": "purine", "fused": True, "fused_stem": ("purine", "嘌呤"), "locant_prefix": "7H-", "prefix_nh_conditional": True, "standard": (PURINE_LABELS, (1, 0, 8, 7, 3, 2, 4, 5, 6))},
    "pteridine":    {"smiles": "c1cnc2ncncc2n1",    "stem_en": "pteridine",    "stem_zh": "蝶啶",   "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
    # 呫吨/噻吨（表 2.8 第 22 项；P-25.3.3 传统编号）：外周两苯环 1-8，中央碳 9（CH2/C=O）、O(或 S) 10，
    # 四个稠合碳 4a/8a/9a/10a（走行方向同 acridine）。未登记时 xanthone 被拆成 benzo[b]chromen-13-one
    # 之类 chromene 体系（该体系无 13 位），并出现 benzo[b]benzo[b]thian 前缀重复。
    "xanthene":     {"smiles": "C1c2ccccc2Oc2ccccc21", "stem_en": "xanthene",     "stem_zh": "氧杂蒽", "naming_class": "xanthene", "fused": True, "standard": (XANTHENE_LABELS, (12, 11, 10, 9, 8, 5, 4, 3, 2, 1, 0, 13, 7, 6))},
    "thioxanthene": {"smiles": "C1c2ccccc2Sc2ccccc21", "stem_en": "thioxanthene", "stem_zh": "噻吨",   "naming_class": "xanthene", "fused": True, "standard": (XANTHENE_LABELS, (12, 11, 10, 9, 8, 5, 4, 3, 2, 1, 0, 13, 7, 6))},
    # cyclopenta[a]phenanthrene（P-25.3.3 传统甾体编号 1-17，无字母位；10/13 为角甲基碳）。
    # 母体模板为 7 对非累积双键的 mancude 型（2=3/4=5/6=7/8=9/11=12/14=15/16=17），1 位为 CH2。
    # 甾体及其加氢衍生物经氢化骨架匹配（_fixed_numbering 的 _Q_H 回退）走此模板，否则外周编号
    # 会产生 4a,6a-dimethyl、-2-en-2-yl 等非甾体定位。
    "cyclopenta[a]phenanthrene": {"smiles": "C1=CCC2C(=C1)C=CC1=C2C=CC2C=CC=C12", "stem_en": "cyclopenta[a]phenanthrene", "stem_zh": "环戊[a]菲", "naming_class": "steroid", "fused": True, "standard": (STEROID_LABELS, (2, 1, 0, 5, 4, 6, 7, 8, 9, 3, 10, 11, 12, 16, 15, 14, 13))},
    # NOTE: carbonyl mothers（benzoquinone / anthraquinone / chromenone /
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
    """附加组分的保留稠合前缀 (en, zh)（P-25.3.2.2.3 保留前缀 / P-25.3.2.2.1 单环烃）；无登记则 None（调用方走通用规则）。"""
    entry = _TEMPLATES.get(sid)
    if entry is not None:
        return entry.get("fused_prefix")
    return fusion_carbocycle_prefix(sid)


def fusion_carbocycle_prefix(sid: str) -> tuple[str, str] | None:
    """单环烃附加组分前缀 (en, zh)（P-25.3.2.2.1）；非该类组分返回 None。"""
    entry = _FUSION_CARBOCYCLES.get(sid)
    return (entry["prefix_en"], entry["prefix_zh"]) if entry else None


def omits_fusion_numbers(sid: str) -> bool:
    """稠合描述符是否省略数字位次：一级单环烃附加组分 benzo 及 P-25.3.2.2.1 组分（P-25.3.8.1）。"""
    return sid == "benzene" or sid in _FUSION_CARBOCYCLES


def match_fusion_carbocycle(info: dict, atom_ids) -> str | None:
    """单环烃骨架精确等于某 P-25.3.2.2.1 附加组分时返回 sid，否则 None（元素签名预过滤 + 骨架子图同构）。"""
    mol = info["mol"]
    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)
    for sid, q in _Q_CYCLO.items():
        if _CYCLO_ELEM[sid] != elem:
            continue
        for m in mol.GetSubstructMatches(q, uniquify=True):
            if set(m) == atoms:
                return sid
    return None


def match_fusion_component(info: dict, atom_ids) -> str | None:
    """稠环拆解的组分匹配（P-25.3.2）：保留 mancude 母体优先，其次单环烃附加组分（P-25.3.2.2.1）。"""
    return match_retained(info, atom_ids, mancude_only=True) or match_fusion_carbocycle(info, atom_ids)


# 查询子结构与元素签名，import 时构建一次。
_Q: dict[str, Mol] = {sid: MolFromSmiles(entry["smiles"]) for sid, entry in _TEMPLATES.items()}


# 固定编号视图（由 _TEMPLATES 条目的 `standard` 字段派生；L4 按 scaffold_id 查）
_STANDARD_LABELS: dict[str, tuple[str, ...]] = {
    sid: e["standard"][0] for sid, e in _TEMPLATES.items() if "standard" in e
}
_STANDARD_ORDERS: dict[str, tuple[int, ...]] = {
    sid: e["standard"][1] for sid, e in _TEMPLATES.items() if "standard" in e
}


def _validate_standard_fields() -> None:
    """import 期校验 `standard` 字段：order 须为模板原子的排列，labels 长度须等于模板原子数。"""
    for sid, entry in _TEMPLATES.items():
        std = entry.get("standard")
        if not std:
            continue
        labels, order = std
        n = _Q[sid].GetNumAtoms()
        if sorted(order) != list(range(n)):
            raise ValueError(f"{sid}: standard order {order} 不是 0..{n - 1} 的排列")
        if len(labels) != n:
            raise ValueError(f"{sid}: standard labels 长度 {len(labels)} != 模板原子数 {n}")


_validate_standard_fields()


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


def mancude_atoms(scaffold_id: str, match) -> frozenset[int]:
    """保留模板的 mancude（Kekulé 双键）位映射到分子后的原子集：该集合内部的 C=C 由母体氢化物名隐含（P-31.1.2），不得再写成 -ene/-yne。"""
    if not match or scaffold_id not in _Q:
        return frozenset()
    return frozenset(match[qi] for qi in _kekule_double_atoms(scaffold_id) if qi < len(match))


def extra_indicated_atoms(mol: Mol, scaffold_id: str, match) -> frozenset[int]:
    """保留母体名未隐含、而分子中该芳香杂环位带 H 的原子（P-58.2.1 须显式标指示氢）：模板同位无 H 而分子有 H，
    如 1H-喹啉-4-酮的 N1、1H-嘧啶-2,4-二酮的 N1/N3。全芳香环系（无饱和位）的指示氢只可能来自此处。"""
    q = _Q.get(scaffold_id)
    # 仅稠合母体（≥2 环）：单环 mancude 杂芳环（吡啶/嘧啶）的 [nH] 输入是内酰胺-内酰亚胺互变异构写法，
    # 其位次由母体名与后缀共同固定，gold 不标指示氢；稠合母体（喹啉/异喹啉）无 =N-H 位，不标则名不可解。
    if q is None or not match or len(match) != q.GetNumAtoms() or len(q.GetRingInfo().AtomRings()) < 2:
        return frozenset()
    return frozenset(
        mi for qi, mi in enumerate(match)
        if mi < mol.GetNumAtoms()
        and (qa := q.GetAtomWithIdx(qi)).GetAtomicNum() != 6
        and qa.GetAtomicNum() == mol.GetAtomWithIdx(mi).GetAtomicNum()
        and qa.GetTotalNumHs() == 0
        and mol.GetAtomWithIdx(mi).GetTotalNumHs() > 0
        and mol.GetAtomWithIdx(mi).GetIsAromatic()
    )


def hydrogenated_atoms(mol: Mol, scaffold_id: str, match) -> frozenset[int]:
    """match（模板原子→分子原子）下被加氢的分子原子集：模板某原子承载不饱和双键（Kekulé）而分子中该位已非芳香者（P-31.2.2：hydro 修饰源于双键的饱和）。用原子芳香性而非键级，桥头位归属不受取代基影响。"""
    if not match or mol is None or scaffold_id not in _Q:
        return frozenset()
    from namepredict.layer4.hydrogenation import HYDRO_MULT_N

    ring_atoms = set(match)
    out: set[int] = set()
    suffix: set[int] = set()  # 环内碳带环外多重键（=O/=N 后缀位）：本身不计入加氢，但其配对位要靠指示氢收尾
    for qi in _kekule_double_atoms(scaffold_id):
        if qi >= len(match):
            continue
        mi = match[qi]
        if mi >= mol.GetNumAtoms():
            continue
        atom = mol.GetAtomWithIdx(mi)
        if all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()) and atom.GetTotalNumHs() > 0:
            out.add(mi)  # 分子中该位已无多重键（原带双键、现饱和）才是加氢位；残留芳香/多重键者未加氢（hybridization 对 NH 会误报 SP2，故查键级）；不带 H 的位（季碳、4,4-二甲基型）加不了 H，不占 hydro 位次，交给指示氢
        elif atom.GetAtomicNum() == 6 and any(
                b.GetBondType() != Chem.BondType.SINGLE and b.GetOtherAtomIdx(mi) not in ring_atoms
                for b in atom.GetBonds()):
            suffix.add(mi)
    # 环杂原子（N/O/S）失去双键后新增的 H 由指示氢承载（P-58.2.1），不计入 hydro 计数：
    # 计入会得到「2,3-dihydro-1H-喹啉」的杂原子位而被写成「1,2,3-trihydro」，且 3 为奇数使 hydro_prefix
    # 整体放弃（奇数不在倍增表内），连正确的 2,3-dihydro 一起丢。仅当剔除后计数合法（偶数倍增）才剔除，
    # 否则保留原集合（如 1,2-二氢吡啶：N1+C2 恰为 2，两者同为 hydro 位）。
    carbons = frozenset(a for a in out if mol.GetAtomWithIdx(a).GetAtomicNum() == 6)
    hydro = carbons if len(carbons) != len(out) and len(carbons) in HYDRO_MULT_N else frozenset(out)
    # 计数仍为奇数：环内带后缀 =O/=N 的位（其 H 被后缀取代）使配对加氢位多出一个，该位改由指示氢承载
    # （P-58.2.1），naphthalen-1-one 遂得「2H」+「3,4-dihydro」而非整体放弃。
    if len(hydro) not in HYDRO_MULT_N and suffix:
        for mi in sorted(hydro):
            if any(mol.GetBondBetweenAtoms(mi, s) is not None for s in suffix):
                hydro = frozenset(hydro - {mi})
                break
    return hydro


def _elem_sig(mol: Mol, atom_ids) -> frozenset:
    """原子集的元素组成签名（(Z, 计数) 冻结集合）。"""
    return frozenset(Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids).items())


_TEMPLATE_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q.items()}

# 单环烃附加组分（P-25.3.2.2.1）：饱和单环烃名删尾 'ne' 得前缀，表示最大数目非累积双键的形式。
# 只作稠合附加零件，故不入 _TEMPLATES：入表会让单环骨架解析成保留名，破坏 P-31 单环通用路径
# （carbocycle 按环大小动态命名）；它们也不是母体组分（P-25.3.2.1.1：单环烃母体用 [n]annulene/苯）。
_FUSION_CARBOCYCLES: dict[str, dict] = {
    "cyclopropane": {"smiles": "C1CC1",     "prefix_en": "cyclopropa", "prefix_zh": "环丙并"},
    "cyclobutane":  {"smiles": "C1CCC1",    "prefix_en": "cyclobuta",  "prefix_zh": "环丁并"},
    "cyclopentane": {"smiles": "C1CCCC1",   "prefix_en": "cyclopenta", "prefix_zh": "环戊并"},
    "cyclohexane":  {"smiles": "C1CCCCC1",  "prefix_en": "cyclohexa",  "prefix_zh": "环己并"},
    "cycloheptane": {"smiles": "C1CCCCCC1", "prefix_en": "cyclohepta", "prefix_zh": "环庚并"},
    "cyclooctane":  {"smiles": "C1CCCCCCC1","prefix_en": "cycloocta",  "prefix_zh": "环辛并"},
}


def _cyclo_component_query(smiles: str) -> Mol:
    """由环状 SMILES 的原子数派生纯碳环骨架查询（SMARTS 缺省键 = 单键或芳香键）。"""
    n = MolFromSmiles(smiles).GetNumAtoms()
    return MolFromSmarts("[#6]1" + "[#6]" * (n - 1) + "1")


_Q_CYCLO: dict[str, Mol] = {sid: _cyclo_component_query(e["smiles"]) for sid, e in _FUSION_CARBOCYCLES.items()}
_CYCLO_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q_CYCLO.items()}


def _spec_from_template(sid: str, entry: dict) -> ScaffoldSpec:
    """由模板条目派生 ScaffoldSpec（n_rings/ring 从 smiles 自动算）。"""
    q = _Q[sid]
    n_rings = len(q.GetRingInfo().AtomRings())
    ring = "carbo" if all(q.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in range(q.GetNumAtoms())) else "hetero"
    standard = entry.get("standard")
    numbering = NumberingPolicy(standard_path=standard[0] if standard else ())
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


def match_retained(info: dict, atom_ids, *, mancude_only: bool = False) -> str | None:
    """返回模板精确覆盖 atom_ids 的保留母体 sid，无命中 None；元素签名预过滤后子图同构，同命中取表序第一个（防御性兜底）；mancude_only 只认 fused 保留名（融合组分）。"""
    hit = _match_with_map(info, atom_ids, mancude_only=mancude_only)
    return hit[0] if hit else None


def _match_with_map(info: dict, atom_ids, *, mancude_only: bool = False) -> tuple[str, tuple[int, ...]] | None:
    """模板精确覆盖 atom_ids 时返回 (sid, match)；match[i] 供 standard_path 把模板原子映射到分子原子。 """
    mol = info["mol"]
  
    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)  
    for sid, q in _Q.items():
        if mancude_only and not _TEMPLATES[sid].get("fused"):
            continue  # 饱和保留名不作稠合组分（P-25.2.1 表 2.8）
        if _TEMPLATE_ELEM[sid] != elem:
            continue
        for m in mol.GetSubstructMatches(q, uniquify=True):  # print(sid,m)
            if set(m) == atoms:  # print("yes")
                return sid, m
    mol_h = memo.by_mol("hydrogenated", _hydrogenated, mol)  # 精确匹配失败后按完全氢化骨架再比对（加氢衍生物，P-25.3.4）；整分子重建对同一 mol 只做一次
    if mol_h is not None:
        for sid, qh in _Q_H.items():
            if mancude_only and not _TEMPLATES[sid].get("fused"):
                continue  # 氢化骨架同样只取 mancude 母体（P-25.3.4）
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
        n_rings = sum(1 for ring in sssr_rings(mol) if set(ring) <= atoms)
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
