"""骨架规格 + 保留 SMILES 模板表 + 环解析（唯一事实来源）。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol, MolFromSmarts, MolFromSmiles
from rdkit import Chem
from namepredict.constants import HW_COMPONENT_PREFIX, C
from namepredict.tools import memo
from namepredict.layer2.parent_skeleton import ParentSkeleton, _ring_count
from namepredict.layer1.ring_systems import kekulized


def _ring_kind(mol, atom_ids) -> str:
    """环组分的元素类型：全碳为 carbo，否则 hetero。"""
    return "carbo" if all(mol.GetAtomWithIdx(a).GetAtomicNum() == C for a in atom_ids) else "hetero"


@dataclass(frozen=True)# ScaffoldSpec 定义
class NumberingPolicy:
    """骨架编号策略：固定编号路径。"""
    standard_path: tuple = ()


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
    locant_prefix: str = ""  # 五元杂环 locant 前缀（P-61.2.4 / P-25.1）
    prefix_nh_conditional: bool = False

    @property
    def identity(self) -> ScaffoldIdentity:
        """返回本规格的 ScaffoldIdentity。"""
        return ScaffoldIdentity(self.id, self.naming_class, self.n_rings, self.ring)


@dataclass(frozen=True)
class ScaffoldIdentity:
    """骨架身份：id、命名类、环数与环类型。"""
    id: str
    naming_class: str
    n_rings: int
    ring: str


FUSED56_LABELS: tuple[str, ...] = ("1", "2", "3", "3a", "4", "5", "6", "7", "7a")  # 保留 fused 母体的固定编号标签（P-25.4）
PURINE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "5", "6", "7", "8", "9")  # purine（嘌呤）保留传统编号：桥头得纯数字位（P-25.3.3）
CARBAZOLE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a", "9", "9a", "9b")  # carbazole（13 原子）：N9，桥头位带 a（P-25.4.1.4）
ACRIDINE_LABELS: tuple[str, ...] = ("1","2","3","4","4a","5","6","7","8","8a","9","9a","10","10a")  # acridine（14 原子）：N10，对位 C9 连取代基
PHENOTHIAZINE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "9", "9a", "10", "10a", "10b")  # phenothiazine（14 原子）：S5、N10
NAPH_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a")
ANTHRACENE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "10", "10a", "5", "6", "7", "8", "8a", "9", "9a")  # anthracene(14 原子): 中环碳得数字位(P-25.4.1)
PHENANTHRENE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "4b", "5", "6", "7", "8", "8a", "9", "10", "10a")  # phenanthrene(14 原子) 传统编号(P-14.4(a))：中环两个 CH 得 9,10，桥头 4a,4b,8a,10a
PYRENE_LABELS: tuple[str, ...] = ("1", "2", "3", "3a", "4", "5", "5a", "6", "7", "8", "8a", "8b", "9", "10", "10a", "10b")  # pyrene(16 原子): 外周 1-10(P-25.3.3.3.1)
XANTHENE_LABELS: tuple[str, ...] = ("1", "2", "3", "4", "4a", "5", "6", "7", "8", "8a", "9", "9a", "10", "10a")  # xanthene/thioxanthene: 中央碳 9、O/S 10
STEROID_LABELS: tuple[str, ...] = tuple(str(i) for i in range(1, 18))  # 甾体传统编号 1-17 全数字(10/13 为角甲基碳)
ADAMANTANE_LABELS: tuple[str, ...] = tuple(str(i) for i in range(1, 11))  # 金刚烷保留编号 1-10：1,3,5,7 为次甲基碳，其余为亚甲基碳

_TEMPLATES: dict[str, dict] = {  # 保留母体 SMILES 模板注册表（唯一事实来源）
    "benzene":     {"smiles": "c1ccccc1",             "stem_en": "benzene",    "stem_zh": "苯",   "naming_class": "mono_carbo", "fused": True, "fused_prefix": ("benzo", "苯并")},  # carbocycles
    "naphthalene": {"smiles": "c1ccc2ccccc2c1",       "stem_en": "naphthalene","stem_zh": "萘",    "naming_class": "naph_family", "fused": True, "fused_prefix": ("naphtho", "萘并")},
    "anthracene":  {"smiles": "c1ccc2cc3ccccc3cc2c1", "stem_en": "anthracene", "stem_zh": "蒽",    "naming_class": "anthra", "fused": True, "fused_prefix": ("anthra", "蒽并"), "standard": (ANTHRACENE_LABELS, (13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))},  # 模板原子序即外周环序，缺此字段蒽位号会错
    "phenanthrene":{"smiles": "c1ccc2c(c1)ccc1ccccc12", "stem_en": "phenanthrene","stem_zh": "菲", "naming_class": "phenanthrene", "fused": True, "fused_prefix": ("phenanthro", "菲并"), "standard": (PHENANTHRENE_LABELS, (9, 10, 11, 12, 13, 3, 2, 1, 0, 5, 4, 6, 7, 8))},  # 模板原子序即按 1,2,3,4,4a,4b,5,6,7,8,8a,9,10,10a 的外周行走序
    "pyrene":      {"smiles": "c1cc2ccc3cccc4ccc(c1)c2c34", "stem_en": "pyrene",  "stem_zh": "芘", "naming_class": "pyrene", "fused": True, "standard": (PYRENE_LABELS, (6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 14, 3, 4, 5, 15))},
    "indene":      {"smiles": "C1=CCc2ccccc21", "stem_en": "1H-indene", "stem_zh": "1H-茚", "naming_class": "fused56", "fused": True, "fused_stem": ("indene", "茚"), "standard": (FUSED56_LABELS, (2, 1, 0, 8, 7, 6, 5, 4, 3))},  # 茚（PIN 1H-indene），5+6 稠合碳环，并入 fused56
    "chrysene":    {"smiles": "c1ccc2c(c1)ccc1c3ccccc3ccc21", "stem_en": "chrysene", "stem_zh": "屈", "naming_class": "chrysene", "fused": True, "fused_prefix": ("chryseno", "䓛并")},  # 䓛（PIN chrysene，中文用「屈」），四环稠烃
    "picene":      {"smiles": "c1ccc2c(c1)ccc1c2ccc2c3ccccc3ccc21", "stem_en": "picene", "stem_zh": "苉", "naming_class": "picene", "fused": True, "fused_prefix": ("piceno", "苉并")},  # 苉（表 2.7 第 6 位保留名），五环稠烃，按 P-25.3.3 编号
    "furan":       {"smiles": "c1ccoc1",    "stem_en": "furan",       "stem_zh": "呋喃",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("furo", "呋喃并")},  # monocyclic heteroarenes
    "thiophene":   {"smiles": "c1ccsc1",    "stem_en": "thiophene",   "stem_zh": "噻吩",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("thieno", "噻吩并")},
    "pyrrole":     {"smiles": "c1cc[nH]c1", "stem_en": "pyrrole",     "stem_zh": "吡咯",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyridine":    {"smiles": "n1ccccc1",   "stem_en": "pyridine",    "stem_zh": "吡啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrido", "吡啶并")},
    "pyridazine":  {"smiles": "c1ccnnc1",   "stem_en": "pyridazine",  "stem_zh": "哒嗪",   "naming_class": "monohetero", "fused": True},
    "pyrimidine":  {"smiles": "c1cncnc1",   "stem_en": "pyrimidine",  "stem_zh": "嘧啶",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrimido", "嘧啶并")},
    "pyrazine":    {"smiles": "c1cnccn1",   "stem_en": "pyrazine",    "stem_zh": "吡嗪",   "naming_class": "monohetero", "fused": True},
    "pyran":       {"smiles": "O1C=CC=CC1", "stem_en": "pyran",       "stem_zh": "吡喃",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrano", "吡喃并")},  # 吡喃（P-25.1 表 2.8 保留名，6 元含氧 mancude 母体）
    "imidazole":   {"smiles": "c1cnc[nH]1", "stem_en": "imidazole",   "stem_zh": "咪唑",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("imidazo", "咪唑并"), "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "pyrazole":    {"smiles": "c1ccn[nH]1", "stem_en": "pyrazole",    "stem_zh": "吡唑",   "naming_class": "monohetero", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "oxazole":     {"smiles": "c1cocn1",    "stem_en": "oxazole",     "stem_zh": "噁唑",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3]oxazolo", "[1,3]噁唑并"), "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5"), (2, 3, 4, 0, 1))},  # 1,3-二唑编号固定（N 得 1,3 位）
    "thiazole":    {"smiles": "c1cscn1",    "stem_en": "thiazole",    "stem_zh": "噻唑",   "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3]thiazolo", "[1,3]噻唑并"), "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5"), (2, 3, 4, 0, 1))},
    "isoxazole":   {"smiles": "c1ccno1",    "stem_en": "1,2-oxazole",  "stem_zh": "1,2-噁唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2]oxazolo", "[1,2]噁唑并"), "locant_prefix": "1,2-"},  # 其余保留名杂芳环（P-25.1 表 2.8，前缀带方括号）
    "triazole":    {"smiles": "c1nc[nH]n1", "stem_en": "1,2,4-triazole", "stem_zh": "1,2,4-三唑", "naming_class": "monohetero","fused": True, "fused_prefix": ("[1,2,4]triazolo", "[1,2,4]三唑并"), "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "tetrazole":   {"smiles": "c1nnn[nH]1", "stem_en": "tetrazole", "stem_zh": "四唑", "naming_class": "monohetero","fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},  # 指示氢条件化：环含未取代 NH 时才注入 1H-
    "triazine":    {"smiles": "c1ncncn1",   "stem_en": "1,3,5-triazine", "stem_zh": "1,3,5-三嗪", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3,5]triazino", "[1,3,5]三嗪并"), "locant_prefix": "1,3,5-"},
    "pyrrolidine": {"smiles": "C1CCNC1",  "stem_en": "pyrrolidine", "stem_zh": "吡咯烷", "naming_class": "monohetero", "fused": False},  # saturated monohetero rings（P-22.2.2）
    "piperidine":  {"smiles": "C1CCNCC1", "stem_en": "piperidine",  "stem_zh": "哌啶",   "naming_class": "monohetero", "fused": False},
    "morpholine":  {"smiles": "C1COCCN1", "stem_en": "morpholine",  "stem_zh": "吗啉",   "naming_class": "monohetero", "fused": False},
    "thiomorpholine": {"smiles": "C1CSCCN1", "stem_en": "thiomorpholine", "stem_zh": "硫代吗啉", "naming_class": "monohetero", "fused": False},  # 吗啉硫类似物（表 2.8 保留名），未登记则回退成 1,4-thiazinane
    "piperazine":  {"smiles": "C1CNCCN1", "stem_en": "piperazine",  "stem_zh": "哌嗪",   "naming_class": "monohetero", "fused": False},
    "oxolane":     {"smiles": "C1CCOC1",  "stem_en": "oxolane",     "stem_zh": "四氢呋喃", "naming_class": "monohetero", "fused": False},
    "oxane":       {"smiles": "C1CCCOC1", "stem_en": "oxane",       "stem_zh": "氧杂环己烷", "naming_class": "monohetero", "fused": True, "fused_prefix": ("pyrano", "吡喃并")},  # 稠合组分须 mancude（P-25.3.1.2.1/.2.3）：由 oxane 派生应取 pyrano
    "oxirane":     {"smiles": "C1CO1",    "stem_en": "oxirane",     "stem_zh": "环氧乙烷", "naming_class": "monohetero","fused": True,},  # 小环与含硫饱和杂环（P-22.2.2）
    "thiolane":    {"smiles": "C1CCSC1",  "stem_en": "thiolane",    "stem_zh": "四氢噻吩", "naming_class": "monohetero","fused": True, "fused_prefix": ("thieno", "噻吩并")},  # 同上：thiolane → thieno
    "thiane":      {"smiles": "C1CCSCC1", "stem_en": "thiane",      "stem_zh": "四氢噻喃", "naming_class": "monohetero","fused": True, "fused_prefix": ("thiopyrano", "噻喃并")},  # 同上：thiane → thiopyrano（噻喃）
    "dioxolane":   {"smiles": "C1COCO1",  "stem_en": "1,3-dioxolane", "stem_zh": "1,3-二氧戊环", "naming_class": "monohetero", "locant_prefix": "1,3-"},  # 双氧/三氧饱和环（缩醛/缩酮、溶剂类骨架）
    "dioxane":     {"smiles": "C1COCCO1", "stem_en": "1,4-dioxane",   "stem_zh": "1,4-二氧六环", "naming_class": "monohetero", "locant_prefix": "1,4-"},
    "trioxane":    {"smiles": "C1OCOCO1", "stem_en": "1,3,5-trioxane", "stem_zh": "1,3,5-三氧六环", "naming_class": "monohetero", "locant_prefix": "1,3,5-"},
    "dithiole":    {"smiles": "S1SC=CC1", "stem_en": "dithiole", "stem_zh": "二硫杂环戊二烯", "naming_class": "monohetero", "fused": False},  # 1,2-二硫杂环戊二烯 mancude 母体；须先于 dithiolane 登记
    "dithiolane12": {"smiles": "C1CSSC1", "stem_en": "dithiolane", "stem_zh": "二硫杂环戊烷", "naming_class": "monohetero", "fused": False},  # 1,2-二硫戊环（P-22.2.2 HW 名；金标不写 1,2- 位次）
    "imidazolidine":{"smiles": "C1NCCN1", "stem_en": "imidazolidine",  "stem_zh": "咪唑烷", "naming_class": "monohetero"},
    "pyrazolidine": {"smiles": "C1CNNC1", "stem_en": "pyrazolidine", "stem_zh": "吡唑烷", "naming_class": "monohetero"},  # 表 2.3 保留名：饱和吡唑环须用它而非氢化 pyrazole
    "dihydrofuran":  {"smiles": "C1C=CCO1",   "stem_en": "2,5-dihydrofuran", "stem_zh": "2,5-二氢呋喃", "naming_class": "monohetero", "standard": (("1", "2", "3", "4", "5"), (4, 0, 1, 2, 3))},  # 部分不饱和 5/6 元杂环（P-22.2.2 加氢前缀）；字面位次即固定编号(P-14.4(a)/(b))
    "dihydropyran":  {"smiles": "C1=COCCC1",  "stem_en": "3,4-dihydro-2H-pyran", "stem_zh": "3,4-二氢-2H-吡喃", "naming_class": "monohetero", "standard": (("1", "2", "3", "4", "5", "6"), (2, 3, 4, 5, 0, 1))},
    "dihydropyrrole":{"smiles": "C1C=CCN1",   "stem_en": "2,5-dihydro-1H-pyrrole", "stem_zh": "2,5-二氢-1H-吡咯", "naming_class": "monohetero", "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (("1", "2", "3", "4", "5"), (4, 0, 1, 2, 3))},
    "dihydroimidazole":{"smiles": "C1=NCCN1", "stem_en": "4,5-dihydro-1H-imidazole", "stem_zh": "4,5-二氢-1H-咪唑", "naming_class": "monohetero", "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (("1", "2", "3", "4", "5"), (4, 0, 1, 2, 3))},
    "dihydrothiazole":{"smiles": "C1=NCCS1",  "stem_en": "4,5-dihydro-1,3-thiazole", "stem_zh": "4,5-二氢-1,3-噻唑", "naming_class": "monohetero", "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5"), (4, 0, 1, 2, 3))},
    "indole":         {"smiles": "c1ccc2[nH]ccc2c1", "stem_en": "1H-indole",      "stem_zh": "吲哚",     "naming_class": "fused56", "fused": True, "fused_stem": ("indole", "吲哚"), "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},  # fused 5+6（9 原子）标准编号，杂原子走远离桥头方向
    "indazole":       {"smiles": "c1ccc2cn[nH]c2c1", "stem_en": "indazole",       "stem_zh": "吲唑",     "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (FUSED56_LABELS, (6, 5, 4, 3, 2, 1, 0, 8, 7))},
    "benzimidazole":  {"smiles": "c1ccc2[nH]cnc2c1", "stem_en": "benzimidazole",  "stem_zh": "苯并咪唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True},
    "benzofuran":     {"smiles": "c1ccc2occc2c1",    "stem_en": "benzofuran",     "stem_zh": "苯并呋喃", "naming_class": "fused56", "fused": True, "locant_prefix": "1-", "fused_prefix": ("[1]benzofuro", "[1]苯并呋喃并"), "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},  # 附加组分前缀：benzofuran → [1]benzofuro（保留前缀，非「去尾 e 加 o」通用式）
    "benzofuran2":    {"smiles": "c1ccc2cocc2c1",    "stem_en": "benzofuran",     "stem_zh": "苯并呋喃", "naming_class": "fused56", "fused": True, "locant_prefix": "2-", "fused_prefix": ("[2]benzofuro", "[2]苯并呋喃并"), "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},  # 异苯并呋喃（PIN 2-benzofuran，P-25.1 表 2.8）：O 居五元环中央不与桥头相邻，与 1- 异构区分
    "benzothiophene": {"smiles": "c1ccc2sccc2c1",    "stem_en": "benzothiophene", "stem_zh": "苯并噻吩", "naming_class": "fused56", "fused": True, "locant_prefix": "1-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzothiazole":  {"smiles": "c1ccc2scnc2c1",    "stem_en": "benzothiazole",  "stem_zh": "苯并噻唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzoxazole":    {"smiles": "c1ccc2ocnc2c1",    "stem_en": "benzoxazole",    "stem_zh": "苯并噁唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},
    "benzotriazole":  {"smiles": "c1ccc2[nH]nnc2c1",    "stem_en": "benzotriazole",  "stem_zh": "苯并三唑", "naming_class": "fused56", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True, "standard": (FUSED56_LABELS, (4, 5, 6, 7, 8, 0, 1, 2, 3))},  # 苯并三唑（P-25.1 表 2.8）：N1/N2/N3，桥头 3a/7a
    "carbazole":      {"smiles": "c1ccc2[nH]c3ccccc3c2c1", "stem_en": "carbazole", "stem_zh": "咔唑", "naming_class": "carbazole", "fused": True, "locant_prefix": "9H-", "prefix_nh_conditional": True, "standard": (CARBAZOLE_LABELS, (6, 7, 8, 9, 5, 12, 0, 1, 2, 11, 4, 3, 10))},  # 三环 6+5+6（13 原子）：两个苯环融合吡咯，N9 邻位桥头 4a/9a。
    "acridine":       {"smiles": "c1ccc2nc3ccccc3cc2c1",   "stem_en": "acridine",     "stem_zh": "吖啶",   "naming_class": "acridine", "fused": True, "standard": (ACRIDINE_LABELS, (9, 8, 7, 6, 5, 2, 1, 0, 13, 12, 11, 10, 4, 3))},  # 三环 6+6+6（14 原子）：中间吡啶/含 S 环两侧苯环融合。
    "phenothiazine":  {"smiles": "c1ccc2Sc3ccccc3Nc2c1",   "stem_en": "phenothiazine", "stem_zh": "吩噻嗪", "naming_class": "phenothiazine", "fused": True, "locant_prefix": "10H-", "prefix_nh_conditional": True, "standard": (PHENOTHIAZINE_LABELS, (13, 0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 5))},
    "benzodioxole":   {"smiles": "c1ccc2OCOc2c1",          "stem_en": "benzodioxole", "stem_zh": "苯并二氧杂环戊烯", "naming_class": "benzodioxole", "fused": True, "locant_prefix": "1,3-", "standard": (FUSED56_LABELS, (6, 5, 4, 3, 2, 1, 0, 8, 7))},  # 双环 5+6（9 原子）：苯环并二氧戊环，O1/C2/O3。
    "quinoline":    {"smiles": "c1ccc2ncccc2c1", "stem_en": "quinoline",    "stem_zh": "喹啉",   "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},  # fused 6+6；naph_family（10 原子）标准编号
    "isoquinoline": {"smiles": "c1nccc2ccccc21", "stem_en": "isoquinoline", "stem_zh": "异喹啉", "naming_class": "naph_family", "fused": True},
    "quinazoline":  {"smiles": "c1ccc2ncncc2c1", "stem_en": "quinazoline",  "stem_zh": "喹唑啉", "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
    "quinoxaline":  {"smiles": "c1ccc2nccnc2c1", "stem_en": "quinoxaline",  "stem_zh": "喹喔啉", "naming_class": "naph_family", "fused": True},
    "cinnoline":    {"smiles": "c1ccc2nnccc2c1", "stem_en": "cinnoline",    "stem_zh": "噌啉", "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},  # 噌啉（PIN cinnoline），母体名已固定位次
    "dioxine": {"smiles": "C1=COC=CO1", "stem_en": "1,4-dioxine", "stem_zh": "1,4-二噁英", "naming_class": "monohetero", "fused": True,   "fused_prefix": ("[1,4]dioxino", "[1,4]二噁英并"), "locant_prefix": "1,4-",  "standard": (("1","2","3","4","5","6"), (2,3,4,5,0,1))},
    "chromene":     {"smiles": "C1=COc2ccccc2C1", "stem_en": "chromene",   "stem_zh": "色烯",  "naming_class": "naph_family", "fused": True, "fused_prefix": ("chromeno", "色烯并"), "standard": (NAPH_LABELS, (2, 1, 0, 9, 8, 7, 6, 5, 4, 3))},  # 苯并吡喃（保留名 chromene/isochromene）
    "isochromene":  {"smiles": "C1=CC2=CC=CC=C2CO1", "stem_en": "isochromene", "stem_zh": "异色烯", "naming_class": "naph_family", "fused": True, "fused_prefix": ("isochromeno", "异色烯并"), "standard": (NAPH_LABELS, (8, 9, 0, 1, 2, 3, 4, 5, 6, 7))},
    "purine":       {"smiles": "c1ncc2[nH]cnc2n1", "stem_en": "7H-purine",     "stem_zh": "嘌呤",   "naming_class": "purine", "fused": True, "fused_stem": ("purine", "嘌呤"), "locant_prefix": "7H-", "prefix_nh_conditional": True, "standard": (PURINE_LABELS, (1, 0, 8, 7, 3, 2, 4, 5, 6))},  # 保留名 purine（嘌呤）/ pteridine（蝶啶）
    "indolizine":   {"smiles": "c1ccn2ccccc12", "stem_en": "indolizine", "stem_zh": "中氮茚", "naming_class": "indolizine", "fused": True, "fused_prefix": ("indolizino", "中氮茚并")},  # 表 2.8 第 20 位保留名，5+6 稠环，编号按 P-25.3.3（N 得数字位 4）
    "pyrrolizine":  {"smiles": "C1=CC2=CC=CN2C1", "stem_en": "pyrrolizine", "stem_zh": "吡咯嗪", "naming_class": "pyrrolizine", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True, "fused_prefix": ("pyrrolizino", "吡咯嗪并")},  # 表 2.8 第 21 位保留名 1H-pyrrolizine，5+5 稠环，按 P-25.3.3 编号
    "pteridine":    {"smiles": "c1cnc2ncncc2n1",    "stem_en": "pteridine",    "stem_zh": "蝶啶",   "naming_class": "naph_family", "fused": True, "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},
    "xanthene":     {"smiles": "C1c2ccccc2Oc2ccccc21", "stem_en": "xanthene",     "stem_zh": "氧杂蒽", "naming_class": "xanthene", "fused": True, "standard": (XANTHENE_LABELS, (12, 11, 10, 9, 8, 5, 4, 3, 2, 1, 0, 13, 7, 6))},  # 呫吨/噻吨（表 2.8 第 22 项；P-25.3.3 传统编号）
    "thioxanthene": {"smiles": "C1c2ccccc2Sc2ccccc21", "stem_en": "thioxanthene", "stem_zh": "噻吨",   "naming_class": "xanthene", "fused": True, "standard": (XANTHENE_LABELS, (12, 11, 10, 9, 8, 5, 4, 3, 2, 1, 0, 13, 7, 6))},
    "cyclopenta[a]phenanthrene": {"smiles": "C1=CCC2C(=C1)C=CC1=C2C=CC2C=CC=C12", "stem_en": "cyclopenta[a]phenanthrene", "stem_zh": "环戊[a]菲", "naming_class": "steroid", "fused": True, "standard": (STEROID_LABELS, (2, 1, 0, 5, 4, 6, 7, 8, 9, 3, 10, 11, 12, 16, 15, 14, 13))},  # 环戊[a]菲（P-25.3.3 传统甾体编号 1-17）
     "oxadiazole124": {"smiles": "o1ncnc1", "stem_en": "1,2,4-oxadiazole", "stem_zh": "1,2,4-噁二唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2,4]oxadiazolo", "[1,2,4]噁二唑并"), "locant_prefix": "1,2,4-", "standard": (("1", "2", "3", "4", "5"), (0, 1, 2, 3, 4))},
    "oxadiazole134": {"smiles": "o1cnnc1", "stem_en": "1,3,4-oxadiazole", "stem_zh": "1,3,4-噁二唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3,4]oxadiazolo", "[1,3,4]噁二唑并"), "locant_prefix": "1,3,4-", "standard": (("1", "2", "3", "4", "5"), (0, 1, 2, 3, 4))},
    "oxadiazole125": {"smiles": "o1nccn1", "stem_en": "1,2,5-oxadiazole", "stem_zh": "1,2,5-噁二唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2,5]oxadiazolo", "[1,2,5]噁二唑并"), "locant_prefix": "1,2,5-", "standard": (("1", "2", "3", "4", "5"), (0, 1, 2, 3, 4))},
    "thiadiazole134": {"smiles": "s1cnnc1", "stem_en": "1,3,4-thiadiazole", "stem_zh": "1,3,4-噻二唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3,4]thiadiazolo", "[1,3,4]噻二唑并"), "locant_prefix": "1,3,4-", "standard": (("1", "2", "3", "4", "5"), (0, 1, 2, 3, 4))},
    "thiadiazole124": {"smiles": "s1ncnc1", "stem_en": "1,2,4-thiadiazole", "stem_zh": "1,2,4-噻二唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2,4]thiadiazolo", "[1,2,4]噻二唑并"), "locant_prefix": "1,2,4-", "standard": (("1", "2", "3", "4", "5"), (0, 1, 2, 3, 4))},
    "triazine124": {"smiles": "n1ncncc1", "stem_en": "1,2,4-triazine", "stem_zh": "1,2,4-三嗪", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2,4]triazino", "[1,2,4]三嗪并"), "locant_prefix": "1,2,4-", "standard": (("1", "2", "3", "4", "5", "6"), (0, 1, 2, 3, 4, 5))},
    "tetrazine1245": {"smiles": "n1ncnnc1", "stem_en": "1,2,4,5-tetrazine", "stem_zh": "1,2,4,5-四嗪", "naming_class": "monohetero", "fused": True, "locant_prefix": "1,2,4,5-", "standard": (("1", "2", "3", "4", "5", "6"), (0, 1, 2, 3, 4, 5))},
    "thiazole12": {"smiles": "c1cnsc1", "stem_en": "1,2-thiazole", "stem_zh": "1,2-噻唑", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,2]thiazolo", "[1,2]噻唑并"), "locant_prefix": "1,2-", "standard": (("1", "2", "3", "4", "5"), (3, 2, 1, 0, 4))},  # 异噻唑：S1/N2 相邻（原 c1cncs1 与 thiazole 同构，1,2- 名实不符）
    "azepane": {"smiles": "N1CCCCCC1", "stem_en": "azepane", "stem_zh": "氮杂环庚烷", "naming_class": "monohetero", "fused": True, "fused_stem": ("azepine", "氮杂卓"), "locant_prefix": "1H-", "prefix_nh_conditional": True},  # 作稠合母体须取 mancude 词干 azepine（P-25.3.1.2.2）
    "thiazine13": {"smiles": "S1C=NC=CC1", "stem_en": "1,3-thiazine", "stem_zh": "1,3-噻嗪", "naming_class": "monohetero", "fused": True, "fused_prefix": ("[1,3]thiazino", "[1,3]噻嗪并"), "locant_prefix": "1,3-", "standard": (("1", "2", "3", "4", "5", "6"), (0, 1, 2, 3, 4, 5))},
    "azulene":      {"smiles": "c1ccc2cccc2cc1", "stem_en": "azulene", "stem_zh": "薁", "naming_class": "naph_family", "fused": True, "fused_prefix": ("azuleno", "薁并")},  # 表 2.7 保留名，5+7 稠合碳环（桥头 3a/8a）
    "pentalene":    {"smiles": "C1=CC=C2C=CC=C12", "stem_en": "pentalene", "stem_zh": "戊搭烯", "naming_class": "pentalene", "fused": True, "fused_prefix": ("pentaleno", "戊搭烯并"), "standard": (("1", "2", "3", "3a", "4", "5", "6", "6a"), tuple(range(8)))},  # 表 2.7 保留名 5+5 稠合双环（P-25.1.2.3 多轮烯），桥头 3a/6a
    "phenalene":    {"smiles": "C1ccc2cccc3cccc1c23", "stem_en": "phenalene", "stem_zh": "菲那烯", "naming_class": "phenalene", "standard": (("1", "2", "3", "3a", "4", "5", "6", "6a", "7", "8", "9", "9a", "9b"), tuple(range(13)))},  # 表 2.7 保留名 1H-phenalene（P-25.3.3.3.1 新编号）：内碳 9b，模板原子序即外周 1,2,3,3a…9a
    "phthalazine":  {"smiles": "c1ccc2cnncc2c1", "stem_en": "phthalazine", "stem_zh": "酞嗪", "naming_class": "naph_family", "fused": True, "fused_prefix": ("phthalazino", "酞嗪并"), "standard": (NAPH_LABELS, (4, 5, 6, 7, 8, 9, 0, 1, 2, 3))},  # 2,3-二氮杂萘：N 得 2/3 位，C1 与 C4 分居两环
    "naphthyridine18": {"smiles": "c1cc2cccnc2nc1", "stem_en": "1,8-naphthyridine", "stem_zh": "1,8-萘啶", "naming_class": "naph_family", "fused": True, "fused_prefix": ("[1,8]naphthyridino", "[1,8]萘啶并")},  # N1/N8 分居两环
    "benzodioxine14": {"smiles": "O1C=COc2ccccc12", "stem_en": "1,4-benzodioxine", "stem_zh": "1,4-苯并二噁英", "naming_class": "naph_family", "fused": True, "locant_prefix": "1,4-"},  # 2,3-二氢体即常见 1,4-苯并二噁烷母体
    "benzodioxin12": {"smiles": "O1OC=CC2=CC=CC=C12", "stem_en": "1,2-benzodioxin", "stem_zh": "1,2-苯并二噁英", "naming_class": "naph_family", "fused": True, "locant_prefix": "1,2-"},  # 相邻双氧（O1/O2）
    "benzothiophene2": {"smiles": "c1ccc2cscc2c1", "stem_en": "2-benzothiophene", "stem_zh": "2-苯并噻吩", "naming_class": "fused56", "fused": True, "fused_prefix": ("2-benzothieno", "2-苯并噻吩并")},  # 异苯并噻吩（P-25.1 表 2.8）；与 c1ccc2sccc2c1 的 1- 异构区分
    "pyrrolopyridazine12b": {"smiles": "c1ccn2ncccc12", "stem_en": "pyrrolo[1,2-b]pyridazine", "stem_zh": "吡咯并[1,2-b]哒嗪", "naming_class": "fused56", "fused": True, "fused_prefix": ("pyrrolo[1,2-b]pyridazino", "吡咯并[1,2-b]哒嗪并")},  # 桥头 N（5+6，9 原子）
    # 7H-pyrrolo[2,3-d]pyrimidine 未登记：它会与 pyrimido[5,4-b]indole 争夺母体组分并致金标退化（tiers-25570），实测净收益为负
    "pyrazolopyrimidine54d": {"smiles": "c1n[nH]c2ncncc12", "stem_en": "pyrazolo[5,4-d]pyrimidine", "stem_zh": "吡唑并[5,4-d]嘧啶", "naming_class": "purine", "fused": True, "fused_prefix": ("pyrazolo[5,4-d]pyrimidino", "吡唑并[5,4-d]嘧啶并")},
    "adamantane":   {"smiles": "C1C2CC3CC1CC(C2)C3", "stem_en": "adamantane", "stem_zh": "金刚烷", "naming_class": "adamantane", "fused": False, "standard": (ADAMANTANE_LABELS, (1, 2, 3, 4, 5, 6, 7, 8, 0, 9))},  # 表 2.7 保留名：三环桥烃，非稠合零件故 fused=False；固定编号 1,3,5,7 给次甲基碳（P-25.1 保留名表）
    "benzodiazepine14": {"smiles": "N1C=CN=Cc2ccccc12", "stem_en": "1,4-benzodiazepine", "stem_zh": "1,4-苯并二氮杂卓", "naming_class": "naph_family", "fused": True, "locant_prefix": "1,4-", "fused_prefix": ("[1,4]benzodiazepino", "[1,4]苯并二氮杂卓并")},  # 两个 N 相隔 C2/C3（P-25.1 表 2.8）
    "benzoxazine31": {"smiles": "N1COCc2ccccc12", "stem_en": "3,1-benzoxazine", "stem_zh": "3,1-苯并噁嗪", "naming_class": "naph_family", "fused": True, "locant_prefix": "3,1-", "fused_prefix": ("[3,1]benzoxazino", "[3,1]苯并噁嗪并")},  # 母体式：N1/C2/O3/C4，羰基由 FG 后缀补（2,4-二酮）
    "phenoxazine":  {"smiles": "c1ccc2Nc3ccccc3Oc2c1", "stem_en": "phenoxazine", "stem_zh": "吩噁嗪", "naming_class": "phenothiazine", "fused": True, "locant_prefix": "10H-", "prefix_nh_conditional": True, "fused_prefix": ("phenoxazino", "吩噁嗪并")},  # 吩噻嗪的 O 类似物（表 2.8 第 4 位）；N10 得指示氢
    "benzazepine1": {"smiles": "C1=CC=Cc2ccccc2N1", "stem_en": "1-benzazepine", "stem_zh": "1-苯并氮杂卓", "naming_class": "naph_family", "fused": True, "locant_prefix": "1H-", "prefix_nh_conditional": True, "fused_stem": ("[1]benzazepine", "[1]苯并氮杂卓")},  # 稠合时按 P-25.3.5 引用位次 [1]
    "benzothiazepine15": {"smiles": "S1C=CCNc2ccccc12", "stem_en": "1,5-benzothiazepine", "stem_zh": "1,5-苯并硫氮杂卓", "naming_class": "naph_family", "fused": True, "locant_prefix": "1,5-", "fused_prefix": ("[1,5]benzothiazepino", "[1,5]苯并硫氮杂卓并")},
    "benzothiazine13": {"smiles": "S1C=NCc2ccccc12", "stem_en": "1,3-benzothiazine", "stem_zh": "1,3-苯并噻嗪", "naming_class": "naph_family", "fused": True, "locant_prefix": "1,3-", "fused_prefix": ("[1,3]benzothiazino", "[1,3]苯并噻嗪并")},

}


def _bracketed_stem(entry: dict) -> tuple[str, str]:
    """杂原子位次带方括号的稠合组分名（P-25.3.5 稠合母体须引用完整位次）。"""
    lp = entry.get("locant_prefix") or ""
    bare = lp[:-1] if lp.endswith("-") else ""
    if bare and bare.replace(",", "").isdigit() \
            and not entry["stem_en"][:1].isdigit():  # 1-/1,3-/1,2,4- 型；"1H-" 与词干已带位次（1,2-oxazole）不重复加
        br = f"[{bare}]"
        return br + entry["stem_en"], br + entry["stem_zh"]
    return entry["stem_en"], entry["stem_zh"]


def component_stem(sid: str) -> tuple[str, str] | None:
    """稠合组分词干 """
    if sid.startswith(HW_COMPONENT_PREFIX):  # 生成式 HW 组分（P-25.3.2.1.2）
        from namepredict.layer2.hantzsch_widman import component_names
        return component_names(sid)
    entry = _TEMPLATES.get(sid)
    if not entry or not entry.get("fused"):
        return None
    return entry.get("fused_stem") or _bracketed_stem(entry)


def retained_fusion_prefix(sid: str) -> tuple[str, str] | None:
    """附加组分的保留稠合前缀 (en, zh)（P-25.3.2.2.3）。"""
    entry = _TEMPLATES.get(sid)
    if entry is not None:
        return entry.get("fused_prefix")
    return fusion_carbocycle_prefix(sid)


def fusion_carbocycle_prefix(sid: str) -> tuple[str, str] | None:
    """单环烃附加组分前缀 (en, zh)（P-25.3.2.2.1）。"""
    entry = _FUSION_CARBOCYCLES.get(sid)
    return (entry["prefix_en"], entry["prefix_zh"]) if entry else None


def omits_fusion_numbers(sid: str) -> bool:
    """稠合描述符是否省略数字位次（P-25.3.8.1）。"""
    return sid == "benzene" or sid in _FUSION_CARBOCYCLES


def match_fusion_carbocycle(info: dict, atom_ids) -> str | None:
    """单环烃骨架精确等于 P-25.3.2.2.1 附加组分时返回 sid。"""
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
    """稠环拆解的组分匹配（P-25.3.2）：保留母体优先，其次单环烃，最后生成式 HW 名。"""
    mol = info["mol"]
    from namepredict.layer2.hantzsch_widman import component_key

    # 按原子集记忆：对称笼架多个环集张成同一原子集，不记忆会重复全模板扫描
    return memo.by_key("fusion_component", (id(mol), frozenset(atom_ids)),
                       lambda: match_retained(info, atom_ids, mancude_only=True)
                       or match_fusion_carbocycle(info, atom_ids)
                       or component_key(mol, atom_ids), mol)

_Q: dict[str, Mol] = {sid: MolFromSmiles(entry["smiles"]) for sid, entry in _TEMPLATES.items()}  # 查询子结构与元素签名，import 时构建一次。

_STANDARD_LABELS: dict[str, tuple[str, ...]] = {  # 固定编号视图（由 _TEMPLATES 的 standard 字段派生）
    sid: e["standard"][0] for sid, e in _TEMPLATES.items() if "standard" in e
}
_STANDARD_ORDERS: dict[str, tuple[int, ...]] = {
    sid: e["standard"][1] for sid, e in _TEMPLATES.items() if "standard" in e
}


def _hydrogenated(mol: Mol) -> Mol | None:
    """返回完全氢化的分子副本（清芳香性并全键改单键）。"""
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

_Q_H: dict[str, Mol] = {  # 完全氢化模板，供匹配加氢衍生物（P-25.3.4）
    sid: h for sid, q in _Q.items()
    if _TEMPLATES[sid]["naming_class"] != "mono_carbo" and (h := _hydrogenated(q)) is not None
}


_KEKULE_ATOMS: dict[str, frozenset[int]] = {}


def _kekule_double_atoms(sid: str) -> frozenset[int]:
    """模板 Kekulé 双键端点原子集（缓存，防稠合单键误判）。"""
    hit = _KEKULE_ATOMS.get(sid)
    if hit is not None:
        return hit
    q = _Q.get(sid)
    out: set[int] = set()
    m = kekulized(q) if q is not None else None  # 缓存的是原子集合，分子副本不参与缓存
    if m is not None:
        out = {b.GetBeginAtomIdx() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE}
        out |= {b.GetEndAtomIdx() for b in m.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE}
    _KEKULE_ATOMS[sid] = frozenset(out)
    return _KEKULE_ATOMS[sid]


def mancude_ring_atoms(scaffold_id: str, match) -> frozenset[int]:
    """保留 mancude 母体名的**整个不饱和环**映射到分子后的原子集。 """
    q = _Q.get(scaffold_id)
    if not match or q is None:
        return frozenset()
    dbl = _kekule_double_atoms(scaffold_id)
    keep = {i for ring in q.GetRingInfo().AtomRings() if set(ring) & dbl for i in ring}
    return frozenset(match[qi] for qi in keep if qi < len(match))


def extra_hydrogenated_atoms(mol: Mol, scaffold_id: str, match) -> frozenset[int]:
    """氢数多于保留母体模板同位、或模板双键位已饱和的环原子（P-58.2.1）：如 4H-异喹啉-1,3-二酮的 C4。"""
    q = _Q.get(scaffold_id or "")
    if q is None or not match or len(match) != q.GetNumAtoms():
        return frozenset()
    dbl = _kekule_double_atoms(scaffold_id)
    out: set[int] = set()
    for qi, mi in enumerate(match):
        if mi >= mol.GetNumAtoms() or not mol.GetAtomWithIdx(mi).IsInRing():
            continue
        atom = mol.GetAtomWithIdx(mi)
        if q.GetAtomWithIdx(qi).GetTotalNumHs() < atom.GetTotalNumHs():
            out.add(mi)  # 模板同位氢数增加：加氢位
        elif qi in dbl and atom.GetTotalNumHs() > 0 and all(  # 模板双键位已饱和：H 数可能不增（该位另有取代基）
                b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()):
            out.add(mi)
    return frozenset(out)


def extra_indicated_atoms(mol: Mol, scaffold_id: str, match) -> frozenset[int]:
    """模板同位无 H 而分子有 H 的芳香杂环原子（P-58.2.1 指示氢）。"""
    q = _Q.get(scaffold_id)
    if q is None or not match or len(match) != q.GetNumAtoms() or len(q.GetRingInfo().AtomRings()) < 2:  # 仅稠合母体（≥2 环）；单环杂芳环的 [nH] 是互变异构写法
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
    """被加氢的分子原子集：模板双键位在分子中已饱和者（P-31.2.2）。"""
    if not match or mol is None or scaffold_id not in _Q:
        return frozenset()
    from namepredict.constants import HYDRO_MULT_N

    ring_atoms = set(match)
    out: set[int] = set()
    suffix: set[int] = set()  # 环内碳带环外多重键（后缀位），不计入加氢
    for qi in _kekule_double_atoms(scaffold_id):
        if qi >= len(match):
            continue
        mi = match[qi]
        if mi >= mol.GetNumAtoms():
            continue
        atom = mol.GetAtomWithIdx(mi)
        if all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds()) and atom.GetTotalNumHs() > 0:
            out.add(mi)  # 该位无多重键（原双键现饱和）才算加氢位
        elif atom.GetAtomicNum() == 6 and any(
                b.GetBondType() != Chem.BondType.SINGLE and b.GetOtherAtomIdx(mi) not in ring_atoms
                for b in atom.GetBonds()):
            suffix.add(mi)
    hydro = frozenset(out)  # 环杂原子新增的 H 与碳位同为「加氢位」，一并进 hydro 前缀（P-31.2.2）
    if len(hydro) not in HYDRO_MULT_N and suffix:  # 计数仍为奇数时，改由指示氢承载一个位（P-58.2.1）
        for mi in sorted(hydro):
            if any(mol.GetBondBetweenAtoms(mi, s) is not None for s in suffix):
                hydro = frozenset(hydro - {mi})
                break
    return hydro


def _elem_sig(mol: Mol, atom_ids) -> frozenset:
    """原子集的元素组成签名（(Z, 计数) 冻结集合）。"""
    return frozenset(Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids).items())


_TEMPLATE_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q.items()}

# 原模板是否含多重键（P-31.2 已饱和单环须排除芳香/不饱和母体）。两轮匹配共用原模板的键型。
_HAS_MULTI_BOND: dict[str, bool] = {
    sid: any(b.GetBondType() != Chem.BondType.SINGLE for b in q.GetBonds())
    for sid, q in _Q.items()
}

_FUSION_CARBOCYCLES: dict[str, dict] = {  # 单环烃附加组分（P-25.3.2.2.1），只作稠合零件
    "cyclopropane": {"smiles": "C1CC1",     "prefix_en": "cyclopropa", "prefix_zh": "环丙并"},
    "cyclobutane":  {"smiles": "C1CCC1",    "prefix_en": "cyclobuta",  "prefix_zh": "环丁并"},
    "cyclopentane": {"smiles": "C1CCCC1",   "prefix_en": "cyclopenta", "prefix_zh": "环戊并"},
    "cyclohexane":  {"smiles": "C1CCCCC1",  "prefix_en": "benzo",      "prefix_zh": "苯并"},  # P-25.3.2.2.1「除 benzo 外」：六元环的稠合前缀是保留前缀 benzo，非 cyclohexa
    "cycloheptane": {"smiles": "C1CCCCCC1", "prefix_en": "cyclohepta", "prefix_zh": "环庚并"},
    "cyclooctane":  {"smiles": "C1CCCCCCC1","prefix_en": "cycloocta",  "prefix_zh": "环辛并"},
}


def _cyclo_component_query(smiles: str) -> Mol:
    """由环状 SMILES 原子数派生纯碳环骨架查询（SMARTS）。"""
    n = MolFromSmiles(smiles).GetNumAtoms()
    # 键级用 ~ 通配：附加组分只论环大小与元素，不饱和碳环（环戊烯/环己烯等）同样可作稠合零件（P-25.3.2.2.1）
    return MolFromSmarts("[#6]1" + "~[#6]" * (n - 1) + "~1")


_Q_CYCLO: dict[str, Mol] = {sid: _cyclo_component_query(e["smiles"]) for sid, e in _FUSION_CARBOCYCLES.items()}
_CYCLO_ELEM: dict[str, frozenset] = {sid: _elem_sig(q, range(q.GetNumAtoms())) for sid, q in _Q_CYCLO.items()}


def _spec_from_template(sid: str, entry: dict) -> ScaffoldSpec:
    """由模板条目派生 ScaffoldSpec（环数/环型自动算）。"""
    q = _Q[sid]
    n_rings = len(q.GetRingInfo().AtomRings())
    ring = _ring_kind(q, range(q.GetNumAtoms()))
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


def get_spec(spec_id: str) -> ScaffoldSpec | None:
    """按 id 查 ScaffoldSpec（无则 None）。"""
    return _BY_ID.get(spec_id)


def all_specs() -> tuple[ScaffoldSpec, ...]:
    """返回全部 ScaffoldSpec（由 _TEMPLATES 派生）。"""
    return _ALL_SPECS


def numbering_scaffold_facts(spec_id: str | None, atom_count: int) -> dict | None:
    """物化纯母体 facts（标签数与原子数相符时返回）。"""
    spec = get_spec(spec_id or "")
    if spec is None:
        return None
    labels = spec.numbering.standard_path
    return None if not labels or len(labels) != atom_count else {
        "scaffold_id": spec.id, "labels": labels,
    }


def match_retained(info: dict, atom_ids, *, mancude_only: bool = False) -> str | None:
    """返回模板精确覆盖 atom_ids 的保留母体 sid（否则 None）。"""
    hit = _match_with_map(info, atom_ids, mancude_only=mancude_only)
    return hit[0] if hit else None


def _induced_bond_count(m: Mol, atoms: frozenset[int]) -> int:
    """原子集诱导子图的键数（只数两端都在集内的键）。"""
    return sum(1 for b in m.GetBonds()
               if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms)


def _is_induced_match(m: Mol, q: Mol, atoms: frozenset[int]) -> bool:
    """匹配须为诱导覆盖：集内键数等于模板键数。

    子图同构容忍目标多出的键，环闭合键会因此隐形；若不校验，
    环系的真子图模板（如萘嵌入三环骨架）会把多出的环静默丢掉。
    """
    return _induced_bond_count(m, atoms) == q.GetNumBonds()


def _isolated_saturated_ring(mol: Mol, atoms: frozenset[int]) -> bool:
    """原子集是否为不与他环稠合的单环，且 Kekulé 视图环内无重键（RDKit 会误判芳香）。"""
    from namepredict.layer2.hantzsch_widman import isolated_ring
    if not isolated_ring(mol, atoms):
        return False
    kek = kekulized(mol) or mol
    return all(b.GetBondType() == Chem.BondType.SINGLE for b in kek.GetBonds()
               if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms)


def _match_with_map(info: dict, atom_ids, *, mancude_only: bool = False) -> tuple[str, tuple[int, ...]] | None:
    """模板精确覆盖 atom_ids 时返回 (sid, match)。"""
    mol = info["mol"]

    atoms = frozenset(atom_ids)
    elem = _elem_sig(mol, atom_ids)
    sat_ring = _isolated_saturated_ring(mol, atoms)  # 已饱和单环须取饱和母体氢化物名（P-31.2）

    def _scan(base: Mol, qmap) -> tuple[str, tuple[int, ...]] | None:
        """在 base 分子上按 qmap 的顺序找精确覆盖 atoms 的模板。"""
        for sid, q in qmap.items():
            if mancude_only and not _TEMPLATES[sid].get("fused"):
                continue  # 饱和保留名不作稠合组分（P-25.2.1 表 2.8）
            if sat_ring and _HAS_MULTI_BOND[sid]:
                continue
            if _TEMPLATE_ELEM[sid] != elem:
                continue
            for m in base.GetSubstructMatches(q, uniquify=True):
                if set(m) == atoms and _is_induced_match(base, q, atoms):
                    return sid, m
        return None

    hit = _scan(mol, _Q)
    if hit is not None:
        return hit
    mol_h = memo.by_mol("hydrogenated", _hydrogenated, mol)  # 精确匹配失败后按完全氢化骨架再比对（P-25.3.4）
    if mol_h is None:
        return None
    inhit = _scan(mol_h, _Q_H)
    if inhit is None:
        return None
    # 孤立单环在匹配原子集内仍有非单键：按饱和保留名命名会整段丢掉不饱和度，改走生成式 HW
    from namepredict.layer2.hantzsch_widman import isolated_ring
    sid = inhit[0]
    if not _HAS_MULTI_BOND[sid] and isolated_ring(mol, atoms) and any(
            b.GetBondType() != Chem.BondType.SINGLE
            and b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
            for b in mol.GetBonds()):
        return None
    return inhit


def locant_prefix(spec_id: str | None) -> tuple[str, str, bool]:
    """返回 scaffold 的 locant 前缀三元组。"""
    spec = get_spec(spec_id or "")
    if spec is None or not spec.locant_prefix:
        return "", "", False
    return spec.locant_prefix, spec.locant_prefix, spec.prefix_nh_conditional


def standard_chain(spec_id: str | None, match: tuple[int, ...] | None) -> list[int] | None:
    """把 fused 模板固定编号映射到分子（无标准顺序返回 None）。"""
    if not match:
        return None
    order = _STANDARD_ORDERS.get(spec_id or "")
    if not order or len(order) != len(match):
        return None
    return [match[t] for t in order]

def _generic_carbocycle(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    n_rings = _ring_count(info["mol"], skeleton)
    if n_rings >= 2:
        return ScaffoldIdentity("fused_hetero", "fused_hetero", n_rings, "hetero")
    return ScaffoldIdentity("carbocycle", "carbocycle", 1, "carbo")
   


def resolve_ring_scaffold(info: dict, skeleton: ParentSkeleton) -> ScaffoldIdentity | None:
    """解析骨架的 scaffold 身份"""
    sid = match_retained(info, skeleton.atom_ids)  # 按模板子图同构匹配 scaffold id
    if sid:
        spec = get_spec(sid)
        if spec:
            return spec.identity
    from namepredict.layer2.hantzsch_widman import identity as hw_identity
    hw = hw_identity(info, skeleton)  # 未命中模板的 3-10 元杂单环走生成式 HW（P-22.2.2）
    if hw is not None:
        return hw
    # print("no_sid")
    return _generic_carbocycle(info, skeleton)
