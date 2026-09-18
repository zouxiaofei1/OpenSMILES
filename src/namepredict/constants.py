"""跨层共享常量：原子序数、倍数前缀与名称文本规范化。"""

from __future__ import annotations

# ── 原子序数 ──────────────────────────
H  = 1
Li = 3
B  = 5
C  = 6
N  = 7
O  = 8
F  = 9
Na = 11
Al = 13
Si = 14
P  = 15
S  = 16
Cl = 17
K  = 19
Ga = 31
Ge = 32
As = 33
Se = 34
Br = 35
In = 49
Sn = 50
Sb = 51
Te = 52
I  = 53
Tl = 81
Pb = 82
Bi = 83

# ── 常用集合 ──────────────────────────
HALO_Z = frozenset({F, Cl, Br, I})
RING_HETERO = frozenset({N, O, S})  # 环内杂原子：其单碳酰基按环酮命名（内酰胺/内酯/硫代内酯、N-酰基环胺），见 layer1.analyzer._is_ketone_carbon。
HALO_ZH = {F: "氟", Cl: "氯", Br: "溴", I: "碘"}
HALIDE_EN = {F: "fluoride", Cl: "chloride", Br: "bromide", I: "iodide"}
N_PREFIX_KINDS = frozenset({"n_alkyl", "n_block"})  # N-取代基 kind（P-62.2.2.1）：走 N- 前缀、位次以 N 标注或隐含省略，不参与数字位次通道。

P25_SENIOR = (N, F, Cl, Br, I, O, S, Se, Te, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)  # 杂原子优先序（稠环母体组分选择 P-25.3.2.4 / 稠环与杂环编号 P-25.3.3.1.2(b)）两条序列同源不同序，勿混用：P25 用于"选哪个组分当母体"，P145 用于"哪个杂原子得低位次"。
P145_SENIOR = (F, Cl, Br, I, O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)

_MULT_EN_10 = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",  # 倍数前缀（P-14.2 Table 1.4；数量词与母链碳数同源，支持到 99）en 数量词（deca/undeca/…/icosa/henicosa…）同时供 layer5.stems 生成长链母链词干，故下沉到本层：en_num_term 是唯一来源，stems 取词干仅去其尾 'a'。
               6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
_MULT_ZH_10 = {1: "", 2: "二", 3: "三", 4: "四", 5: "五",
               6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}
_EN_UNIT = {1: "hen", 2: "do", 3: "tri", 4: "tetra", 5: "penta",
            6: "hexa", 7: "hepta", 8: "octa", 9: "nona"}
_EN_DECADE = {1: "deca", 2: "icosa", 3: "triaconta", 4: "tetraconta",
              5: "pentaconta", 6: "hexaconta", 7: "heptaconta",
              8: "octaconta", 9: "nonaconta"}
ZH_DIGITS = "一二三四五六七八九"


def en_num_term(n: int) -> str | None:
    """英文数值词干（末带 'a'）：1–99 按个位+十位组合生成。"""
    if n < 1 or n > 99:
        return None
    if n <= 10:
        return _MULT_EN_10[n]
    if n == 11:
        return "undeca"
    tens, ones = divmod(n, 10)
    if tens == 1:  # 12–19：个位 + deca
        return f"{_EN_UNIT[ones]}deca"
    dec = _EN_DECADE[tens]
    if ones == 0:
        return dec
    unit = _EN_UNIT[ones]
    if tens == 2 and unit[-1] in "aeiou":  # 20s 前导 i 省略：docosa/tricosa…
        dec = "cosa"
    return f"{unit}{dec}"


def zh_numeral(n: int) -> str | None:
    """中文数值词（倍数用，1 返空）：≤10 查表，11–99 按十一/二十二组合。"""
    if n < 1 or n > 99:
        return None
    if n <= 10:
        return _MULT_ZH_10[n]
    tens, ones = divmod(n, 10)
    head = "十" if tens == 1 else f"{ZH_DIGITS[tens-1]}十"
    return head if ones == 0 else f"{head}{ZH_DIGITS[ones-1]}"


MULT_EN = {n: (en_num_term(n) or "") for n in range(1, 100)}
MULT_ZH = {n: (zh_numeral(n) or "") for n in range(1, 100)}

def zh_bridge_root(name: str) -> str:
    """桥后缀（氨基/氧基/硫基）前的中文烃基名去尾「基」。"""
    return name[:-1] if name.endswith("基") else name

AMIDO_RETAINED = {  # P-66.1.1.4.3
    "acetyl": ("acetamido", "乙酰氨基"),
    "formyl": ("formamido", "甲酰胺基"),
    "benzoyl": ("benzamido", "苯甲酰胺基"),
}
AMIDO_RETAINED_EN = frozenset(v[0] for v in AMIDO_RETAINED.values())

# ── L0 电荷归一 / 盐解离 ────────────────────
DONOR_KIND = ("carboxyl", "phospho")  # 允许作为强酸供体的酸类：羧酸 + 磷酸。磷酸供体在「酰胺 O⁻ 受体」场景下才有产出（见 preprocessor 的二次互变归一，gold 把 N=C([O-]) 写成酰胺、把 P-OH 写成 oxidophosphoryl）；sulfo 实测 0 收益，关闭以免扩大 blast radius。
ACCEPTOR_Z = frozenset({O, N})  # 弱受体允许的元素：O（酚氧/烯醇氧/酰胺氧）、N（去质子化氮）；保守可只留 {O}。
ACID_CENTERS = {          # 中心元素 → 酸类名；键序即供体搬运的酸强度序
    C: "carboxyl",        # C(=O)OH
    P: "phospho",         # P(=O)OH
    S: "sulfo",           # S(=O)nOH
}
ALKALI_EN = {Li: "lithium", Na: "sodium", K: "potassium"}  # 原子序数 → 英文金属名（IUPAC 官能团类盐）
METAL_ZH = {"lithium": "锂", "sodium": "钠", "potassium": "钾"}

# ── L2/L5 单核母体氢化物（P-15.4.1）────────────
MONONUCLEAR_HYDRIDES: dict[str, tuple] = {
    "oxidane":    ("O", "氧化烷", ("hydroxy",   "羟基"),     ("oxy",      "氧基"), "氧基"),
    "azane":      ("N", "氮烷",   ("amino",     "氨基"),     ("amino",    "氨基"), "氨基"),
    "sulfane":    ("S", "硫烷",   ("sulfanyl",  "硫基"),     ("sulfanyl", "硫基"), "硫基"),
    "sulfinyl":   ("S", "亚磺酰", ("sulfinyl",  "亚磺酰基"), None,                 "基亚磺酰基"),
    "sulfonyl":   ("S", "磺酰",   ("sulfonyl",  "磺酰基"),   None,                 "磺酰基"),
    "imine":      ("N", "亚胺",   ("imino",     "亚氨基"),   None,                 "亚氨基"),
    "phosphoryl": ("P", "磷酰",   ("phosphoryl","磷酰基"),   None,                 "磷酰基"),
    "phosphanyl": ("P", "磷烷基", ("phosphanyl","磷烷基"),   None,                 "磷烷基"),
}
MONONUCLEAR_BY_ELEMENT = {N: "azane", O: "oxidane", P: "phosphoryl", S: "sulfane"}  # 原子序数 → 该元素的默认 free 名（表 2.1 每元素一行），锚点自由基据此起步
SULFUR_STEM_BY_OXO = {0: "sulfane", 1: "sulfinyl", 2: "sulfonyl"}  # P-63.2.2：S 锚点的 =O 数 → free 名，=O 数必须落进母体名，否则 S(=O)/S(=O)(=O) 与硫醚同形（氧被整段丢弃、净多 2H）。
PHOSPHORUS_STEM_BY_OXO = {0: "phosphanyl", 1: "phosphoryl"}  # P-67.1.4.1.1.2/6：P 锚点的 =O 数 → free 名，否则 P(=O) 与膦同形。
NITROGEN_STEM_BY_FREE_DOUBLE = {False: "azane", True: "imine"}  # P-66.1.1：N 锚点的自由价键级 → free 名，漏掉双键会把亚胺写成胺（净多 2H）。
MONONUCLEAR_ZERO_YL = {(en, v[1]): v[2] for en, v in MONONUCLEAR_HYDRIDES.items()}         # (free_en, free_zh) → 零价去氢名
MONONUCLEAR_BRIDGE = {(en, v[1]): v[3] for en, v in MONONUCLEAR_HYDRIDES.items() if v[3]}  # (free_en, free_zh) → 桥后缀（仅 O/N/S 三行）
MONONUCLEAR_YL = {en: (v[1], (v[3] or v[2])[0], v[4]) for en, v in MONONUCLEAR_HYDRIDES.items()}  # free_en → (free_zh, 去氢 yl_en, 组装名中文尾)
PHOSPHORYL_STEMS = tuple(en for en, v in MONONUCLEAR_HYDRIDES.items() if v[0] == "P")  # 替代碳词干的 P 酰基词干（避免 P 被当碳中心）

# ── L3 取代基词表 ──────────────────────
SIMPLE_ALKOXY_NO_PAREN = frozenset({  # 简单保留烷氧基作前缀不加括号
    "methoxy", "ethoxy", "propoxy", "butoxy", "phenoxy", "isopropoxy",
})
NAME_KIND = {  # 锚定表/保留叶子的名称暗含非烷基 kind，使 L5 以 `kind` 键（iso 稠合、聚茴香醚、卤代苯）触发。
    "fluoro": "halo", "chloro": "halo", "bromo": "halo", "iodo": "halo",
   
   
}
CLAIM_KIND = {"amine_n": "n_block",
              "ring_c": "alkyl", "chain_c": "alkyl"}  # claim 槽位 → 取代基 kind
# O-侧臂 kind：连在 parent 的 O 上的侧链作 O 侧烷基，由 L5 酯/含氧酸整名消费
ESTER_O_SIDE_KINDS = frozenset({"ester", "phosphate", "phosphonate", "sulfate", "sulfonate"})
# 中心原子自任母体的含氧酸 kind（L2/L5 共用）；碳锚定的磺酸族不在内
OXO_CENTER_KINDS = frozenset({"phosphate", "phosphonate", "sulfate"})

# ── L4 位次 / 编号 ──────────────────────
HYDRO_MULT_N = frozenset({2, 4, 6, 8, 10, 12, 14, 16, 18, 20})  # 加氢前缀覆盖的氢原子数（P-31.2.2 以偶数倍增前缀表示双键的饱和，位次数为加氢原子数）；数量词本身取自本层（唯一来源），此处只表达 L4 的适用域，域外放弃而非给错名。
TRADITIONAL_NUMBERING_IDS = frozenset({  # P-25.3.3：这些保留骨架按传统编号，不走 P-25.3.3.1 的优选取向自动编号。xanthene 及其硫属类似物（xanthene/thioxanthene）与 cyclopenta[a]phenanthrene（甾体 1-17）已按_TEMPLATES 的 standard 字段登记传统编号，故一并列入。
    "anthracene", "phenanthrene", "acridine", "carbazole", "purine",
    "xanthene", "thioxanthene", "cyclopenta[a]phenanthrene",
})
RS_HI = frozenset({"R", "M", "r"})   # 编号优先级高的 CIP 描述符（P-91.2）
RS_LO = frozenset({"S", "P", "s"})   # 与之成对、取较低位次的次位描述符

# ── L5 组装词表 ──────────────────────
BIS_EN = {2: "bis", 3: "tris", 4: "tetrakis"}  # P-16.3.2 复合前缀倍增（bis/tris，非 di/tri）
BIS_ZH = {2: "双", 3: "三", 4: "四"}
BRIDGE_SUFFIX_EN = ("oxy", "sulfanyl", "amino")   # O/S/N 桥后缀（P-63.2.2.1）：平铺 [-yl]oxy 与合一 [-yloxy] 均属合法
BRIDGE_SUFFIX_ZH = ("氧基", "硫基", "氨基")        # 与 BRIDGE_SUFFIX_EN 同序同位
DIATOMIC_BRIDGE_YL = ("diazenyl", "disulfanyl")   # 双原子桥合一保留前缀：-N=N-R / -S-S-R 收成 R-diazenyl / R-disulfanyl
SIMPLE_BRIDGE_YL_NO_PAREN = frozenset(f"phenyl{s}" for s in DIATOMIC_BRIDGE_YL)  # 裸苯基前端 + 双原子桥：作前缀免括
BRIDGE_FUSION_YL: dict[tuple[str, str], tuple[tuple[str, ...], str, str]] = {
    # (中心单核氢化物词干, 前端名尾 en) → (前端名尾 zh 候选, 合一前缀 en, 合一前缀 zh)
    ("azane", "imino"): (("亚氨基",), "diazenyl", "二氮烯基"),            # P-68.3.1.3：diazenyl 名优先于 azo 名
    ("sulfane", "sulfanyl"): (("硫基",), "disulfanyl", "二硫代基"),  # P-68.3.1.4：-S-S- 的保留前缀
}
BRIDGE_DIATOMIC_ZH = tuple(v[2] for v in BRIDGE_FUSION_YL.values())  # 双原子桥合一前缀的中文名尾
# 可拆桥后缀（含双原子桥）：前端复合时前端加围栏、桥留括号外
BRIDGE_SPLIT_SUFFIX_EN = BRIDGE_SUFFIX_EN + DIATOMIC_BRIDGE_YL
BRIDGE_SPLIT_SUFFIX_ZH = BRIDGE_SUFFIX_ZH + BRIDGE_DIATOMIC_ZH
BRIDGE_YL_SUFFIX = tuple((f"yl{en}", en) for en in BRIDGE_SUFFIX_EN)      # -yl 型桥基：方括号闭在 -yl 后、桥后缀放括号外（[(2R)环己基]氧基）
BRIDGE_ZH_YL_SUFFIX = tuple((f"基{zh}", zh) for zh in BRIDGE_SUFFIX_ZH)
EXO_RING_SUF: dict[str, tuple] = {  # 环外主基系统名后缀表（group_class → 后缀规格），六个环外 worker 共用 `_exocyclic_ring_names` 一条管线：singular (en, zh) 单取代后缀；plural (en, zh) | None 多取代后缀基底（前拼 MULT_EN/MULT_ZH 倍数词）；None = 该主基无多取代系统名
    "acid":     (("carboxylic acid", "羧酸"), ("carboxylic acid", "羧酸")),
    "aldehyde": (("carbaldehyde", "甲醛"),    ("carbaldehyde", "甲醛")),
    "ester":    (("carboxylate", "羧酸"),     ("carboxylate", "羧酸")),
    "amide":    (("carboxamide", "甲酰胺"),   ("carboxamide", "甲酰胺")),
    "nitrile":  (("carbonitrile", "甲腈"),    ("carbonitrile", "甲腈")),
    "acyl":     (("carbonyl", "羰基"),        ("carbonyl", "羰基")),
    "acyl_halide": (("carbonyl", "甲酰"),     None),  # 卤素词由 chain_engine 按 hal_z 动态拼接
}
AZANE_PAREN_SUF = ("benzoyl", "carbonyl", "acetyl")  # azane 单取代基内层组加括号的尾缀白名单
ALKOXY_YLOXY_EN = (  
  ("enyloxy", "enoxy"),  ("ethyloxy", "ethoxy"), ("propyloxy", "propoxy"), ("butyloxy", "butoxy"),
    ("phenyloxy", "phenoxy"),
)
ALKOXY_YLOXY_ZH = (
    ("乙基氧基", "乙氧基"), ("丙基氧基", "丙氧基"), ("丁基氧基", "丁氧基"),
    ("苯基氧基", "苯氧基"),
)
CHAIN_RETAINED = {  # C1/C2 开链英文 IUPAC 保留名（formic/acetic…）；C3+ 系统名由词干生成（_chain_plain 回落），中文无保留名。
    "acid": {1: ("formic acid", "甲酸"), 2: ("acetic acid", "乙酸")},
    "acyl": {1: ("formyl", "甲酰基"), 2: ("acetyl", "乙酰基")},
    "aldehyde": {1: ("formaldehyde", "甲醛"), 2: ("acetaldehyde", "乙醛")},
    "amide": {1: ("formamide", "甲酰胺"), 2: ("acetamide", "乙酰胺")},
    "nitrile": {1: ("formonitrile", "甲腈"), 2: ("acetonitrile", "乙腈")},
    "ester": {1: ("formate", "甲酸"), 2: ("acetate", "乙酸")},
}
RETAINED_FUSION_ALIASES: dict[str, tuple[str, str]] = {  # 稠合组装名 → 保留名（P-25.1.1）；只收位次形态与稠合名完全相同的条目，命中即整名替换
    "benzo[c]furan":       ("2-benzofuran",  "2-苯并呋喃"),    # 异苯并呋喃：O 占 2 位，1/3 位为环碳
    "benzo[c]pyrrole":     ("isoindole",     "异吲哚"),        # N 占 2 位
    "benzo[b]benzofuran":  ("dibenzofuran",  "二苯并呋喃"),      # 兼容无方括号写法
    "benzo[b][1]benzofuran": ("dibenzofuran", "二苯并呋喃"),     # 组分名按 P-25.3.5 引用 [1]
    "benzo[b]quinoxaline": ("phenazine",     "菲嗪"),
    "benzo[a]indene":      ("fluorene",      "芴"),
    "benzo[d]1,2-oxazole": ("1,2-benzoxazole", "1,2-苯并噁唑"),
    "benzo[b]anthracene":  ("tetracene",     "并四苯"),
}
HS_NUMBER = "甲乙丙丁戊己庚辛壬癸"
