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
Mg = 12
Al = 13
Si = 14
P  = 15
S  = 16
Cl = 17
K  = 19
Ca = 20
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
Po = 84
At = 85

# ── 常用集合 ──────────────────────────
HALO_Z = frozenset({F, Cl, Br, I})
HALO_ZH = {F: "氟", Cl: "氯", Br: "溴", I: "碘"}
HALIDE_EN = {F: "fluoride", Cl: "chloride", Br: "bromide", I: "iodide"}
N_PREFIX_KINDS = frozenset({"n_alkyl", "n_block"})  # N-取代基 kind（P-62.2.2.1）：走 N- 前缀、位次以 N 标注或隐含省略，不参与数字位次通道。
N_LOCANT_KINDS = frozenset({"urea", "thiourea", "guanidine"})  # 保留名母体：N 取代基按 1/3 数字位次引用，位次不可省（P-66.3.1 脲/硫脲/胍）

# ── 同位素（P-8 同位素修饰化合物）──────────
ISO_NUCLIDE_PROP = "isoNuclide"  # 非氢核素属性，值形如 "13C"/"18F"/"123I"，挂核素原子自身
ISO_H_PROPS = ((2, "iso2H", "deuterio", "氘代"), (3, "iso3H", "tritio", "氚代"))  # (质量数, 重原子属性名, 英文词干, 中文词干)

P25_SENIOR = (N, F, Cl, Br, I, O, S, Se, Te, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)  # 杂原子优先序（稠环母体组分选择 P-25.3.2.4 / 稠环与杂环编号 P-25.3.3.1.2(b)）两条序列同源不同序，勿混用：P25 用于"选哪个组分当母体"，P145 用于"哪个杂原子得低位次"。
P145_SENIOR = (F, Cl, Br, I, O, S, Se, Te, N, P, As, Sb, Bi, Si, Ge, Sn, Pb, B, Al, Ga, In, Tl)

# ── Hantzsch-Widman 杂单环（P-22.2.2）────────────
HW_ID = "hw_mono"  # 生成式 HW 杂单环骨架 id：不登记进 _TEMPLATES，靠生成器产出词干
HW_CLASS = "heterocycle"  # 其 naming_class（与 monohetero 分列，避免改动既有模板行为）
HW_COMPONENT_PREFIX = "hw:"  # 稠合组分编码前缀：'hw:' + 环序元素符号串

HW_PREFIX_EN = {9: "fluora", 17: "chlora", 35: "broma", 53: "ioda",  # Table 2.4 'a' 前缀（按优先性递减）
                8: "oxa", 16: "thia", 34: "selena", 52: "tellura",
                7: "aza", 15: "phospha", 33: "arsa", 51: "stiba", 83: "bisma",
                14: "sila", 32: "germa", 50: "stanna", 82: "plumba",
                5: "bora", 13: "aluma", 31: "galla", 49: "indiga", 81: "thalla"}
HW_PREFIX_ZH = {9: "氟杂", 17: "氯杂", 35: "溴杂", 53: "碘杂",  # 表 3-3 中文前缀（长式）
                8: "氧杂", 16: "硫杂", 34: "硒杂", 52: "碲杂",
                7: "氮杂", 15: "磷杂", 33: "砷杂", 51: "锑杂", 83: "铋杂",
                14: "硅杂", 32: "锗杂", 50: "锡杂", 82: "铅杂",
                5: "硼杂", 13: "铝杂", 31: "镓杂", 49: "铟杂", 81: "铊杂"}
HW_ZH_SHORT = {(5, 8): "噁", (5, 16): "噻", (6, 16): "噻"}  # 含氮环的短式（噁唑/噻唑/噻嗪），键为 (环大小, 原子序数)
HW_MAX_VALENCE = {1: 1, 5: 3, 6: 4, 7: 3, 8: 2, 9: 1, 13: 3, 14: 4, 15: 3, 16: 2, 17: 1,  # Table 2.4 键数，供 mancude 参照
                  31: 3, 32: 4, 33: 3, 34: 2, 35: 1, 49: 3, 50: 4, 51: 3, 52: 2, 53: 1, 81: 3, 82: 4, 83: 3}
HW_VOWELS = frozenset("aeiou")  # 首尾元音省略判据（P-22.2.2.1.1）
HW_SIX_A = frozenset({O, S, Se, Te, Bi})  # Table 2.5 六元环 A 组
HW_SIX_B = frozenset({N, Si, Ge, Sn, Pb})  # Table 2.5 六元环 B 组（余者归 C 组）
HW_UNSAT_SIX = {"A": "ine", "B": "ine", "C": "inine"}  # 六元不饱和词干（P-22.2.2.1.6）
HW_SAT_SIX = {"A": "ane", "B": "inane", "C": "inane"}  # 六元饱和词干
HW_UNSAT_TAIL = {7: "epine", 8: "ocine", 9: "onine", 10: "ecine"}  # Table 2.5 七至十元不饱和
HW_SAT_TAIL = {7: "epane", 8: "ocane", 9: "onane", 10: "ecane"}
HW_ZH_RING = {3: "环丙", 4: "环丁", 5: "环戊", 6: "环己",  # 中文基干环前缀（表 3-4）
              7: "环庚", 8: "环辛", 9: "环壬", 10: "环癸"}

_MULT_EN_10 = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",  # 倍数前缀（P-14.2 Table 1.4；数量词与母链碳数同源，支持到 9999）en 数量词（deca/undeca/…/icosa/henicosa…）同时供 layer5.stems 生成长链母链词干，故下沉到本层：en_num_term 是唯一来源，stems 取词干仅去其尾 'a'。
               6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
_MULT_ZH_10 = {1: "", 2: "二", 3: "三", 4: "四", 5: "五",
               6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}
_EN_UNIT = {1: "hen", 2: "do", 3: "tri", 4: "tetra", 5: "penta",
            6: "hexa", 7: "hepta", 8: "octa", 9: "nona"}
_EN_DECADE = {1: "deca", 2: "icosa", 3: "triaconta", 4: "tetraconta",
              5: "pentaconta", 6: "hexaconta", 7: "heptaconta",
              8: "octaconta", 9: "nonaconta"}
_EN_HECTO = {1: "hecta", 2: "dicta", 3: "tricta", 4: "tetracta",  # Table 1.4 百位（101 = hen+hecta，100 = hecta）
             5: "pentacta", 6: "hexacta", 7: "heptacta", 8: "octacta", 9: "nonacta"}
_EN_KILO = {1: "kilia", 2: "dilia", 3: "trilia", 4: "tetralia",  # Table 1.4 千位（1001 = hen+kilia，1000 = kilia）
            5: "pentalia", 6: "hexalia", 7: "heptalia", 8: "octalia", 9: "nonalia"}
ZH_DIGITS = "一二三四五六七八九"


def en_num_term(n: int) -> str | None:
    """英文数值词干（末带 'a'）：1–9999 按个→十→百→千反序拼基本词（P-14.2.1.2）。"""
    if n < 1 or n > 9999:
        return None
    if n <= 10:
        return _MULT_EN_10[n]
    ones, tens = n % 10, (n // 10) % 10
    parts: list[str] = []
    if n % 100 == 11:  # P-14.2.1.1.1：11 独用 undeca（111 → undecahecta）
        parts.append("undeca")
    else:
        if ones:
            parts.append(_EN_UNIT[ones])
        if tens:
            dec = _EN_DECADE[tens]
            if tens == 2 and parts and parts[-1][-1] in "aeiou":  # P-14.2.1.2：icosa 在元音后省 i（22 → docosa）
                dec = "cosa"
            parts.append(dec)
    hundreds, thousands = (n // 100) % 10, n // 1000
    if hundreds:
        parts.append(_EN_HECTO[hundreds])
    if thousands:
        parts.append(_EN_KILO[thousands])
    return "".join(parts)


def zh_numeral(n: int) -> str | None:
    """中文数值词（倍数用，1 返空）：≤10 查表，11–99 按十位组合，100+ 拼百/千。"""
    if n < 1 or n > 9999:
        return None
    if n <= 10:
        return _MULT_ZH_10[n]
    if n < 100:
        tens, ones = divmod(n, 10)
        head = "十" if tens == 1 else f"{ZH_DIGITS[tens-1]}十"
        return head if ones == 0 else f"{head}{ZH_DIGITS[ones-1]}"
    thousands, rest = divmod(n, 1000)
    hundreds, rest = divmod(rest, 100)
    tens, ones = divmod(rest, 10)
    parts: list[str] = []
    if thousands:
        parts.append(f"{ZH_DIGITS[thousands-1]}千")
    if hundreds:
        parts.append(f"{ZH_DIGITS[hundreds-1]}百")
    elif thousands and (tens or ones):
        parts.append("零")  # 高位与低位间的空位补「零」（1001 → 一千零一）
    if tens:
        parts.append(f"{ZH_DIGITS[tens-1]}十")
    elif parts and parts[-1] != "零" and ones:
        parts.append("零")
    if ones:
        parts.append(ZH_DIGITS[ones-1])
    name = "".join(parts)
    if n in (100, 1000) and name.startswith("一"):  # 表 3-2：整百/整千用「百烷/千烷」；复合数保留前导「一」（151 → 一百五十一）
        return name[1:]
    return name


MULT_EN = {n: (en_num_term(n) or "") for n in range(1, 10000)}
MULT_ZH = {n: (zh_numeral(n) or "") for n in range(1, 10000)}

def zh_bridge_root(name: str) -> str:
    """桥后缀（氨基/氧基/硫基）前的中文烃基名去尾「基」。"""
    if name.endswith("羰基"):
        return name  # 羰基/甲氧羰基等酰基名保留「基」：金标作羰基氨基而非羰氨基
    return name[:-1] if name.endswith("基") else name

AMIDO_RETAINED = {  # P-66.1.1.4.3
    "acetyl": ("acetamido", "乙酰氨基"),
    "formyl": ("formamido", "甲酰胺基"),
    "benzoyl": ("benzamido", "苯甲酰胺基"),
}
AMIDO_RETAINED_EN = frozenset(v[0] for v in AMIDO_RETAINED.values())

# ── L0 简单分子 / 单质查表（规范 SMILES → (en, zh)）────
# 单元素与极简单分子无母体链/环，走不到 L2–L5，故在命名入口直接查表。名称取 IUPAC 保留名
# （P-21：water/ammonia/氢卤酸二元名）。键为 L0 电荷归一后的规范 SMILES，取值限制重原子 ≤ 3。
SIMPLE_MAX_HEAVY = 3  # 查表覆盖的重原子上限；超出者仍走常规管线
ELEMENT_METAL_NAMES = {  # 单质金属：元素符号 → (en, zh)；中性单原子，键在 SIMPLE_MOLECULES 补成 "[X]"
    "Li": ("lithium", "锂"), "Be": ("beryllium", "铍"), "Na": ("sodium", "钠"),
    "Mg": ("magnesium", "镁"), "Al": ("aluminium", "铝"), "K": ("potassium", "钾"),
    "Ca": ("calcium", "钙"), "Sc": ("scandium", "钪"), "Ti": ("titanium", "钛"),
    "V": ("vanadium", "钒"), "Cr": ("chromium", "铬"), "Mn": ("manganese", "锰"),
    "Fe": ("iron", "铁"), "Co": ("cobalt", "钴"), "Ni": ("nickel", "镍"),
    "Cu": ("copper", "铜"), "Zn": ("zinc", "锌"), "Ga": ("gallium", "镓"),
    "Ge": ("germanium", "锗"), "As": ("arsenic", "砷"), "Rb": ("rubidium", "铷"),
    "Sr": ("strontium", "锶"), "Y": ("yttrium", "钇"), "Zr": ("zirconium", "锆"),
    "Nb": ("niobium", "铌"), "Mo": ("molybdenum", "钼"), "Tc": ("technetium", "锝"),
    "Ru": ("ruthenium", "钌"), "Rh": ("rhodium", "铑"), "Pd": ("palladium", "钯"),
    "Ag": ("silver", "银"), "Cd": ("cadmium", "镉"), "In": ("indium", "铟"),
    "Sn": ("tin", "锡"), "Sb": ("antimony", "锑"), "Cs": ("caesium", "铯"),
    "Ba": ("barium", "钡"), "La": ("lanthanum", "镧"), "Ce": ("cerium", "铈"),
    "Hf": ("hafnium", "铪"), "Ta": ("tantalum", "钽"), "W": ("tungsten", "钨"),
    "Re": ("rhenium", "铼"), "Os": ("osmium", "锇"), "Ir": ("iridium", "铱"),
    "Pt": ("platinum", "铂"), "Au": ("gold", "金"), "Hg": ("mercury", "汞"),
    "Tl": ("thallium", "铊"), "Pb": ("lead", "铅"), "Bi": ("bismuth", "铋"),
    "Th": ("thorium", "钍"), "U": ("uranium", "铀"),
}
SIMPLE_MOLECULES: dict[str, tuple[str, str]] = {
    # 单元素分子（H₂/O₂/O₃）；N₂/F₂/Cl₂/Br₂/I₂/H₂O₂/N₂H₄ 由 L2 均一杂原子链命名，不在此覆盖
    "[H][H]": ("hydrogen", "氢"), "O=O": ("oxygen", "氧"),
    "O=[O+][O-]": ("ozone", "臭氧"),
    # 母体氢化物（P-21 第 15/16/17 族保留名与二元名）
    "O": ("water", "水"), "S": ("hydrogen sulfide", "硫化氢"), "P": ("phosphine", "膦"),
    "N": ("ammonia", "氨"), "[SeH2]": ("hydrogen selenide", "硒化氢"),
    "[TeH2]": ("hydrogen telluride", "碲化氢"), "[AsH3]": ("arsine", "胂"),
    "[SbH3]": ("stibine", "锑化氢"), "[SiH4]": ("silane", "硅烷"),
    "B": ("borane", "硼烷"),
    "F": ("hydrogen fluoride", "氟化氢"), "Cl": ("hydrogen chloride", "氯化氢"),
    "Br": ("hydrogen bromide", "溴化氢"), "I": ("hydrogen iodide", "碘化氢"),
    # 简单氧化物 / 硫化物 / 氮化物
    "NO": ("hydroxylamine", "羟胺"), "[N]=O": ("nitric oxide", "一氧化氮"),
    "[N-]=[N+]=O": ("nitrous oxide", "一氧化二氮"), "O=[N+][O-]": ("nitrogen dioxide", "二氧化氮"),
    "[C-]#[O+]": ("carbon monoxide", "一氧化碳"), "O=C=O": ("carbon dioxide", "二氧化碳"),
    "S=C=S": ("carbon disulfide", "二硫化碳"), "O=C=S": ("carbonyl sulfide", "氧硫化碳"),
    "O=S=O": ("sulfur dioxide", "二氧化硫"),
    # 简单卤素化合物
    "OCl": ("hypochlorous acid", "次氯酸"), "NCl": ("chloramine", "氯胺"),
    "ClOCl": ("dichlorine monoxide", "一氧化二氯"), "O=[Cl+][O-]": ("chlorine dioxide", "二氧化氯"),
    **{f"[{sym}]": nm for sym, nm in ELEMENT_METAL_NAMES.items()},
}

# ── L0 电荷归一 / 盐解离 ────────────────────
DONOR_KIND = ("carboxyl", "phospho", "sulfo")  # 允许作为强酸供体的酸类：羧酸/磷酸/磺酸。磷酸供体在「酰胺 O⁻ 受体」场景下才有产出（见 preprocessor 的二次互变归一，gold 把 N=C([O-]) 写成酰胺、把 P-OH 写成 oxidophosphoryl）；磺酸供体在「净负离子的胺受体」场景下才有产出（sulfonatooxy 两性离子式）。
ACCEPTOR_Z = frozenset({O, N})  # 弱受体允许的元素：O（酚氧/烯醇氧/酰胺氧）、N（去质子化氮）；保守可只留 {O}。
ACID_CENTERS = {          # 中心元素 → 酸类名；键序即供体搬运的酸强度序
    C: "carboxyl",        # C(=O)OH
    P: "phospho",         # P(=O)OH
    S: "sulfo",           # S(=O)nOH
}
METAL_ION_EN = {Li: "lithium", Na: "sodium", Mg: "magnesium", K: "potassium",
                Ca: "calcium"}  # 反离子金属阳离子（P-71.2）：原子序数 → 英文金属名
METAL_ION_ZH = {"lithium": "锂", "sodium": "钠", "potassium": "钾",
                "magnesium": "镁", "calcium": "钙"}
HALIDE_ZH = {F: "氟化物", Cl: "氯化物", Br: "溴化物", I: "碘化物"}  # 卤素阴离子 X⁻（P-71.2）
HALIDE_HX_EN = {F: "hydrofluoride", Cl: "hydrochloride",   # 中性卤化氢加合物 HX（P-71.3）
                Br: "hydrobromide", I: "hydroiodide"}
HALIDE_HX_ZH = {F: "氢氟酸盐", Cl: "盐酸盐", Br: "氢溴酸盐", I: "氢碘酸盐"}

# ── P-14.1 λ 约定 / P-21 母体氢化物 ────────────
STANDARD_BONDING_NUMBERS = {  # 原子序数 → 标准键数（P-14.1.2 表 1.3）；键数偏离此值才用 λ 记号，标准价一律不标
    B: 3, Al: 3, Ga: 3, In: 3, Tl: 3,
    C: 4, Si: 4, Ge: 4, Sn: 4, Pb: 4,
    N: 3, P: 3, As: 3, Sb: 3, Bi: 3,
    O: 2, S: 2, Se: 2, Te: 2, Po: 2,
    F: 1, Cl: 1, Br: 1, I: 1, At: 1,
}
PARENT_HYDRIDE_STEMS: dict[int, tuple[str, str]] = {  # 原子序数 → (en, zh) 母体氢化物词干（P-21 表 2.1）；碳不在此表，甲烷走既有碳链引擎
    B:  ("borane", "硼烷"),       Al: ("alumane", "铝烷"),   Ga: ("gallane", "镓烷"),
    In: ("indigane", "铟烷"),     Tl: ("thallane", "铊烷"),   # alane 禁用（与丙氨酸 alanine 冲突）、indane/aluminane 被烃环占用，均按 P-21.1.1 改用此名
    Si: ("silane", "硅烷"),       Ge: ("germane", "锗烷"),    Sn: ("stannane", "锡烷"),  Pb: ("plumbane", "铅烷"),
    N:  ("azane", "氮烷"),        P:  ("phosphane", "磷烷"),  As: ("arsane", "砷烷"),
    Sb: ("stibane", "锑烷"),      Bi: ("bismuthane", "铋烷"),
    O:  ("oxidane", "氧烷"),      S:  ("sulfane", "硫烷"),    Se: ("selane", "硒烷"),
    Te: ("tellane", "碲烷"),      Po: ("polane", "钋烷"),
    F:  ("fluorane", "氟烷"),     Cl: ("chlorane", "氯烷"),   Br: ("bromane", "溴烷"),
    I:  ("iodane", "碘烷"),       At: ("astane", "砹烷"),
}
# 母体氢化物 → 去氢取代基形式（P-21 表 2.1 / P-29）：(free_en, free_zh, yl_en, yl_zh, ylidene_en, ylidene_zh)
# 单键连母体取 -yl（sulfanyl），双键连母体取 -ylidene（sulfanylidene），λ 段原样保留
HYDRIDE_YL_FORMS: tuple[tuple[str, str, str, str, str, str], ...] = (
    ("borane", "硼烷", "boranyl", "硼基", "boranylidene", "硼亚基"),
    ("alumane", "铝烷", "alumanyl", "铝基", "alumanylidene", "铝亚基"),
    ("gallane", "镓烷", "gallanyl", "镓基", "gallanylidene", "镓亚基"),
    ("indigane", "铟烷", "indiganyl", "铟基", "indiganylidene", "铟亚基"),
    ("thallane", "铊烷", "thallanyl", "铊基", "thallanylidene", "铊亚基"),
    ("silane", "硅烷", "silyl", "甲硅烷基", "silylidene", "亚甲硅基"),
    ("germane", "锗烷", "germyl", "甲锗烷基", "germylidene", "亚甲锗烷基"),
    ("stannane", "锡烷", "stannyl", "甲锡烷基", "stannylidene", "亚甲锡烷基"),
    ("plumbane", "铅烷", "plumbyl", "甲铅烷基", "plumbylidene", "亚甲铅烷基"),
    ("phosphane", "磷烷", "phosphanyl", "磷烷基", "phosphanylidene", "磷亚甲基"),
    ("arsane", "砷烷", "arsanyl", "砷烷基", "arsanylidene", "砷亚甲基"),
    ("stibane", "锑烷", "stibanyl", "锑烷基", "stibanylidene", "锑亚甲基"),
    ("bismuthane", "铋烷", "bismuthanyl", "铋烷基", "bismuthanylidene", "铋亚甲基"),
    ("oxidane", "氧烷", "oxidanyl", "氧基", "oxidanylidene", "氧亚基"),
    ("sulfane", "硫烷", "sulfanyl", "硫代", "sulfanylidene", "硫亚基"),
    ("selane", "硒烷", "selanyl", "硒代", "selanylidene", "硒亚基"),
    ("tellane", "碲烷", "tellanyl", "碲代", "tellanylidene", "碲亚基"),
    ("polane", "钋烷", "polanyl", "钋代", "polanylidene", "钋亚基"),
    ("fluorane", "氟烷", "fluoranyl", "氟代", "fluoranylidene", "氟亚基"),
    ("chlorane", "氯烷", "chloranyl", "氯代", "chloranylidene", "氯亚基"),
    ("bromane", "溴烷", "bromanyl", "溴代", "bromanylidene", "溴亚基"),
    ("iodane", "碘烷", "iodanyl", "碘氧基", "iodanylidene", "碘亚基"),
    ("astane", "砹烷", "astanyl", "砹代", "astanylidene", "砹亚基"),
)

# 多核母体氢化物中文倍数表见文件末 HYDRIDE_MULT_ZH（依赖 HS_NUMBER，须在其后定义）

# ── L2/L5 单核母体氢化物（P-15.4.1）────────────
MONONUCLEAR_HYDRIDES: dict[str, tuple] = {
    "oxidane":    ("O", "氧化烷", ("hydroxy",   "羟基"),     ("oxy",      "氧基"), "氧基"),
    "azane":      ("N", "氮烷",   ("amino",     "氨基"),     ("amino",    "氨基"), "氨基"),
    "sulfane":    ("S", "硫烷",   ("sulfanyl",  "硫基"),     ("sulfanyl", "硫基"), "硫基"),
    # 自由价落在阳离子上（P-73.1.1）：词干取阳离子名，去氢名即 铵基 一族前缀
    "azanium":    ("N", "铵",     ("azaniumyl",    "铵基"),     None, "铵基"),
    "oxidanium":  ("O", "氧鎓",   ("oxidaniumyl",  "氧鎓基"),   None, "氧鎓基"),
    "phosphanium":("P", "鏻",     ("phosphaniumyl", "鏻基"),    None, "鏻基"),
    "sulfanium":  ("S", "硫鎓",   ("sulfonio",     "硫鎓基"),   None, "硫鎓基"),
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
PHOSPHORYL_STEMS = tuple(en for en, v in MONONUCLEAR_HYDRIDES.items() if v[0] == "P" and en != "phosphanium")  # 替代碳词干的 P 酰基词干（避免 P 被当碳中心）；阳离子词干不按酰基拼接

# ── L2 单核母体阳离子（P-73.1.1）────────────
# 表 2.1 去词尾 'ne' 的系统名（P-73.1.1.2：-ane 换成 -ium），本次只登记元素支持集内的行
MONONUCLEAR_CATION_SYSTEMATIC = {
    B:  "boranium",      # borane
    C:  "methanium",     # methane（保留名，泛 'ane' 命名法的碳基准）
    N:  "azanium",       # azane
    O:  "oxidanium",     # oxidane
    F:  "fluoranium",    # fluorane
    P:  "phosphanium",   # phosphane
    S:  "sulfanium",     # sulfane
    Cl: "chloranium",    # chlorane
    Br: "bromanium",     # bromane
    I:  "iodanium",      # iodane
}
# 表 7.3 第 15/16/17 族单核母体阳离子的保留名（仅一般命名；中文取此表，P 按测试集作「鏻」）
CATION_RETAINED_NAMES = {
    N:  ("ammonium",    "铵"),
    P:  ("phosphonium", "鏻"),
    O:  ("oxonium",     "氧鎓"),
    S:  ("sulfonium",   "硫鎓"),
    F:  ("fluoronium",  "氟鎓"),
    Cl: ("chloronium",  "氯鎓"),
    Br: ("bromonium",   "溴鎓"),
    I:  ("iodonium",    "碘鎓"),
}


CATION_FREE_STEMS = {z: en for z, en in MONONUCLEAR_CATION_SYSTEMATIC.items()
                     if en in MONONUCLEAR_HYDRIDES}  # 带电锚点的自由价词干（有氢化物行的元素）
CATION_STEMS = frozenset(CATION_FREE_STEMS.values())  # 阳离子词干名（P-73.1.1），L5 据此保留烃基尾「基」
CATION_YL_STEMS = frozenset(v[1] for en, v in MONONUCLEAR_YL.items()
                            if en in CATION_STEMS)  # 阳离子母体派生的去氢前缀（azaniumyl/oxidaniumyl/…）；为复合前缀，倍数用 bis(azaniumyl)（P-16.3.2）


def cation_parent_names(z: int) -> tuple[str, str] | None:
    """单核母体阳离子的 (en, zh)：英文取表 2.1 系统名，中文取表 7.3 保留名；缺一即无。"""
    en = MONONUCLEAR_CATION_SYSTEMATIC.get(z)
    retained = CATION_RETAINED_NAMES.get(z)
    if en is None or retained is None:
        return None
    return en, retained[1]

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
# 中心原子自任母体的含氧酸 kind（L1/L2/L5 共用）；碳锚定的磺酸族不在内
OXO_CENTER_KINDS = frozenset({"phosphate", "phosphonate", "sulfate", "boronic"})

# ── L4 位次 / 编号 ──────────────────────
HYDRO_MULT_N = frozenset({2, 4, 6, 8, 10, 12, 14, 16, 18, 20})  # 加氢前缀覆盖的氢原子数（P-31.2.2 以偶数倍增前缀表示双键的饱和，位次数为加氢原子数）；数量词本身取自本层（唯一来源），此处只表达 L4 的适用域，域外放弃而非给错名。
TRADITIONAL_NUMBERING_IDS = frozenset({  # P-25.3.3：这些保留骨架按传统编号，不走 P-25.3.3.1 的优选取向自动编号。xanthene 及其硫属类似物（xanthene/thioxanthene）与 cyclopenta[a]phenanthrene（甾体 1-17）已按_TEMPLATES 的 standard 字段登记传统编号，故一并列入。
    "anthracene", "phenanthrene", "acridine", "carbazole", "purine",
    "xanthene", "thioxanthene", "cyclopenta[a]phenanthrene",
})
RS_HI = frozenset({"R", "M", "r"})   # 编号优先级高的 CIP 描述符（P-91.2）
RS_LO = frozenset({"S", "P", "s"})   # 与之成对、取较低位次的次位描述符

# ── L5 组装词表 ──────────────────────
BIS_EN = {2: "bis", 3: "tris",  # P-16.3.2 复合前缀倍增（bis/tris，非 di/tri）；4+ 按 P-14.2.2 于基本词后加 kis（tetrakis…）
          **{n: f"{en_num_term(n)}kis" for n in range(4, 10000)}}
BIS_ZH = {2: "双", **{n: zh_numeral(n) for n in range(3, 10000)}}  # 与 BIS_EN 同域；2 用「双」（单双叁肆），3+ 用中文数字
BIS_EN_SET = frozenset(BIS_EN.values())  # 供 `mult in …` 常数判据：万级 dict 的 .values() 线性扫描会拖慢管线
BRIDGE_SUFFIX_EN = ("oxy", "sulfanyl", "amino")   # O/S/N 桥后缀（P-63.2.2.1）：平铺 [-yl]oxy 与合一 [-yloxy] 均属合法
BRIDGE_SUFFIX_ZH = ("氧基", "硫基", "氨基")        # 与 BRIDGE_SUFFIX_EN 同序同位
DIATOMIC_BRIDGE_YL = ("diazenyl", "disulfanyl")   # 双原子桥合一保留前缀：-N=N-R / -S-S-R 收成 R-diazenyl / R-disulfanyl
SIMPLE_BRIDGE_YL_NO_PAREN = frozenset(f"phenyl{s}" for s in DIATOMIC_BRIDGE_YL)  # 裸苯基前端 + 双原子桥：作前缀免括
BRIDGE_FUSION_YL: dict[tuple[str, str], tuple[tuple[str, ...], str, str]] = {
    # (中心单核氢化物词干, 前端名尾 en) → (前端名尾 zh 候选, 合一前缀 en, 合一前缀 zh)
    ("azane", "imino"): (("亚氨基",), "diazenyl", "二氮烯基"),            # P-68.3.1.3：diazenyl 名优先于 azo 名
    ("sulfane", "disulfanyl"): (("二硫代基", "二硫代"), "trisulfanyl", "三硫代"),  # P-68.3.1.4：-S-S-S- 三硫链（R-二硫代 → R-三硫代）
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
    "benzo[d]1,2-thiazole": ("1,2-benzothiazole", "1,2-苯并噻唑"),  # 硫属类似物，位次形态与 oxazole 行同构
    "pyrido[3,2-b]pyridine": ("1,5-naphthyridine", "1,5-萘啶"),   # 双氮稠合：两 N 分处两环且位次一致
    "pyrido[2,1-a]isoquinoline": ("benzo[a]quinolizine", "苯并[a]喹嗪"),
    "benzo[b]anthracene":  ("tetracene",     "并四苯"),
    "benzo[c]thiophene":   ("2-benzothiophene", "2-苯并噻吩"),  # 异苯并噻吩：S 占 2 位（≠ benzo[b]thiophene 的 1- 异构）
    "benzo[d]azepine":     ("3-benzazepine", "3-苯并氮杂卓"),    # N 占 3 位；稠合位次与亚甲基分布均与稠合名一一对应
}
HS_NUMBER = "甲乙丙丁戊己庚辛壬癸"
# 多核母体氢化物中文倍数用「甲乙丙丁…」而非「二三…」（第3章 3.2.2：diazane=乙氮烷、trisulfane=丙硫烷、pentasilane=戊硅烷）
HYDRIDE_MULT_ZH = {n: HS_NUMBER[n - 1] for n in range(1, 11)}


def hydride_chain_stem(z: int, n: int) -> tuple[str, str] | None:
    """n 核均一母体氢化物裸词干（倍数词 + 去尾词干），供链引擎拼 ene/yne（P-21.2.2）。"""
    forms = PARENT_HYDRIDE_STEMS.get(z)
    if forms is None or n < 2:
        return None
    mult_en, mult_zh = en_num_term(n), (HYDRIDE_MULT_ZH.get(n) or zh_numeral(n))
    return (f"{mult_en}{forms[0][:-3]}", f"{mult_zh}{forms[1][:-1]}")
