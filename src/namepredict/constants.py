"""跨层共享常量：原子序数、倍数前缀与名称文本规范化。"""

from __future__ import annotations
import re

# ── 原子序数 ───────────────────────────────────────────────────
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

# ── 常用集合 ───────────────────────────────────────────────────
HALO_Z = frozenset({F, Cl, Br, I})
RING_HETERO = frozenset({N, O, S})  # 环内杂原子：其单碳酰基按环酮命名（内酰胺/内酯/硫代内酯、N-酰基环胺），见 layer1.analyzer._is_ketone_carbon。
HALO_EN = {F: "fluoro", Cl: "chloro", Br: "bromo", I: "iodo"}
HALO_ZH = {F: "氟", Cl: "氯", Br: "溴", I: "碘"}
HALIDE_EN = {F: "fluoride", Cl: "chloride", Br: "bromide", I: "iodide"}
N_PREFIX_KINDS = frozenset({"n_alkyl", "n_phenyl", "n_benzyl", "n_block"})  # N-取代基 kind（P-62.2.2.1）：走 N- 前缀、位次以 N 标注或隐含省略，不参与数字位次通道。

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
_ZH_DIGITS = "一二三四五六七八九"


def en_num_term(n: int) -> str | None:
    """英文数值词干（数量词/链词共用，末带 'a'）：≤10 查表，11=undeca，12–99 按个位(hen/do)+十位组合。"""
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
    head = "十" if tens == 1 else f"{_ZH_DIGITS[tens-1]}十"
    return head if ones == 0 else f"{head}{_ZH_DIGITS[ones-1]}"


MULT_EN = {n: (en_num_term(n) or "") for n in range(1, 100)}
MULT_ZH = {n: (zh_numeral(n) or "") for n in range(1, 100)}

def zh_bridge_root(name: str) -> str:
    """桥后缀（氨基/氧基/硫基）前的中文烃基名去尾「基」：甲基→甲、叔丁基→叔丁、环己基→环己、丙-2-基→丙-2-（tiers gold 口径）。"""
    return name[:-1] if name.endswith("基") else name

AMIDO_RETAINED = {  # P-66.1.1.4.3
    "acetyl": ("acetamido", "乙酰氨基"),
    "formyl": ("formamido", "甲酰胺基"),
    "benzoyl": ("benzamido", "苯甲酰胺基"),
}
AMIDO_RETAINED_EN = frozenset(v[0] for v in AMIDO_RETAINED.values())  

# ── 文本规范化 ──────────────────────────────────────────────────
_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    """规范化英文名：小写、去重空白并统一连字符/逗号/括号。"""
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = s.replace("[", "(").replace("]", ")")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    """规范化中文名：去首尾空白并统一括号种类（方/圆等价，仅括注外观不同不判分）。"""
    s = (name or "").strip()
    s = s.replace("[", "(").replace("]", ")")
    return s


def nospace(name: str) -> str:
    """删去全部空白字符 """
    return "".join((name or "").split())
