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
HALO_EN = {F: "fluoro", Cl: "chloro", Br: "bromo", I: "iodo"}
HALO_ZH = {F: "氟", Cl: "氯", Br: "溴", I: "碘"}

# ── 倍数前缀 ────────────────────────────────────────────────────
MULT_EN = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",
           6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
MULT_ZH = {1: "", 2: "二", 3: "三", 4: "四", 5: "五",
           6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}

# ── 文本规范化 ──────────────────────────────────────────────────
_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    """规范化英文名：小写、去重空白并统一连字符/逗号。"""
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    """规范化中文名：去除首尾空白。"""
    return (name or "").strip()
