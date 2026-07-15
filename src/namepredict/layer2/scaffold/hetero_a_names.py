"""a-prefix and cycloalkane size tables for skeletal replacement stems.

Local 3–12 size maps (L2 must not import L5). Unknown Z/size → None.
"""
from __future__ import annotations

# Z → EN a-prefix (IUPAC skeletal replacement)
A_PREFIX_EN: dict[int, str] = {
    8: "oxa",
    16: "thia",
    7: "aza",
    15: "phospha",
    34: "selena",
    5: "bora",
    14: "sila",
}

# Z → ZH element root (杂 applied by stem assembler)
A_ROOT_ZH: dict[int, str] = {
    8: "氧",
    16: "硫",
    7: "氮",
    15: "磷",
    34: "硒",
    5: "硼",
    14: "硅",
}

# ring size → full cycloalkane name (total ring atoms, incl. hetero)
CYCLO_EN: dict[int, str] = {
    3: "cyclopropane",
    4: "cyclobutane",
    5: "cyclopentane",
    6: "cyclohexane",
    7: "cycloheptane",
    8: "cyclooctane",
    9: "cyclononane",
    10: "cyclodecane",
    11: "cycloundecane",
    12: "cyclododecane",
}

CYCLO_ZH: dict[int, str] = {
    3: "环丙烷",
    4: "环丁烷",
    5: "环戊烷",
    6: "环己烷",
    7: "环庚烷",
    8: "环辛烷",
    9: "环壬烷",
    10: "环癸烷",
    11: "环十一烷",
    12: "环十二烷",
}

MULT_EN: dict[int, str] = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta"}
MULT_ZH: dict[int, str] = {1: "", 2: "二", 3: "三", 4: "四", 5: "五"}


def a_prefix_en(z: int) -> str | None:
    return A_PREFIX_EN.get(z)


def a_root_zh(z: int) -> str | None:
    return A_ROOT_ZH.get(z)


def cyclo_ane_en(size: int) -> str | None:
    return CYCLO_EN.get(size)


def cyclo_ane_zh(size: int) -> str | None:
    return CYCLO_ZH.get(size)
