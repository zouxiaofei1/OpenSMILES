"""locant 串排序键: 数字位次与字母位次("4a")统一排序。"""
from __future__ import annotations

import re


def locant_key(x) -> tuple[int, str]:
    """locant → 排序键: 数字按数值、字母尾作次级键，保证 "4" < "4a" < "5" < "10"。"""
    m = re.match(r"(\d+)([a-z]*)", str(x))
    return (int(m.group(1)), m.group(2)) if m else (0, "")


def locant_str_sort(locs) -> list:
    """按 locant_key 排序 locant 集合(兼容 int 与 "4a" 混合)。"""
    return sorted(locs, key=locant_key)
