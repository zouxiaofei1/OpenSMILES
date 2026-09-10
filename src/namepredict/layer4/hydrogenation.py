"""L4 加氢程度前缀（P-31.2.2）：按编号后的 locant 表达 'hydro' 修饰，完全氢化省略位次（P-14.3.4.5）。"""
from __future__ import annotations

# 倍增前缀表：键为加氢原子数（P-31.2.2 以偶数倍增前缀表示双键的饱和，位次数为加氢原子数）。
_MULT_EN = {2: "di", 4: "tetra", 6: "hexa", 8: "octa", 10: "deca",
            12: "dodeca", 14: "tetradeca", 16: "hexadeca", 18: "octadeca", 20: "icosa"}
_MULT_ZH = {2: "二", 4: "四", 6: "六", 8: "八", 10: "十",
            12: "十二", 14: "十四", 16: "十六", 18: "十八", 20: "二十"}


def _locant_sort_key(loc: str) -> tuple[int, str]:
    """位次排序键：数字升序，同数字的桥头字母位（4a）紧随其后。"""
    digits = "".join(c for c in loc if c.isdigit())
    letters = "".join(c for c in loc if not c.isdigit())
    return (int(digits) if digits else 0, letters)


def hydro_prefix(chain, labels, hydro_atoms) -> tuple[str, str]:
    """返回加氢前缀 (en, zh)：位次取加氢原子 locant（按位次升序），倍增数为加氢原子数；环原子全加氢时省略位次（P-14.3.4.5，如 decahydronaphthalene）。无加氢或数量超表返回空串对。"""
    chain = list(chain or ())
    hydro = set(hydro_atoms or ())
    if not chain or not hydro:
        return "", ""
    n = len(hydro)
    if n not in _MULT_EN:
        return "", ""
    if set(chain) <= hydro:  # 完全氢化：省略全部位次（P-14.3.4.5）
        return f"{_MULT_EN[n]}hydro", f"{_MULT_ZH[n]}氢"
    use_labels = bool(labels) and len(labels) == len(chain)
    locs = [labels[chain.index(a)] if use_labels else str(chain.index(a) + 1)
            for a in hydro if a in chain]
    if len(locs) != n:  # 加氢原子不全在编号链内：位次无法完整表达，放弃而非给错名
        return "", ""
    locs.sort(key=_locant_sort_key)
    loc = ",".join(locs)
    return f"{loc}-{_MULT_EN[n]}hydro", f"{loc}-{_MULT_ZH[n]}氢"
