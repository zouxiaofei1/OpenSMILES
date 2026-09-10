"""L4 加氢程度前缀（P-31.2.2）：按编号后的 locant 表达 'hydro' 修饰，完全氢化省略位次（P-14.3.4.5）。"""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer4.locant_key import locant_key

# 本前缀覆盖的加氢原子数（P-31.2.2 以偶数倍增前缀表示双键的饱和，位次数为加氢原子数）；
# 数量词本身取自 constants（唯一来源），此处只表达 L4 的适用域，域外放弃而非给错名。
HYDRO_MULT_N = frozenset({2, 4, 6, 8, 10, 12, 14, 16, 18, 20})


def hydro_prefix(chain, labels, hydro_atoms) -> tuple[str, str]:
    """返回加氢前缀 (en, zh)：位次取加氢原子 locant（按位次升序），倍增数为加氢原子数；环原子全加氢时省略位次（P-14.3.4.5，如 decahydronaphthalene）。无加氢或数量超表返回空串对。"""
    chain = list(chain or ())
    hydro = set(hydro_atoms or ())
    if not chain or not hydro:
        return "", ""
    n = len(hydro)
    if n not in HYDRO_MULT_N:
        return "", ""
    if set(chain) <= hydro:  # 完全氢化：省略全部位次（P-14.3.4.5）
        return f"{MULT_EN[n]}hydro", f"{MULT_ZH[n]}氢"
    use_labels = bool(labels) and len(labels) == len(chain)
    locs = [labels[chain.index(a)] if use_labels else str(chain.index(a) + 1)
            for a in hydro if a in chain]
    if len(locs) != n:  # 加氢原子不全在编号链内：位次无法完整表达，放弃而非给错名
        return "", ""
    locs.sort(key=locant_key)
    loc = ",".join(locs)
    return f"{loc}-{MULT_EN[n]}hydro", f"{loc}-{MULT_ZH[n]}氢"
