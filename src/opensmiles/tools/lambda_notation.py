"""λ 约定（P-14.1）：非标准键数计算与记号渲染；标准价不标 λ。"""
from __future__ import annotations

from rdkit.Chem import Atom, Mol

from opensmiles.constants import STANDARD_BONDING_NUMBERS


def bonding_number(atom: Atom) -> int:
    """骨架原子的键数 n（P-14.1.1）；GetTotalValence 对芳香键已按交替单双键计。"""
    return int(atom.GetTotalValence())


def is_nonstandard(atom: Atom) -> bool:
    """该原子的键数是否偏离标准值（P-14.1.2 表 1.3）。"""
    std = STANDARD_BONDING_NUMBERS.get(atom.GetAtomicNum())
    return std is not None and bonding_number(atom) != std


def is_lambda_marked(atom: Atom) -> bool:
    """该原子是否应标 λn：键数偏离标准值且电中性（P-14.1.2；带电原子归 -ium/-ide 路径）。"""
    return atom.GetFormalCharge() == 0 and is_nonstandard(atom)


def nonstandard_bonding(mol: Mol) -> dict[int, int]:
    """全分子的非标准键数原子：{原子 idx: 键数}；带电原子归 -ium/-ide 路径，不标 λ。"""
    return {a.GetIdx(): bonding_number(a) for a in mol.GetAtoms()
            if is_lambda_marked(a)}


def lambda_mark(n: int, *, en: bool = False) -> str:
    """λ 记号本体：中文 λ5、英文 lambda5（测试集口径，P-14.1.3）。"""
    return f"{'lambda' if en else 'λ'}{n}"


def put_lambda(locant: str, n: int, *, en: bool = False) -> str:
    """位次与 λ 连写，中间无连字符（P-15.4.1.3）：1 + 6 → 1λ6。"""
    return f"{locant}{lambda_mark(n, en=en)}"


def locant_lambda_str(items, *, en: bool = False) -> str:
    """位次串：[(位次, 键数|None)] 升序 → 1λ6,2；键数为 None 的位次不带 λ。"""
    parts = (put_lambda(str(loc), n, en=en) if n else str(loc) for loc, n in items)
    return ",".join(parts)


def delta_mark(c: int) -> str:
    """δ 记号（P-25.7.2）：该原子的连续形式双键数 c → δ2，紧随 λn 之后。"""
    return f"δ{c}"
