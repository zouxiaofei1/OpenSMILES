"""L1 官能团的局部环境 SMARTS 表（唯一事实来源）。
每条 = (FG 键, SMARTS)，首原子即中心原子，同名多条取并集。
FG 键与 fg_registry.FG_SPECS 一致；新增官能团只加表项不改检测代码。
"""
from __future__ import annotations

from rdkit import Chem

from namepredict.constants import (
    Al, As, B, Bi, Br, C, Cl, F, Ga, Ge, I, In, N, O, P, PARENT_SENIOR_ATOMS, Pb, Po, S, Sb, Se,
    Si, Sn, STANDARD_BONDING_NUMBERS, Te, Tl,
)

_ACID_O = "[$([#8;H1]),$([#8;-1;X1;H0])]"  # 酸性氧：与单一碳相连的羟基氧，或羧酸盐阴离子氧
_NOT_ACID = "!$([#6;X3](=[#8;X1])%s)" % _ACID_O  # 非羧酸碳
_NOT_ACYCLIC_ESTER = "!$([#6;X3]~[#8;X2;H0;+0;!R]~[#6])"  # 非环酯烷氧基氧（内酯的环内氧不算，留给酮分支）
_NOT_HALO = "!$([#6;X3]~[#9,#17,#35,#53])"  # 无卤素邻居
_ONE_C = "!$([#6;X3](~[#6])~[#6])"  # 至多一个碳邻居
_RING_HET = "~[#7,#8,#16;R]"  # 环内杂原子邻居（内酰胺/内酯/硫代内酯、N-酰基环胺）
_O_PHOS = "(-[$([#8;X2;H1;+0]),$([#8;X1;-1;H0]),$([#8;X2;H0;+0]~[#6,#15])])"
_O_ARM = "(-[#8;X2;H0;+0]~[#6])"      # O-R 臂（硫酸酯/磺酸酯的 O 侧）
_O_S_ARM = "(-[#8;X2;H0;+0]~[#16])"   # S-O-S 桥臂（多硫酸链；P-O-P 桥见 _O_PHOS 的磷臂）
_HALO_ARM = "(-[#9,#17,#35,#53])"     # 卤素臂（磺酰卤）
_N_ARM = "(-[#7;!R])"                 # 非环氮臂（磺酰胺；环内 S-N 走磺酰前缀）
_C_ARM = "(-[#6])"                    # 直连碳臂（碳骨架成员）
_SNY_ARM = "(-[#6,#7,#16])"           # 膦酸的第三臂（C/S/N 均可）

CATION_Z = (N, P, O, S, F, Cl, Br, I)


def _cation_local(z: int) -> str:
    """单核阳离子中心的局部模式：非环、半径 1 内无负形式电荷原子。"""
    return "[#%d;+1;!R;!$([#%d;+1]~[-1])]" % (z, z)


_HETERANE_V = (3, 4, 5, 6, 7)  # 非标准键数上界；v 为 RDKit 总价（含氢），即 P-14.1.1 的键数
_HETERANE_Z = tuple(z for z in STANDARD_BONDING_NUMBERS if z != C)


HETERANE_CHAIN_Z = tuple(z for z in STANDARD_BONDING_NUMBERS if z not in (C, B))  # 硼烷按 P-21.2.2 排除

HYDRIDE_STD_Z = frozenset({B, Al, Ga, In, Tl, Si, Ge, Sn, Pb, P, As, Sb, Bi})


def _oxo_guard(z: int) -> str:
    """排除带 =O（含氧酸）与非硫的 =S（磷硫酰）中心。"""
    bans = ("[#8]",) if z == 16 else ("[#8]", "[#16]")
    return "".join(f";!$([#{z}]={b})" for b in bans)


def _heterane_local(z: int, v: int) -> str:
    """杂原子烃中心的局部模式：中性、键数为 v、不带 oxo/thioxo。"""
    return f"[#{z};v{v};+0{_oxo_guard(z)}]"


def _standard_heterane_local(z: int) -> str:
    """标准价单核母体中心：非环、中性、无 oxo、邻居只含 C/H。"""
    return (f"[#{z};v{STANDARD_BONDING_NUMBERS[z]};+0;!R{_oxo_guard(z)}"
            f";!$([#{z}]~[!#6;!#1])]")


HETERANE_MIN_CHAIN = {N: 3, O: 3, S: 3, Se: 3, Te: 3, Po: 3}  # 取代链须 ≥3 连的元素（P-68.4.0）
HETERANE_CHAIN_MIN_DEFAULT = 2  # 其余元素（Si/Ge/Sn/Pb、P/As/Sb/Bi 等）：≥2 连即可带取代


def _is_parent_hydride(mol, comp: set[int], z: int) -> bool:
    """该同元素集合是否为「由氢饱和」的母体氢化物。"""
    for i in comp:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetAtomicNum() not in (1, z):
                return False
    return True


def is_standard_parent_hydride_center(mol, idx: int) -> bool:
    """标准价单杂原子是否为母体中心（P-44.1.2）。"""
    z = mol.GetAtomWithIdx(idx).GetAtomicNum()
    if z not in HYDRIDE_STD_Z:
        return False
    if sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == z) != 1:
        return False
    present = {a.GetAtomicNum() for a in mol.GetAtoms()}
    return next((s for s in PARENT_SENIOR_ATOMS if s in present), None) == z


def heterane_chain_ok(mol, idx: int) -> bool:
    """该原子能否作杂原子链中心（按元素最小链长）。"""
    from namepredict.tools.lambda_notation import is_nonstandard

    if is_nonstandard(mol.GetAtomWithIdx(idx)):
        return True  # 非标准键数即单核杂原子烃（P-14.1.3），不受链长度约束
    z = mol.GetAtomWithIdx(idx).GetAtomicNum()
    comp, stack = {idx}, [idx]
    while stack:
        for n in mol.GetAtomWithIdx(stack.pop()).GetNeighbors():
            if n.GetAtomicNum() == z and n.GetFormalCharge() == 0 and n.GetIdx() not in comp:
                comp.add(n.GetIdx())
                stack.append(n.GetIdx())
    if len(comp) == 1:
        return is_standard_parent_hydride_center(mol, idx)  # 标准价单核母体氢化物（P-44.1.2）
    if _is_parent_hydride(mol, comp, z):
        return True
    return len(comp) >= HETERANE_MIN_CHAIN.get(z, HETERANE_CHAIN_MIN_DEFAULT)


def _heterane_chain_local(z: int) -> tuple[str, ...]:
    """杂原子链中心的局部模式（P-21.2.2，按元素分档）。"""
    base = f"[#{z};+0;!R{_oxo_guard(z)}"  # !R：仅限无环链，环内杂原子走环系路径
    nz = f"[#{z};+0]"  # 链上邻居也须中性：叠氮 N=[N+]=[N-] 一类带电链不是母体氢化物
    pats = [f"{base};!$([#{z}]~[!#1;!#{z}]);$([#{z}]~{nz})]"]      # 纯母体氢化物：≥2 连
    pats.append(f"{base};$([#{z}]~{nz}~{nz})]" if z in HETERANE_MIN_CHAIN
                else f"{base};$([#{z}]~{nz})]")                     # 带取代/支链：≥3 连或 ≥2 连
    return tuple(pats)


FG_SMARTS: tuple[tuple[str, str], ...] = (
    ("radical", "[!#1]~[#0]"),  # 自由基：哑原子的重原子邻居（中心为该重原子，非哑原子）
    ("acyl", "[#6;X3;!$([#6]-[!#6;!#0;!#1])](=[#8;X1])(-[#6])~[#0]"),  # 酰基残基：哑原子所连的酰基头羰基碳；双键邻居不参与杂原子判定
    ("acid", f"[#6;X3](=[#8;X1]){_ACID_O}"),
    ("ester", f"[#6;X3;{_NOT_ACID}](=[#8;X1])[#8;X2;H0;+0;!R]~[#6]"),
    ("ester", f"[#6;X3;{_NOT_ACID}](=[#8;X1])[#16;X2;H0;+0;!R]~[#6]"),
    ("acyl_halide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#8;X1])~[#9,#17,#35,#53]"),
    ("amide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#8;X1])-[#7;X2,X3;!R;!$([#7;X2,X3]~[#7]~[#6;X3](=[#8,#16])~[#7])]"),
    ("amide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#16;X1])-[#7;X2,X3;!R;!$([#7;X2,X3]~[!#6;!#8]);!$([#7;X3]~[#8;X2]~[#6])]"),
    ("amide", f"[#6;X3;!R;!$([#6](~[#7])(~[#7])~[#7]);{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#7;X2,X3])-[#7;X2,X3;!R;!$([#7;X2,X3]~[!#6;!#8]);!$([#7;X3]~[#8;X2]~[#6])]"),
    ("aldehyde", f"[#6;X3;H1,H2;{_ONE_C};{_NOT_ACID};{_NOT_ACYCLIC_ESTER};{_NOT_HALO}](=[#8;X1])"),
    ("nitrile", "[#6;+0;X2]#[#7;X1;+0]"),
    ("ketone", f"[#6;X3;{_NOT_ACID}](=[#8;X1])(~[#6])~[#6]"),
    ("ketone", f"[#6;X3;{_NOT_ACID};{_ONE_C}](=[#8;X1])(~[#6]){_RING_HET}"),
    ("ketone", f"[#6;X3;{_NOT_ACID};R;{_ONE_C};{_NOT_ACYCLIC_ESTER}](=[#8;X1])~[#7,#8,#16]"),
    ("thione", "[#6;X3]=[#16;X1]"),
    ("oxoacid", "[#15;H0;+0](=[#8;X1;H0;+0])" + _O_PHOS * 3),                    # 磷酸 P(=O)(O)3
    ("oxoacid", "[#15;H0;+0](=[#8;X1;H0;+0])" + _O_PHOS * 2 + _SNY_ARM),         # 膦酸 P(=O)(O)2-X
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS + _C_ARM),   # 磺酸/磺酸盐
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _HALO_ARM + _C_ARM),  # 磺酰卤
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _N_ARM + _C_ARM),     # 磺酰胺
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_ARM + _C_ARM),     # 磺酸酯
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS + _O_ARM),    # 硫酸氢酯/硫酸酯
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS * 2),        # 硫酸/硫酸根/二硫酸链：两臂皆酸式氧或 O-臂
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_PHOS + _O_S_ARM),  # 二硫酸链端节：酸式/酯臂 + S-O-S 桥
    ("oxoacid", "[#16;H0;+0](=[#8;X1;H0;+0])(=[#8;X1;H0;+0])" + _O_S_ARM * 2),       # 多硫酸链中节：两桥
    ("alcohol", f"[#8;X2;H1]~[#6;!$([#6;X3](=[#8;X1]){_ACID_O})]"),
    ("alcohol", f"[#8;X1;-1;H0]~[#6;!$([#6;X3](=[#8;X1]){_ACID_O});!$([#6]=[#7])]"),
    ("thiol", "[#16;X2;H1]~[#6]"),
    ("thiol", "[#16;X1;-1;H0]~[#6]"),  # 去质子硫负离子（硫醇盐）：后缀转 -thiolate
    ("amine", "[#7;!R;!a;+0;H2,H3,H4;!$([#7](~[#6])~[#6])]~[#6]"),
    ("amine", "[#7;!R;!a;+0;H1;!$([#7](~[#6])(~[#6])~[#6])](~[#6])~[#6]"),
    ("amine", "[#7;!R;!a;+0;H0;!$([#7](~[#6])(~[#6])(~[#6])~[#6])](~[#6])(~[#6])~[#6]"),
    *[("cation", _cation_local(z)) for z in CATION_Z],  # 单核母体阳离子（P-73.1.1）：每元素一条，中心即阳离子原子
    ("azanide", "[#7;-1;!a;!$([#7;-1]~[+1])]"),
    *[("heterane", _heterane_local(z, v)) for z in _HETERANE_Z for v in _HETERANE_V
      if v > STANDARD_BONDING_NUMBERS[z]],
    *[("heterane", p) for z in HETERANE_CHAIN_Z for p in _heterane_chain_local(z)],
    *[("heterane", _standard_heterane_local(z)) for z in sorted(HYDRIDE_STD_Z)],
)


def _compile() -> dict[str, tuple[Chem.Mol, ...]]:
    """按 FG 键归并编译整表，SMARTS 写错时立即失败。"""
    out: dict[str, list[Chem.Mol]] = {}
    for fg, smarts in FG_SMARTS:
        pat = Chem.MolFromSmarts(smarts)
        if pat is None:
            raise ValueError("FG_SMARTS 无法解析: %s %s" % (fg, smarts))
        out.setdefault(fg, []).append(pat)
    return {fg: tuple(v) for fg, v in out.items()}


_PATTERNS_BY_FG = _compile()


def match_local_fg(mol: Chem.Mol) -> dict[str, list[tuple[int, ...]]]:
    """返回 {FG 键: 匹配元组升序列表}；元组首元素为中心原子索引。"""
    out: dict[str, list[tuple[int, ...]]] = {}
    for fg, pats in _PATTERNS_BY_FG.items():
        by_center: dict[int, tuple[int, ...]] = {}
        for pat in pats:
            for m in mol.GetSubstructMatches(pat):
                by_center.setdefault(m[0], m)
        out[fg] = [by_center[i] for i in sorted(by_center)]
    return out
