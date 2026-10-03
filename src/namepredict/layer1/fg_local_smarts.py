"""L1 官能团的局部环境 SMARTS 表（唯一事实来源）。
每条 = (FG 键, SMARTS)，模式首原子即该官能团的中心原子，同名多条取并集。
FG 键与 fg_registry.FG_SPECS 一致；新增官能团只加表项，不改检测代码。
非局部判据（酰基头的 heads 联动、磷酸的臂回接与整分子纯度）留在 analyzer 后置。
"""
from __future__ import annotations

from rdkit import Chem

from namepredict.constants import B, Br, C, Cl, F, I, N, O, P, Po, S, Se, STANDARD_BONDING_NUMBERS, Te

# 酸性氧：与单一碳相连的羟基氧，或羧酸盐阴离子氧
_ACID_O = "[$([#8;H1]),$([#8;-1;X1;H0])]"
# 非羧酸碳
_NOT_ACID = "!$([#6;X3](=[#8;X1])%s)" % _ACID_O
# 非环酯烷氧基氧（内酯的环内氧不算，留给酮分支）
_NOT_ACYCLIC_ESTER = "!$([#6;X3]~[#8;X2;H0;+0;!R]~[#6])"
# 无卤素邻居
_NOT_HALO = "!$([#6;X3]~[#9,#17,#35,#53])"
# 至多一个碳邻居
_ONE_C = "!$([#6;X3](~[#6])~[#6])"
# 环内杂原子邻居（内酰胺/内酯/硫代内酯、N-酰基环胺）
_RING_HET = "~[#7,#8,#16;R]"
# 含氧酸中心（P/S）的单键氧三态：OH / O⁻ / O-R（臂根限碳或磷，P-O-P 桥氧即多聚磷酸；排除哑原子臂）
_O_PHOS = "(-[$([#8;X2;H1;+0]),$([#8;X1;-1;H0]),$([#8;X2;H0;+0]~[#6,#15])])"
_O_ARM = "(-[#8;X2;H0;+0]~[#6])"      # O-R 臂（硫酸酯/磺酸酯的 O 侧）
_O_S_ARM = "(-[#8;X2;H0;+0]~[#16])"   # S-O-S 桥臂（多硫酸链；P-O-P 桥见 _O_PHOS 的磷臂）
_HALO_ARM = "(-[#9,#17,#35,#53])"     # 卤素臂（磺酰卤）
_N_ARM = "(-[#7;!R])"                 # 非环氮臂（磺酰胺；环内 S-N 走磺酰前缀）
_C_ARM = "(-[#6])"                    # 直连碳臂（碳骨架成员）
_SNY_ARM = "(-[#6,#7,#16])"           # 膦酸的第三臂（C/S/N 均可）

# 单核母体阳离子的中心元素（P-73.1.1.1 表 7.3 的 15/16/17 族）
CATION_Z = (N, P, O, S, F, Cl, Br, I)


def _cation_local(z: int) -> str:
    """单核阳离子中心的局部模式：非环、半径 1 内无负形式电荷原子。"""
    # 负电荷邻居闸排除硝基/叠氮/N-氧化物/异氰/高氯酸根等"阳离子寄居在别的基团里"的结构
    return "[#%d;+1;!R;!$([#%d;+1]~[-1])]" % (z, z)


# 非碳母体氢化物（P-21 / P-14.1）：中性杂原子且键数偏离标准值。标准价的自任母体情形
# （如膦 P(III)、硫醚 S(II)）由 P-41 类 41 排在碳之后，不在此列，故只取 v > 标准值。
_HETERANE_V = (3, 4, 5, 6, 7)  # 非标准键数上界；v 为 RDKit 总价（含氢），即 P-14.1.1 的键数
_HETERANE_Z = tuple(z for z in STANDARD_BONDING_NUMBERS if z != C)


HETERANE_CHAIN_Z = tuple(z for z in STANDARD_BONDING_NUMBERS if z not in (C, B))  # 硼烷按 P-21.2.2 排除


def _oxo_guard(z: int) -> str:
    """排除带 =O 的中心（归含氧酸通路）；非硫元素另排除 =S（归磷硫酰通路），硫自身的 S=S 属链不饱和。"""
    bans = ("[#8]",) if z == 16 else ("[#8]", "[#16]")
    return "".join(f";!$([#{z}]={b})" for b in bans)


def _heterane_local(z: int, v: int) -> str:
    """杂原子烃中心的局部模式：中性、键数为 v、不带 oxo/thioxo。"""
    return f"[#{z};v{v};+0{_oxo_guard(z)}]"


HETERANE_MIN_CHAIN = {N: 3, O: 3, S: 3, Se: 3, Te: 3, Po: 3}  # 取代链须 ≥3 连的元素：硫族 1–2 连续原子不作母体氢化物（P-68.4.0），N 参照 P-68.3.1.4.1 的 triazane 起
HETERANE_CHAIN_MIN_DEFAULT = 2  # 其余元素（Si/Ge/Sn/Pb、P/As/Sb/Bi 等）：≥2 连即可带取代


def _is_parent_hydride(mol, comp: set[int], z: int) -> bool:
    """该同元素集合是否为「由氢饱和」的母体氢化物（除链内原子外只连氢，P-21.2.2）。"""
    for i in comp:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetAtomicNum() not in (1, z):
                return False
    return True


def heterane_chain_ok(mol, idx: int) -> bool:
    """该原子能否作杂原子链中心：单原子认非标准键数；多原子按元素最小链长（P-21.2.2 / P-68.4.0）。

    取代链须达元素最小链长（硫族与 N 为 3 连），纯母体氢化物只需 2 连。
    """
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
        return False  # 标准键数的孤立杂原子不自任母体（硫醚/醚等归 P-41 类 41）
    if _is_parent_hydride(mol, comp, z):
        return True
    return len(comp) >= HETERANE_MIN_CHAIN.get(z, HETERANE_CHAIN_MIN_DEFAULT)


def _heterane_chain_local(z: int) -> tuple[str, ...]:
    """杂原子链中心（P-21.2.2）：存在同元素链即命中；取代链的最小链长按元素分档。

    写成单原子递归查询而非 `A~A`：对称双原子模式会被 GetSubstructMatches 的 uniquify 折叠成一条。
    """
    base = f"[#{z};+0;!R{_oxo_guard(z)}"  # !R：P-21.2.2 仅限无环链，环内杂原子走 Hantzsch-Widman/环系路径
    nz = f"[#{z};+0]"  # 链上邻居也须中性：叠氮 N=[N+]=[N-] 一类带电链不是母体氢化物
    pats = [f"{base};!$([#{z}]~[!#1;!#{z}]);$([#{z}]~{nz})]"]      # 纯母体氢化物：≥2 连
    pats.append(f"{base};$([#{z}]~{nz}~{nz})]" if z in HETERANE_MIN_CHAIN
                else f"{base};$([#{z}]~{nz})]")                     # 带取代/支链：≥3 连或 ≥2 连
    return tuple(pats)


FG_SMARTS: tuple[tuple[str, str], ...] = (
    # 自由基：哑原子的重原子邻居（中心为该重原子，非哑原子）
    ("radical", "[!#1]~[#0]"),
    # 酰基残基：哑原子所连的酰基头羰基碳；双键邻居不参与杂原子判定
    ("acyl", "[#6;X3;!$([#6]-[!#6;!#0;!#1])](=[#8;X1])(-[#6])~[#0]"),
    ("acid", f"[#6;X3](=[#8;X1]){_ACID_O}"),
    ("ester", f"[#6;X3;{_NOT_ACID}](=[#8;X1])[#8;X2;H0;+0;!R]~[#6]"),
    # P-65.6.3.3.7.1 硫代羧酸 S-酯（thioate）：酯氧换为 S，同属 P-41 类 9，词尾由 L5 按 S 侧切换
    ("ester", f"[#6;X3;{_NOT_ACID}](=[#8;X1])[#16;X2;H0;+0;!R]~[#6]"),
    ("acyl_halide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#8;X1])~[#9,#17,#35,#53]"),
    # 酰胺 N 除羰基碳/H 外只容 C 或不再连碳的 O（哑原子、杂原子均拒）
    # N 允许 X2：-C(=O)-N=C< 的亚胺型酰胺 N（gold 取 carboxamide 后缀 + N-亚基前缀，P-66.1.1.1.1.3）
    ("amide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#8;X1])-[#7;X2,X3;!R;!$([#7;X2,X3]~[!#6;!#8]);!$([#7;X3]~[#8;X2]~[#6])]"),
    # P-43 类 16/17：硫代酰胺 (C=S) 与脒 (C=N) 与酰胺同组，词尾由 L5 按双键杂原子切换
    ("amide", f"[#6;X3;{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#16;X1])-[#7;X2,X3;!R;!$([#7;X2,X3]~[!#6;!#8]);!$([#7;X3]~[#8;X2]~[#6])]"),
    ("amide", f"[#6;X3;!R;!$([#6](~[#7])(~[#7])~[#7]);{_NOT_ACID};{_NOT_ACYCLIC_ESTER}](=[#7;X2,X3])-[#7;X2,X3;!R;!$([#7;X2,X3]~[!#6;!#8]);!$([#7;X3]~[#8;X2]~[#6])]"),
    ("aldehyde", f"[#6;X3;H1,H2;{_ONE_C};{_NOT_ACID};{_NOT_ACYCLIC_ESTER};{_NOT_HALO}](=[#8;X1])"),
    ("nitrile", "[#6;+0;X2]#[#7;X1;+0]"),
    # 酮按碳邻居数分三条：两个碳 / 单碳须连环内杂原子（P-66.1.1）/ 环内零碳
    # 环内单碳羰基是内酰胺/环酮/内酯/硫代内酯，同样作酮；环内零碳者只排除非内酯型酯
    ("ketone", f"[#6;X3;{_NOT_ACID}](=[#8;X1])(~[#6])~[#6]"),
    ("ketone", f"[#6;X3;{_NOT_ACID};{_ONE_C}](=[#8;X1])(~[#6]){_RING_HET}"),
    ("ketone", f"[#6;X3;{_NOT_ACID};R;{_ONE_C};{_NOT_ACYCLIC_ESTER}](=[#8;X1])~[#7,#8,#16]"),
    # 含氧酸：中心 P/S 中性无氢、恰一个/两个末端双键氧，余下为酸式氧或有机臂
    # 每种"元素+双键氧数+臂型"一条；oxo_kind 由 analyzer 依臂角色归一（同一来源）
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
    ("thiol", "[#16;X2;H1]~[#6]"),
    # 胺三条对应原 _amine_degree 的取代度分支；键级不限（~），与 _carbon_neighbor_count 一致
    # 中心限 +0：铵/亚胺鎓的 N+ 是母体阳离子（P-73.1.1），不是胺
    ("amine", "[#7;!R;!a;+0;H2,H3,H4;!$([#7](~[#6])~[#6])]~[#6]"),
    ("amine", "[#7;!R;!a;+0;H1;!$([#7](~[#6])(~[#6])~[#6])](~[#6])~[#6]"),
    ("amine", "[#7;!R;!a;+0;H0;!$([#7](~[#6])(~[#6])(~[#6])~[#6])](~[#6])(~[#6])~[#6]"),
    # 单核母体阳离子（P-73.1.1）：每元素一条，中心即阳离子原子
    *[("cation", _cation_local(z)) for z in CATION_Z],
    # 非碳母体氢化物（P-21）：每 (元素, 非标准键数) 一条，中心即杂原子
    *[("heterane", _heterane_local(z, v)) for z in _HETERANE_Z for v in _HETERANE_V
      if v > STANDARD_BONDING_NUMBERS[z]],
    # 杂原子链母体（P-21.2.2）：每元素一条（或多条分档），命中即该原子属于一条均一杂原子链
    *[("heterane", p) for z in HETERANE_CHAIN_Z for p in _heterane_chain_local(z)],
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
