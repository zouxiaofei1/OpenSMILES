"""L1 官能团分析器：枚举分子中各类官能团条目并汇总为分析结果 dict。"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol
from collections import deque
from namepredict.tools import memo
from namepredict.constants import (
    C,P, CARBONYL_COMPOSITES, FG_PARTS_KEY, H, N, O, RING_HETERO, S,
)
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer1._carbonyl_common import (
    _alkoxy_c_of,
    _amide_n_info,
    _double_bonded_o_idxs,
    _ester_alkoxy_of as _ester_alkoxy_of_common,
    _has_acid_o_neighbor,
    _has_double_bonded_o,
    _is_anhydride_bridge_o,
    _is_single_c_oh,
)

def _heavy(mol: Mol, a) -> list:
    """非氢邻居索引。"""
    return [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]


def _arm_component(mol: Mol, start: int, core: set[int]) -> set[int]:
    """从臂根 start 取不穿过 core 的连通重原子组分。"""
    comp: set[int] = {start}
    dq: deque[int] = deque([start])
    while dq:
        i = dq.popleft()
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            if nb.GetAtomicNum() == 1:
                continue
            j = nb.GetIdx()
            if j in core or j in comp:
                continue
            comp.add(j)
            dq.append(j)
    return comp


def _one_phosphate(mol: Mol, p_idx: int) -> dict | None:
    """判定单个 P 是否为磷酸中心，是则返回 {p_idx,n_oh,n_om,n_arms}，否则 None。"""
    p = mol.GetAtomWithIdx(p_idx)
    if p.GetTotalNumHs() != 0 or p.GetFormalCharge() != 0:
        return None
    nei = _heavy(mol, p)
    if len(nei) != 4 or any(mol.GetAtomWithIdx(i).GetAtomicNum() != O for i in nei):
        return None
    single_o: list[int] = []
    dbl_o: list[int] = []
    for i in nei:
        bt = mol.GetBondBetweenAtoms(p_idx, i).GetBondType()
        if bt == BondType.DOUBLE:
            dbl_o.append(i)
        elif bt == BondType.SINGLE:
            single_o.append(i)
        else:
            return None
    if len(dbl_o) != 1 or len(single_o) != 3:
        return None
    da = mol.GetAtomWithIdx(dbl_o[0])
    if da.GetFormalCharge() != 0 or da.GetTotalNumHs() != 0 or _heavy(mol, da) != [p_idx]:
        return None
    core = {p_idx, *nei}
    n_oh = n_om = n_arms = 0
    arm_all: set[int] = set()
    for o_idx in single_o:
        a = mol.GetAtomWithIdx(o_idx)
        heavy = _heavy(mol, a)
        if a.GetFormalCharge() == 0 and a.GetTotalNumHs() >= 1 and heavy == [p_idx]:
            n_oh += 1
            continue
        if a.GetFormalCharge() == -1 and a.GetTotalNumHs() == 0 and heavy == [p_idx]:
            n_om += 1
            continue
        if a.GetFormalCharge() != 0 or a.GetTotalNumHs() != 0 or p_idx not in heavy:  # O–R：中性无 H、除 P 外另连 1 个重原子
            return None
        others = [j for j in heavy if j != p_idx]
        if len(others) != 1 or mol.GetAtomWithIdx(others[0]).GetAtomicNum() != C:
            return None
        comp = _arm_component(mol, others[0], core)
        attaches = [j for i in comp for j in _heavy(mol, mol.GetAtomWithIdx(i)) if j in core]  # 组分只贴 1 个 core 原子（桥 O）；不能连到 P 或其它 O
        if not attaches or len(set(attaches)) != 1 or attaches[0] != o_idx:
            return None
        arm_all |= comp
        n_arms += 1
    if n_oh + n_om + n_arms != 3:
        return None
    all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}  # 整分子纯度：重原子 = core ∪ 臂（排除臂间成环、P–O–P 焦磷酸等）
    if all_heavy != (core | arm_all):
        return None
    return {"p_idx": p_idx, "n_oh": n_oh, "n_om": n_om, "n_arms": n_arms}


def phosphate_entries(mol: Mol) -> list[dict]:
    """分子中全部磷酸中心（P(=O)(O)₃）的条目列表。"""
    result = [e for a in mol.GetAtoms() if a.GetAtomicNum() == P  and (e := _one_phosphate(mol, a.GetIdx())) is not None]
    return result
def _acyl_hal_of(carbon) -> tuple[int, int] | None:
    """惰性导入酰卤检测并返回碳上的卤素邻居信息。"""
    from namepredict.layer1.acyl_halide import acyl_hal_of
    return acyl_hal_of(carbon)

def _is_carboxyl_carbon(atom) -> bool:
    """判断碳是否为羧基碳（羰基双键氧 + 酸性氧邻居）。"""
    if atom.GetAtomicNum() != C:
        return False
    return _has_double_bonded_o(atom) and _has_acid_o_neighbor(atom)

def _carbon_neighbor_count(atom) -> int:
    """统计原子连有的碳邻居个数。"""
    return len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == C])

def _is_amide_carbon(atom) -> bool:
    """判断碳是否为酰胺羰基碳（排除酸、酯）。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom) or _ester_alkoxy_of(atom) is not None:
        return False
    return _amide_n_info(atom) is not None

def _has_ring_hetero_neighbor(atom) -> bool:
    """判断原子是否连有环内杂原子（N/O/S）：O/S 与羰基同环即经杂原子闭合成内酯/硫代内酯，环内 N 的酰基则按酮命名（P-66.1.1，N-酰基环胺），三者同由 _is_ketone_carbon 作环酮。"""
    return any(n.GetAtomicNum() in RING_HETERO and n.IsInRing() for n in atom.GetNeighbors())

def _is_ketone_carbon(atom) -> bool:
    """判断碳是否为酮羰基碳（非酸、非酰胺、非酯）；双碳邻居，或单碳邻居 + 环内杂原子（N-酰基环胺 → ethanone 型母体；环内 O/S → 内酯/硫代内酯按杂环 -one 命名），或环内零碳邻居（环脲/环碳酸酯型）。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom):
        return False
    n_c = _carbon_neighbor_count(atom)
    if n_c == 1:  # 环外单碳羰基须连环内杂原子（N-酰基环胺）才作酮；环内单碳羰基是内酰胺/环酮/内酯/硫代内酯，同样作酮。
        if not _has_ring_hetero_neighbor(atom) or _is_aldehyde_carbon(atom):
            return False
    elif n_c == 0:  # 环内零碳邻居羰基作环酮；环内非内酯型酯与开链者（脲/CO2）不作。
        if not atom.IsInRing() or (_ester_alkoxy_of(atom) is not None and not _is_lactone_carbon(atom)):
            return False
        return True
    elif n_c != 2:
        return False
    return True

def _anhydride_o_of(carbon) -> int | None:
    """返回羰基碳上的酸酐桥氧索引。"""
    for n in carbon.GetNeighbors():
        if _is_anhydride_bridge_o(n):
            return n.GetIdx()
    return None


def _is_ester_alkoxy_o(oxygen, carbonyl) -> bool:
    """判断 O 是否为酯烷氧基氧（排除酸酐桥）。"""
    if oxygen.GetAtomicNum() != O or oxygen.GetTotalNumHs() != 0:
        return False
    if _is_anhydride_bridge_o(oxygen):
        return False
    return _alkoxy_c_of(oxygen, carbonyl) is not None

def _ester_alkoxy_of(carbon) -> tuple[int, int] | None:
    """在碳上查找酯烷氧基侧并返回 (o_idx, alkoxy_c_idx)。"""
    return _ester_alkoxy_of_common(carbon, _is_ester_alkoxy_o)

def _is_lactone_carbon(atom) -> bool:
    """判断碳是否为环内酯（内酯）羰基碳：酯氧在环内时并入环母体作 -one 后缀（2H-chromen-2-one / 2-benzofuran-1-one / 1,3-dioxolan-2-one），不按酯的 -oate 命名。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    alkoxy = _ester_alkoxy_of(atom)
    if alkoxy is None:
        return False
    return atom.GetOwningMol().GetAtomWithIdx(alkoxy[0]).IsInRing()

def _is_ester_carbon(atom) -> bool:
    """判断碳是否为酯羰基碳（有酸性氧与烷氧基侧；环内酯除外，见 _is_lactone_carbon）。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom) or _ester_alkoxy_of(atom) is None:
        return False
    return not _is_lactone_carbon(atom)

def _ald_blocked(atom) -> bool:
    """判断醛碳是否被酯、酰卤、酸酐或酰胺占用。"""
    if _ester_alkoxy_of(atom) is not None or _acyl_hal_of(atom) is not None:
        return True
    return _anhydride_o_of(atom) is not None

def _is_aldehyde_carbon(atom) -> bool:
    """判断碳是否为醛羰基碳（单碳邻居、带 H 且未被阻断）；环内羰基不再作醛——内酰胺/环脲等由 _is_ketone_carbon 作酮、以 -one 后缀表达（P-66.6.1），当作醛会把喹唑啉-4-酮错拼成 …醛。"""
    if atom.GetAtomicNum() != C or atom.GetTotalDegree() < 3:
        return False
    if not _has_double_bonded_o(atom) or _has_acid_o_neighbor(atom):
        return False
    if _carbon_neighbor_count(atom) > 1:
        return False
    if atom.GetTotalNumHs() < 1:  # 无 H 的羰基是 N-酰基/环酮/内酰胺等，不作醛（环上外环 -CHO 带 H，照旧作醛）
        return False
    return not _ald_blocked(atom)

def _is_hydroxyl_oxygen(atom) -> bool:
    """判断 O 是否为醇羟基氧（排除羧酸羟基）。"""
    if not _is_single_c_oh(atom):
        return False
    if _is_carboxyl_carbon(_carbon_neighbor(atom)):
        return False
    return True

def _carbon_neighbor(atom):
    """返回原子连有的第一个碳邻居。"""
    return next(n for n in atom.GetNeighbors() if n.GetAtomicNum() == C)

def _hydroxyl_entry(atom) -> dict:
    """组装单个羟基条目 dict（氧为中心，所连碳为周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": [_carbon_neighbor(atom).GetIdx()]}

def _hydroxyl_entries(mol: Mol) -> list[dict]:
    """收集分子中所有羟基条目的列表。"""
    return [_hydroxyl_entry(a) for a in mol.GetAtoms() if _is_hydroxyl_oxygen(a)]

def _thiol_entries(mol: Mol) -> list[dict]:
    """收集分子中所有硫醇条目的列表。"""
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != S or atom.GetTotalNumHs() < 1:
            continue
        if _carbon_neighbor_count(atom) != 1:
            continue
        out.append({"center_idx": atom.GetIdx(), "surr_idx": [_carbon_neighbor(atom).GetIdx()]})
    return out

def _amine_degree(atom) -> int | None:
    """返回胺 N 取代度（1/2/3），非胺返回 None；排除芳香/酰胺/环内 N（吡咯烷等的环 N 是环杂原子，非胺官能团）。"""
    if atom.GetAtomicNum() != N  or atom.GetIsAromatic():
        return None
    if atom.IsInRing():
        return None
    n_c, n_h = _carbon_neighbor_count(atom), atom.GetTotalNumHs()
    if n_c == 1 and n_h >= 2:
        return 1
    return 2 if n_c == 2 and n_h == 1 else (3 if n_c == 3 and n_h == 0 else None)

def _amine_entry(atom) -> dict:
    """组装胺条目 dict（氮为中心，全部碳臂为周边）；取代度即碳臂数，由下游取 len(surr_idx)。"""
    return {"center_idx": atom.GetIdx(),
            "surr_idx": [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == C]}

def _amine_entries(mol: Mol) -> list[dict]:
    """收集分子中所有胺条目的列表。"""
    return [_amine_entry(a) for a in mol.GetAtoms() if _amine_degree(a) is not None]

def _carboxyl_entry(atom) -> dict:
    """组装单个羧基条目 dict（羧基碳为中心，两个氧为周边）；是否阴离子由下游按氧的形式电荷现判。"""
    return {"center_idx": atom.GetIdx(),
            "surr_idx": [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == O]}

def _carboxyl_entries(mol: Mol) -> list[dict]:
    """收集分子中所有羧基条目的列表。"""
    return [_carboxyl_entry(a) for a in mol.GetAtoms() if _is_carboxyl_carbon(a)]

def _carbonyl_entry(atom) -> dict:
    """组装单羰基条目 dict（羰基碳为中心，羰基氧为周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": _double_bonded_o_idxs(atom)}

def _ketone_entries(mol: Mol) -> list[dict]:
    """收集分子中所有酮条目的列表。"""
    return [_carbonyl_entry(a) for a in mol.GetAtoms() if _is_ketone_carbon(a)]

def _amide_entry(atom) -> dict:
    """组装单个酰胺条目 dict（羰基碳为中心，羰基氧与酰胺氮为周边）。"""
    n_idx, _ = _amide_n_info(atom)
    return {"center_idx": atom.GetIdx(), "surr_idx": [*_double_bonded_o_idxs(atom), n_idx]}

def _amide_entries(mol: Mol) -> list[dict]:
    """收集分子中所有酰胺条目的列表。"""
    return [_amide_entry(a) for a in mol.GetAtoms() if _is_amide_carbon(a)]

def _aldehyde_entries(mol: Mol) -> list[dict]:
    """收集分子中所有醛条目的列表。"""
    return [_carbonyl_entry(a) for a in mol.GetAtoms() if _is_aldehyde_carbon(a)]

def _acyl_chloride_entries(mol: Mol) -> list[dict]:
    """酰卤条目（F/Cl/Br/I）；键保留以兼容 L2/L3。"""
    from namepredict.layer1.acyl_halide import acyl_halide_entries
    return acyl_halide_entries(mol)

def _ester_entry(atom) -> dict:
    """组装单个酯条目 dict（羰基碳为中心，羰基氧与酯氧为周边）；烷氧基臂属取代基侧，不入 FG 本体。"""
    o_idx, _ = _ester_alkoxy_of(atom)
    return {"center_idx": atom.GetIdx(), "surr_idx": [*_double_bonded_o_idxs(atom), o_idx]}

def _ester_entries(mol: Mol) -> list[dict]:
    """收集分子中所有酯条目的列表。"""
    return [_ester_entry(a) for a in mol.GetAtoms() if _is_ester_carbon(a)]

def _anhydride_entry(atom) -> dict:
    """组装单个酸酐条目 dict（桥氧为中心，两个羰基碳与两个羰基氧为周边），并归一化两碳顺序。"""
    mol = atom.GetOwningMol()
    o_idx = _anhydride_o_of(atom)
    c1, c2 = sorted((atom.GetIdx(), True))
    return {"center_idx": o_idx, "surr_idx": [c1, c2, *_double_bonded_o_idxs(atom),
                                              *_double_bonded_o_idxs(mol.GetAtomWithIdx(other))]}

def _is_anhydride_carbon(atom) -> bool:
    """判断碳是否为酸酐羰基碳（排除羧酸）。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom):
        return False
    if _has_acid_o_neighbor(atom):
        return False
    return _anhydride_o_of(atom) is not None

def _anhydride_entries(mol: Mol) -> list[dict]:
    """收集分子中所有去重后的酸酐条目列表。"""
    seen: set[int] = set()
    out: list[dict] = []
    for atom in mol.GetAtoms():
        if not _is_anhydride_carbon(atom):
            continue
        e = _anhydride_entry(atom)
    return out

def _is_cc_double(bond) -> bool:
    """判断键是否为 C=C 双键（排除芳香键）。"""
    if bond.GetBondType() != BondType.DOUBLE or bond.GetIsAromatic():
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == C and b.GetAtomicNum() == C

def _is_cc_triple(bond) -> bool:
    """判断键是否为 C≡C 三键。"""
    if bond.GetBondType() != BondType.TRIPLE:
        return False
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    return a.GetAtomicNum() == C and b.GetAtomicNum() == C

def _bond_entry(bond) -> dict:
    """将 C=C/C≡C 键组装为条目 dict（两碳索引有序）。"""
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    return {"c1": min(a, b), "c2": max(a, b)}

def _filter_bond_entries(mol: Mol, pred, entry_fn) -> list[dict]:
    """按谓词过滤键并映射为条目列表。"""
    return [entry_fn(bond) for bond in mol.GetBonds() if pred(bond)]

def _double_bond_entries(mol: Mol) -> list[dict]:
    """收集分子中所有 C=C 双键条目。"""
    return _filter_bond_entries(mol, _is_cc_double, _bond_entry)

def _triple_bond_entries(mol: Mol) -> list[dict]:
    """收集分子中所有 C≡C 三键条目。"""
    return _filter_bond_entries(mol, _is_cc_triple, _bond_entry)

def _is_cn_triple(bond) -> bool:
    """判断键是否为 C≡N 三键。"""
    if bond.GetBondType() != BondType.TRIPLE:
        return False
    z = {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
    return z == {6, 7}

def _nitrile_entry(bond) -> dict:
    """将腈三键组装为条目 dict（腈碳为中心，氮为周边）。"""
    a, b = bond.GetBeginAtom(), bond.GetEndAtom()
    c = a if a.GetAtomicNum() == C else b
    n = b if a.GetAtomicNum() == C else a
    return {"center_idx": c.GetIdx(), "surr_idx": [n.GetIdx()]}

def _nitrile_entries(mol: Mol) -> list[dict]:
    """收集分子中所有腈条目的列表。"""
    return _filter_bond_entries(mol, _is_cn_triple, _nitrile_entry)

def _ring_entry(atom_ids: tuple) -> dict:
    """将环原子索引元组组装为环条目 dict。"""
    return {"atom_ids": atom_ids}

def _ring_entries(mol: Mol) -> list[dict]:
    """收集分子中所有环条目的列表；环感知对分子恒定，按 mol 记忆（同 ring_systems._sssr）。"""
    return memo.by_mol("ring_entries", lambda m: [_ring_entry(r) for r in m.GetRingInfo().AtomRings()], mol)

def _ring_meta(mol: Mol) -> dict:
    """汇总环事实：环条目、环系与数量统计。"""
    from namepredict.layer1.ring_systems import build_ring_systems
    rings = _ring_entries(mol)
    systems = build_ring_systems(mol)
    return {
        "rings": rings, "n_rings": len(rings), "has_ring": bool(rings),
        "ring_systems": systems, "n_ring_systems": len(systems),
    }

def _carbon_ids(mol: Mol) -> list[int]:
    """返回分子中所有碳原子的索引列表。"""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == C]

def _is_anchored(atom) -> bool:
    """原子是否直接连着 `*` 自由基锚点（原子序 0 邻居）。"""
    return any(n.GetAtomicNum() == 0 for n in atom.GetNeighbors())


def _is_acyl_head(mol: Mol, atom) -> bool:
    """锚定羰基碳是否为酰基头：带 =O、恰好 1 个单键碳邻居、无其它重邻居（P-65.1.7.2 酸衍生）；环酮/内酯/酰胺因双碳或杂原子邻居被排除。"""
    if atom.GetAtomicNum() != C or not _has_double_bonded_o(atom) or not _is_anchored(atom):
        return False
    carbs: list = []
    hetero = False
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() == 0:
            continue
        b = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
        if b is not None and b.GetBondType() == BondType.DOUBLE:
            continue  # 羰基 =O（或醛 C=C 等其它双键；真实酰基头无第二个双键）
        if nb.GetAtomicNum() == C:
            carbs.append(nb)
        elif nb.GetAtomicNum() != H:
            hetero = True
    if hetero or len(carbs) != 1:
        return False
    return True


def _acyl_entries(mol: Mol) -> list[dict]:
    """虚拟原子邻居中判为酰基头的羰基碳条目（取代 radical 成为主基团，避免醛→酮误降级）。"""
    return [_carbonyl_entry(n)
            for a in mol.GetAtoms() if a.GetAtomicNum() == 0
            for n in a.GetNeighbors() if n.GetAtomicNum() != H and _is_acyl_head(mol, n)]


def _radical_entries(mol: Mol, exclude: frozenset[int] = frozenset()) -> list[dict]:
    """虚拟原子（原子序 0）邻居碳的 P-41 自由基位点；"""
    out: list[dict] = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 0:
            continue
        for n in a.GetNeighbors():
            if n.GetAtomicNum() != H and n.GetIdx() not in exclude:
                out.append({"center_idx": n.GetIdx(), "surr_idx": []})
    return out


def _bond_lists(mol: Mol) -> dict:
    """不饱和度键列表：C=C / C≡C 是结构事实而非官能团，不入 FG 通道（_fg_parts 只承载 P-41 官能团条目）。"""
    return {"double_bonds": _double_bond_entries(mol), "triple_bonds": _triple_bond_entries(mol)}

_SUPPRESSIBLE = {**CARBONYL_COMPOSITES, "nitriles": "nitrile"}  # 可被更高优先级 FG 整体压制的组合 FG：组合羰基 + 腈
_LEAF_DEMOTED = ("carboxyls", "nitriles")  # 降级为"前缀叶"的组合 FG：整组碳排除出主链（P-61.1.3 carboxy/cyano）。其余组合 FG（酯/酰胺/醛/酰卤/酸酐）降级为"氧代"——羰基碳留在链内，仅 O 作 oxo/formyl 前缀，由 L3 锚定叶识别。


def _arbitrate_parts(parts: dict) -> tuple[dict, frozenset[str]]:
    """P-41 主基团仲裁：已有更高优先级 FG 时，组合 FG 整组退出主基团。叶型降级（羧酸/腈）的条目保留在清单中并标记 demoted（供 L2 排除出主链），其余组合 FG 整组清空。返回 (仲裁后 parts, 降级的 occurrence id 集)。"""
    p41 = {sp.fg: sp.p41 for sp in FG_SPECS if sp.p41}
    present = {fg for fg, key in FG_PARTS_KEY.items() if parts.get(key)}
    out = dict(parts)
    demoted: set[str] = set()
    for key, fg in _SUPPRESSIBLE.items():
        if out[key] and any(p41[h] < p41[fg] for h in present if h != fg):
            if key in _LEAF_DEMOTED:
                demoted |= {f"{key}:{i}" for i in range(len(out[key]))}
            else:
                out[key] = []
    return out, frozenset(demoted)


def _detect_parts(mol: Mol) -> dict:
    """检测（未仲裁）分子中各类官能团条目；键集由 fg_registry.FG_SPECS 的 list_key 派生，唯一事实来源。"""
    acyls = _acyl_entries(mol)
    heads = frozenset(e["center_idx"] for e in acyls)
    return {"carboxyls": _carboxyl_entries(mol), "hydroxyls": _hydroxyl_entries(mol),
        "esters": _ester_entries(mol), "amides": _amide_entries(mol),
        "ketones": _ketone_entries(mol), "radicals": _radical_entries(mol, heads),
        "acyls": acyls,
        "aldehydes": [e for e in _aldehyde_entries(mol) if e["center_idx"] not in heads],
        "amines": _amine_entries(mol),
        "nitriles": _nitrile_entries(mol),
        "acyl_chlorides": _acyl_chloride_entries(mol),
        "anhydrides": _anhydride_entries(mol), "thiols": _thiol_entries(mol),
        "phosphates": phosphate_entries(mol)}

def _collect_fgs(mol: Mol) -> dict:
    """聚合官能团条目并构建带类型清单（单一 FG 出口：清单承载全部 FG 事实，存在性由清单内容判定；不饱和度独立于 FG 通道）。"""
    from namepredict.layer1.functional_group_inventory import build_inventory

    parts, demoted = _arbitrate_parts(_detect_parts(mol))
    return {**_bond_lists(mol), "fg_inventory": build_inventory(parts, mol, demoted)}

def _info(mol: Mol, carbons: list[int]) -> dict:
    """组装分子分析结果 dict（碳信息 + 官能团 + 环事实）。"""
    base = {"mol": mol, "carbon_ids": carbons, "n_carbons": len(carbons)}
    return {**base, **_collect_fgs(mol), **_ring_meta(mol)}

def analyze(mol: Mol) -> dict:
    """分析分子并返回完整的官能团与结构信息 dict。"""
    result = _info(mol, _carbon_ids(mol))
    # print(result,"\n\n\n")
    return result
