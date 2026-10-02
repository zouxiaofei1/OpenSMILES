"""L1 官能团分析器：枚举分子中各类官能团条目并汇总为分析结果 dict。
含氧酸中心（P/S 非碳骨架成键）的非局部判据与 oxo_kind 归一本模块。
羰基原语 _double_bonded_o_idxs/_alkoxy_c_of 供 L2 主基团表达式复用。
"""
from __future__ import annotations

from rdkit.Chem import BondType, Mol
from collections import deque
from namepredict.constants import (
    C, H, N, O, OXO_CENTER_KINDS, S,
)
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer1.fg_local_smarts import match_local_fg

def _heavy(a) -> list:
    """非氢邻居索引。"""
    return [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() != 1]


def _double_bonded_o_idxs(carbon) -> list[int]:
    """返回碳上羰基双键氧的索引列表。"""
    mol = carbon.GetOwningMol()
    out: list[int] = []
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() != O:
            continue
        b = mol.GetBondBetweenAtoms(carbon.GetIdx(), n.GetIdx())
        if b is not None and b.GetBondType() == BondType.DOUBLE:
            out.append(n.GetIdx())
    return out


def _alkoxy_c_of(oxygen, carbonyl) -> int | None:
    """返回氧上除羰基碳外的烷氧基碳索引。"""
    for n in oxygen.GetNeighbors():
        if n.GetAtomicNum() == C and n.GetIdx() != carbonyl.GetIdx():
            return n.GetIdx()
    return None


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


_OXO_KIND_P = {(1, False): "phosphate", (1, True): "phosphonate"}  # P：(双键氧数, 有无直连碳臂) → kind
_OXO_KIND_S_ARM = {"halo": "sulfonyl_chloride", "n": "sulfonamide",   # S 带碳臂时：另一臂角色 → kind
                   "acid": "sulfonic", "o_arm": "sulfonate"}


def _oxo_partition(mol: Mol, z_idx: int) -> dict[str, list[int]]:
    """中心原子邻居按角色分桶：双键氧/羟基氧/阴离子氧/O-臂氧与臂根碳/碳/卤素/氮/硫。

    o_arm 存 O-臂的氧原子（计数用），o_arm_root 存臂根碳（_arm_single_attach 用）——
    两者个数不同，勿只看长度。
    """
    part: dict[str, list[int]] = {k: [] for k in
                                  ("oxo", "oh", "om", "o_arm", "o_arm_root", "c", "hal", "n", "s")}
    for n in mol.GetAtomWithIdx(z_idx).GetNeighbors():
        z, i = n.GetAtomicNum(), n.GetIdx()
        if z in (1, 0):
            continue
        if z == O:
            if mol.GetBondBetweenAtoms(z_idx, i).GetBondType() == BondType.DOUBLE:
                part["oxo"].append(i)
            elif n.GetFormalCharge() == -1:
                part["om"].append(i)
            elif n.GetTotalNumHs() >= 1:
                part["oh"].append(i)
            else:
                part["o_arm"].append(i)
                part["o_arm_root"].extend(j for j in _heavy(n) if j != z_idx)
        elif z == C:
            part["c"].append(i)
        elif z in (9, 17, 35, 53):
            part["hal"].append(i)
        elif z == N:
            part["n"].append(i)
        elif z == S:
            part["s"].append(i)
    return part


def _oxo_kind(mol: Mol, z_idx: int, part: dict | None = None) -> str | None:
    """按中心元素、双键氧数与臂角色归一 oxo_kind（B/P/S 同一张表驱动）。"""
    z = mol.GetAtomWithIdx(z_idx).GetAtomicNum()
    if z not in (5, 15, 16):
        return None
    r = {k: len(v) for k, v in (part or _oxo_partition(mol, z_idx)).items()}
    if z == 5:  # 硼酸 B(OH)2：恰一个碳臂 + 彻底酸式的两个氧（P-68.2.1）；硼酸酯（O-臂）不在此列
        return "boronic" if r["c"] == 1 and r["oh"] + r["om"] == 2 else None
    if z == 15:
        return _OXO_KIND_P.get((r["oxo"], r["c"] > 0))
    if r["oxo"] != 2:
        return None
    if not r["c"]:  # 无直连碳：中心自任母体，两臂皆酸式氧或 O-臂（硫酸/硫酸根/硫酸酯/多硫酸链）
        return "sulfate" if r["oh"] + r["om"] + r["o_arm"] == 2 else None
    arm = ("halo" if r["hal"] else "n" if r["n"] else "acid" if r["oh"] + r["om"] else "o_arm")
    return _OXO_KIND_S_ARM.get(arm)


def _arm_single_attach(mol: Mol, roots: list[int], core: set[int]) -> set[int] | None:
    """臂须单点回接且各臂互不相连：返回臂原子并集，违反返回 None。"""
    arm: set[int] = set()
    for root in roots:
        comp = _arm_component(mol, root, core)
        if comp & arm:
            return None
        arm |= comp
    return arm


def _oxoacid_entry(mol: Mol, core: tuple[int, ...]) -> dict | None:
    """校验含氧酸候选的非局部条件并归一 oxo_kind、锚点与负载。"""
    z_idx = core[0]
    part = _oxo_partition(mol, z_idx)  # 一次分桶供 kind 判定与臂根取用
    kind = _oxo_kind(mol, z_idx, part)
    if kind is None:
        return None
    core_set = set(core)
    c_roots, o_roots, oh, om = part["c"], part["o_arm_root"], part["oh"], part["om"]
    if kind in OXO_CENTER_KINDS:  # 中心为母体：臂单点回接 + 整分子纯度（排除臂间成环）
        arm = _arm_single_attach(mol, c_roots + o_roots, core_set)
        if arm is None:
            return None
        all_heavy = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1}
        if all_heavy != (core_set | arm):
            return None
        return {"oxo_z": z_idx, "center_idx": z_idx, "oxo_kind": kind,
                "n_oh": len(oh), "n_om": len(om)}
    if len(c_roots) != 1:  # 碳锚定：直连碳臂唯一，该碳即母体碳骨架的入口
        return None
    return {"oxo_z": z_idx, "center_idx": c_roots[0], "oxo_kind": kind, "n_oh": len(oh), "n_om": len(om),
            "surr_idx": list(om)}  # surr_idx 只收阴离子氧：L2 据此补 anion 标志


# oxo_kind → P-41 类别键：磺酰胺属酰胺类 11，其余落 oxoacid 类 9
_OXO_CLASS_BY_KIND = {"sulfonamide": "sulfonamide"}


def _oxo_bridge_arms(mol: Mol, z_idx: int) -> int:
    """中心 P/S 的 O-桥氧数（P-67.2 双核/多核磷酸与多硫酸的链内链节数；非 P/S 中心恒 0）。"""
    z = mol.GetAtomWithIdx(z_idx)
    znum = z.GetAtomicNum()
    if znum not in (15, 16):
        return 0
    n = 0
    for nb in z.GetNeighbors():
        if nb.GetAtomicNum() != O:
            continue
        if mol.GetBondBetweenAtoms(z_idx, nb.GetIdx()).GetBondType() == BondType.DOUBLE:
            continue
        if any(x.GetAtomicNum() == znum and x.GetIdx() != z_idx for x in nb.GetNeighbors()):
            n += 1
    return n


def oxoacid_entries(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> list[dict]:
    """分子中全部含氧酸中心（P/S 非碳骨架成键）的条目列表。"""
    if matches is None:
        matches = match_local_fg(mol).get("oxoacid", [])
    by_z: dict[int, dict] = {}
    for m in matches:
        e = _oxoacid_entry(mol, m)
        if e is not None:
            by_z[e["oxo_z"]] = e
    entries = [by_z[k] for k in sorted(by_z)]
    arms = {e["oxo_z"]: _oxo_bridge_arms(mol, e["oxo_z"]) for e in entries}  # 桥臂数每条只算一次
    # 缩合含氧酸（P-O-P / S-O-S）须链上留有全酸式末端（同一中心 ≥2 个酸式氧）才按功能母体识别（P-67.2.1）
    if not any(int(e["n_oh"]) + int(e["n_om"]) >= 2 for e in entries):
        # 无全酸式末端时也不能整链退为取代基：至少留酸式氧最多的链节作母体（否则母体落到甲烷）
        top_acid = max((int(e["n_oh"]) + int(e["n_om"]) for e in entries), default=0)
        entries = [e for e in entries
                   if not arms[e["oxo_z"]]
                   or int(e["n_oh"]) + int(e["n_om"]) == top_acid]
    # 链内中心让位于更少质子化的酸中心（P-41 酸根优先；否则抢走母体会把酸根写成前缀）
    top_om = max((int(e["n_om"]) for e in entries), default=0)
    return [e for e in entries
            if not (int(e["n_om"]) < top_om and arms[e["oxo_z"]])]


def boronic_entries(mol: Mol) -> list[dict]:
    """硼酸中心条目：SMARTS 表未登记 B，按元素直接扫描（B 无局部双键氧可匹配）。"""
    out: list[dict] = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 5:
            continue
        core = (a.GetIdx(),) + tuple(sorted(n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == O))
        e = _oxoacid_entry(mol, core)
        if e is not None:
            out.append(e)
    return out


def oxoacid_lists(mol: Mol, matches: list[tuple[int, ...]] | None = None) -> dict[str, list[dict]]:
    """含氧酸条目按 P-41 类别键归位（同一检测器，类别由 oxo_kind 表决定）。"""
    out: dict[str, list[dict]] = {}
    for e in oxoacid_entries(mol, matches) + boronic_entries(mol):
        out.setdefault(_OXO_CLASS_BY_KIND.get(e["oxo_kind"], "oxoacid"), []).append(e)
    return out

def _surr_idx(atom) -> list[int]:
    """周边原子：中心全部重原子邻居；碳中心不在环内时排除环内邻居。"""
    ring_excl = atom.GetAtomicNum() == C and not atom.IsInRing()
    return [n.GetIdx() for n in atom.GetNeighbors()
            if n.GetAtomicNum() != H and not (ring_excl and n.IsInRing())]

def _fg_entry(atom) -> dict:
    """组装官能团条目 dict（中心原子 + 重原子周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": _surr_idx(atom)}

def _cc_bond_entries(mol: Mol) -> tuple[list[dict], list[dict]]:
    """一遍遍历全部键，按 C=C（排除芳香）与 C≡C 分别组装条目（两碳索引有序）。"""
    dbl: list[dict] = []
    tpl: list[dict] = []
    for b in mol.GetBonds():
        bt = b.GetBondType()
        if bt == BondType.DOUBLE and not b.GetIsAromatic():
            dst = dbl
        elif bt == BondType.TRIPLE:
            dst = tpl
        else:
            continue
        if b.GetBeginAtom().GetAtomicNum() != C or b.GetEndAtom().GetAtomicNum() != C:
            continue
        a, z = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        dst.append({"c1": min(a, z), "c2": max(a, z)})
    return dbl, tpl

def _ring_meta(mol: Mol) -> dict:
    """汇总环事实：环条目、环系与数量统计。"""
    from namepredict.layer1.ring_systems import build_ring_systems, sssr_rings
    rings = [{"atom_ids": r} for r in sssr_rings(mol)]
    systems = build_ring_systems(mol)
    return {
        "rings": rings, "n_rings": len(rings), "has_ring": bool(rings),
        "ring_systems": systems, "n_ring_systems": len(systems),
    }

def _radical_entry(atom) -> dict:
    """组装自由基位点条目 dict（哑原子所连重原子为中心，无周边）。"""
    return {"center_idx": atom.GetIdx(), "surr_idx": []}


_SUPPRESSIBLE = frozenset({"acid", "ester", "acyl_halide", "amide", "nitrile", "aldehyde"})  # 可被更高优先级 FG 整体压制的组合 FG：组合羰基 + 腈
_PRESENCE_SKIP = frozenset({"oxoacid", "sulfonamide"})  # 含氧酸不参与存在性判定（纳入会改写压制结果）
_LEAF_DEMOTED = ("acid", "nitrile")  # 降级为"前缀叶"的组合 FG：整组碳排除出主链（P-61.1.3 carboxy/cyano）。其余组合 FG（酯/酰胺/醛/酰卤）降级为"氧代"——羰基碳留在链内，仅 O 作 oxo/formyl 前缀，由 L3 锚定叶识别。


def _has_negative_atom(mol: Mol) -> bool:
    """分子内是否存在形式电荷为负的原子（P-41 表 4.1 类 4 阴离子）。"""
    return mol is not None and any(a.GetFormalCharge() < 0 for a in mol.GetAtoms())


def _arbitrate_parts(parts: dict, mol: Mol | None = None,
                     has_anion: bool | None = None) -> tuple[dict, frozenset[str]]:
    """P-41 仲裁：更高优先级 FG 使组合 FG 退出，叶型标 demoted。"""
    p41 = {sp.fg: sp.p41 for sp in FG_SPECS if sp.p41}
    present = {fg for fg in p41 if parts.get(fg) and fg not in _PRESENCE_SKIP}
    if has_anion is None:
        has_anion = _has_negative_atom(mol)
    if has_anion:
        present.discard("cation")  # 阴离子（类 4）> 阳离子（类 6）：酸根在场时阳离子不压制酸
    out = dict(parts)
    demoted: set[str] = set()
    for fg in _SUPPRESSIBLE:
        if out[fg] and any(p41[h] < p41[fg] for h in present if h != fg):
            if fg in _LEAF_DEMOTED:
                demoted |= {f"{fg}:{i}" for i in range(len(out[fg]))}
            else:
                out[fg] = []
    return out, frozenset(demoted)


# 局部环境 SMARTS 命中且无非局部判据的 FG 键；顺序沿用历史输出
_LOCAL_ENTRY_FGS = ("acid", "alcohol", "ester", "amide", "ketone", "amine", "thiol",
                    "nitrile", "acyl_halide", "cation")


def _local_entries(mol: Mol, hits: dict) -> dict:
    """按局部环境命中结果组装各 FG 的条目列表（中心原子升序）。"""
    return {fg: [_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get(fg, [])]
            for fg in _LOCAL_ENTRY_FGS}


def _detect_parts(mol: Mol) -> dict:
    """检测（未仲裁）分子中各类官能团条目。"""
    hits = match_local_fg(mol)
    acyl = [_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("acyl", [])]  # 须先于 aldehyde/radical
    heads = frozenset(e["center_idx"] for e in acyl)
    out = _local_entries(mol, hits)
    result = {**out,
        "radical": [_radical_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("radical", []) if t[0] not in heads],
        "acyl": acyl,
        "aldehyde": [e for e in (_fg_entry(mol.GetAtomWithIdx(t[0])) for t in hits.get("aldehyde", [])) if e["center_idx"] not in heads],
        **oxoacid_lists(mol, hits.get("oxoacid", []))}
    return _drop_claimed_cations(result)


def _drop_claimed_cations(parts: dict) -> dict:
    """已被其它 FG 检测器命中中心的原子不再作阳离子：只留未被 fg 检测器检测到的阳离子。"""
    cations = parts.get("cation") or []
    if not cations:
        return parts
    claimed = {e["center_idx"] for fg, entries in parts.items() if fg != "cation" for e in entries}
    return {**parts, "cation": [e for e in cations if e["center_idx"] not in claimed]}

def _collect_fgs(mol: Mol) -> dict:
    """聚合官能团条目并构建带类型清单（FG 唯一出口）。"""
    from namepredict.layer1.functional_group_inventory import build_inventory

    has_anion = _has_negative_atom(mol)  # 负电荷扫描全流程只算一次
    parts, demoted = _arbitrate_parts(_detect_parts(mol), mol, has_anion)
    double_bonds, triple_bonds = _cc_bond_entries(mol)
    return {"double_bonds": double_bonds, "triple_bonds": triple_bonds,
            "fg_inventory": build_inventory(parts, mol, demoted, has_anion)}

def analyze(mol: Mol) -> dict:
    """分析分子并返回完整的官能团与结构信息 dict（碳信息 + 官能团 + 环事实）。"""
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == C]
    base = {"mol": mol, "carbon_ids": carbons, "n_carbons": len(carbons)}
    return {**base, **_collect_fgs(mol), **_ring_meta(mol)}
