"""L4 定向编号引擎：按 P-14.4 规则筛出链/环原子顺序候选。"""
from __future__ import annotations

from namepredict.constants import (
    C, P145_SENIOR, RS_HI, RS_LO, TRADITIONAL_NUMBERING_IDS,
)
from namepredict.tools import memo
from namepredict.tools.re import alpha_order_key
from namepredict.layer1.ring_systems import sssr_rings
from namepredict.layer4.indicated_hydrogen import saturated_ring_atoms
from namepredict.layer4.locant_calc import locant_key


# ── 候选生成 ──────────────────────────

def _numbered(chain: list[int]) -> dict[int, int]:
    """把链原子顺序映射为 {原子: 位次}。"""
    return {a: i + 1 for i, a in enumerate(chain)}


def _to_chain(cand: dict[int, int]) -> list[int]:
    """按位次升序还原原子顺序列表。"""
    return sorted(cand, key=cand.get)


def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    """返回 (基团字母序键, 位次) 排序对，用于字母序平局。"""
    return sorted(
        (alpha_order_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents if s["attach_idx"] in chain
    )


def _ring_cands(chain: list[int]) -> list[dict[int, int]]:
    """生成环的全部旋转/翻转编号候选。"""
    n = len(chain)
    out: list[dict[int, int]] = []
    for i in range(n):
        fwd = chain[i:] + chain[:i]
        out.append(_numbered(fwd))
        out.append(_numbered(list(reversed(fwd))))
    return out


# ── 位次集合键 ─────────────────────────

def _locant_set(cand: dict[int, int], atoms: list[int]) -> tuple[int, ...] | None:
    """计算原子集合在候选编号下的排序位次元组。"""
    locs = sorted(cand[a] for a in atoms if a in cand)
    return tuple(locs) if locs else None


def _edge_locants(pos: dict[int, int], n: int, bonds) -> tuple[int, ...] | None:
    """在已建的编号顺序映射下，算一组键占据的边位次（seam 感知，每条键一个数）。"""
    if not bonds:
        return None
    locs = []
    for a, b in bonds:
        if a not in pos or b not in pos:
            continue
        ia, ib = pos[a], pos[b]
        if ia > ib:
            ia, ib = ib, ia
        if ib - ia == 1:
            locs.append(ia + 1)              # 普通相邻边
        elif ia == 0 and ib == n - 1:
            locs.append(n)                   # seam 闭合边：跨编号首尾，记 n
        else:
            return None                      # 端点不沿编号相邻：判据不适用
    return tuple(sorted(locs)) if len(locs) == len(bonds) else None


def _bond_locant_pairs(cand: dict[int, int], bonds, doubles) -> tuple:
    """一次构建编号顺序，同时给出双键/多重键两组位次（P-14.4(e)）。"""
    if not bonds:
        return None, None
    order = sorted(cand, key=cand.get)  # locant 升序 → 环/链遍历顺序
    pos = {a: i for i, a in enumerate(order)}
    n = len(order)
    return _edge_locants(pos, n, bonds), _edge_locants(pos, n, doubles)


def narrow(cands: list, key_fn, *, reverse: bool = False, skip_none: bool = False) -> list:
    """保留 key_fn 键最小(默认)/最大(reverse)的候选。"""
    if len(cands) <= 1:
        return cands
    keys = [key_fn(c) for c in cands]
    if skip_none and any(k is None for k in keys):
        return cands  # 特征全部缺失 → 规则不适用，候选原样返回
    best = (max if reverse else min)(keys)
    return [c for c, k in zip(cands, keys) if k == best]


def narrow_by_senior(cands: list, key_fn, heteros, by_z, *, skip_none: bool = False) -> list:
    """杂原子集最低位次 → 按 P145_SENIOR 逐元素收窄。"""
    cands = narrow(cands, lambda c: key_fn(c, heteros), skip_none=skip_none)  # (a)
    for z in P145_SENIOR:                                                    # (b)
        atoms = by_z.get(z)
        if atoms:
            cands = narrow(cands, lambda c, at=sorted(atoms): key_fn(c, at), skip_none=skip_none)
    return cands


# ── 桥环/螺环编号共享设施（P-23.3.2 / P-24.2.2 / P-14.4 通用段） ───

def candidates(parent: dict, key: str) -> list:
    """L2 下传的并列候选；无对应节点返回空表。"""
    nodes = parent.get(f"{key}s")
    if nodes:
        return list(nodes)
    node = parent.get(key)
    return [node] if node is not None else []


def hetero_atoms(mol, atoms) -> list[int]:
    """骨架杂原子，按给定原子序列顺序。"""
    if mol is None:
        return []
    return [a for a in atoms if mol.GetAtomWithIdx(a).GetAtomicNum() != C]


def by_z(mol, atoms: list[int]) -> dict[int, list[int]]:
    """原子序数 → 该元素的骨架原子列表。"""
    out: dict[int, list[int]] = {}
    for a in atoms:
        out.setdefault(mol.GetAtomWithIdx(a).GetAtomicNum(), []).append(a)
    return out


def alpha_locants(numbering: dict[int, int], chain: list[int], substituents: list) -> tuple:
    """按前修饰基字母序排列的位次元组（P-14.4(g)）。"""
    pairs = sorted(
        (alpha_order_key(s.get("en") or ""), numbering.get(s["attach_idx"], 0))
        for s in substituents if s.get("attach_idx") in chain)
    return tuple(loc for _, loc in pairs)


def pick_equivalent(nodes: list, feat):
    """并列候选特征全同才取首个，否则视为不可判定。"""
    if not nodes:
        return None
    return nodes[0] if len({feat(nd) for nd in nodes}) == 1 else None


def node_feature_key(mol, chain: list[int], parent: dict, substituents: list, extra):
    """候选渲染特征工厂：描述符 + extra(node) + 杂原子/后缀/取代基位次。"""
    def feat(node) -> tuple:
        het = tuple(sorted((mol.GetAtomWithIdx(a).GetAtomicNum(), node.numbering[a])
                           for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != C)) \
            if mol is not None else ()
        suffix = tuple(sorted(node.numbering[a] for a in _principal_atoms(parent) if a in chain))
        subs = tuple(sorted(node.numbering[s["attach_idx"]] for s in (substituents or [])
                            if s.get("attach_idx") in chain))
        return (node.descriptor, extra(node), het, suffix, subs)
    return feat


# ── P-14.4(j)：CIP 平局破（R/M/r 优先） ───


def assign_cip(mol) -> None:
    """强制重算分子 CIP（隐式 H 手性碳先补显式 H）。"""
    memo.by_mol("cip", _assign_cip_uncached, mol)


def _assign_cip_uncached(mol) -> None:
    """实际执行 CIP 重算（无记忆版本，见 assign_cip）。"""
    from rdkit import Chem
    from rdkit.Chem import ChiralType, rdCIPLabeler

    r = Chem.RWMol(mol)
    for a in list(r.GetAtoms()):
        if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED and a.GetTotalNumHs() == 0 and a.GetDegree() < 4:
            r.AddBond(a.GetIdx(), r.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
    m = r.GetMol()
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)
    rdCIPLabeler.AssignCIPLabels(m)
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode"):
            a.ClearProp("_CIPCode")
        b = m.GetAtomWithIdx(a.GetIdx())
        if b.HasProp("_CIPCode"):
            a.SetProp("_CIPCode", b.GetProp("_CIPCode"))


def _chain_rs_codes(mol, chain: list[int]) -> dict[int, str]:
    """取母体链原子上的 CIP 代码映射（仅 R/S）。"""
    if mol is None or not chain:
        return {}
    assign_cip(mol)
    codes: dict[int, str] = {}
    for idx in chain:
        a = mol.GetAtomWithIdx(int(idx))
        if a.HasProp("_CIPCode") and a.GetProp("_CIPCode") in RS_HI | RS_LO:
            codes[int(idx)] = a.GetProp("_CIPCode")
    return codes


def _rs_locant_key(codes: dict[int, str], chain: list[int], labels: list[str] | None = None) -> tuple:
    """P-14.4(j) 排序键：R/M/r 位次在前、S/P/s 在后，小者优先。"""
    use = labels if labels and len(labels) == len(chain) else None
    hi, lo = [], []
    for i, idx in enumerate(chain):
        code = codes.get(int(idx))
        loc = locant_key(use[i] if use else i + 1)
        if code in RS_HI:
            hi.append(loc)
        elif code in RS_LO:
            lo.append(loc)
    return tuple(hi), tuple(lo)


# ── 从 parent dict 提取 P-14.4 特征 ────

def _principal_atoms(parent: dict) -> list[int]:
    """P-14.4(c)：principal 特征基团的附着原子。"""
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts is not None and facts.attachment_atoms else []


def _unsat_bonds(parent: dict) -> tuple[list, list]:
    """（全部多重键、双键）端点对：单数键亦进 all_bonds。"""
    all_bonds, doubles = [], []
    for key in ("double_bond", "triple_bond"):
        v = parent.get(key)
        if v:
            pair = (v[0], v[1])
            all_bonds.append(pair)
            if key == "double_bond":
                doubles.append(pair)
    for key in ("double_bonds", "triple_bonds"):
        for b in parent.get(key) or []:
            pair = (b[0], b[1])
            all_bonds.append(pair)
            if key == "double_bonds":
                doubles.append(pair)
    return all_bonds, doubles


def _is_ring(parent: dict) -> bool:
    """按 scaffold_id 判断 parent 是否为环系。"""
    return bool(parent.get("scaffold_id"))


def _hydro_indicated_atoms(parent: dict, chain: list[int]) -> list[int]:
    """加氢位与指示氢位（饱和环位）的并集：P-14.4(b)(d)(e)(i) 一起最小化。"""
    if not parent.get("hydro_atoms"):  # 饱和度由 ene 词尾表达者（环烯/环炔）不参与
        return []
    mol = parent.get("mol")
    chain_set = set(chain)
    sats = set(saturated_ring_atoms(mol, chain_set)) & chain_set if mol is not None else set()
    return sorted(sats | set(parent.get("hydro_atoms") or ()))


def _nh_sites(heteros: list[int], mol) -> list[int]:
    """P-14.4(b)/P-31.2.2 指示氢位：环内 NH 及其 N-取代等价位；无 NH 则位次不定。"""
    def _n(a):
        return mol.GetAtomWithIdx(a)
    n_atoms = [a for a in heteros if _n(a).GetAtomicNum() == 7]
    has_h = [a for a in n_atoms if _n(a).GetTotalNumHs() > 0]
    if not has_h:
        return []
    subs = [a for a in n_atoms if _n(a).GetTotalNumHs() == 0 and _n(a).GetDegree() == 3]
    return has_h + subs  # NH 与 N-取代位同为母体指示氢位：位次集并列，交由后续 (c)(f) 裁决

def resolve_numbering(parent: dict, substituents: list, node_key: str, feature_fn, ladder_fn) -> list[int] | None:
    """编号裁决公共外壳：候选 → 收窄阶梯 → 并列等价 → 写回节点并按位次升序返回。"""
    nodes = candidates(parent, node_key)
    chain = list(parent.get("chain") or ())
    if not nodes or not chain:
        return None
    mol = parent.get("mol")
    heteros = hetero_atoms(mol, chain)
    nodes = ladder_fn(nodes, parent, substituents, mol, chain, heteros)
    best = pick_equivalent(nodes, node_feature_key(mol, chain, parent, substituents, feature_fn))
    if best is None:
        return None
    parent[node_key] = best  # L5 依它取描述符/上标，须与选中的编号自洽
    return sorted(best.numbering, key=best.numbering.get)


def _narrow_hetero_ring(cands: list[dict], mol, chain: list[int], float_hetero: bool) -> list[dict]:
    """杂环编号 P-22.2.2.1.3/(b)：位次 1 给最先元素→杂原子集→元素序→指示氢 NH 位次最小化。"""
    heteros = hetero_atoms(mol, chain)
    zmap = by_z(mol, heteros)
    first = next((z for z in P145_SENIOR if z in zmap), None)  # (0) 引用序最先者得位次 '1'
    if first is not None:
        atoms = zmap[first]
        cands = narrow(cands, lambda c: min(c[a] for a in atoms), skip_none=True)
    cands = narrow_by_senior(                                           # (a)(b)
        cands, lambda c, at: _locant_set(c, at), heteros, zmap, skip_none=True)
    if not float_hetero:                                                # (b)
        n_active = _nh_sites(heteros, mol)
        if n_active:
            cands = narrow(cands, lambda c: _locant_set(c, sorted(n_active)), skip_none=True)
    return cands


def _template_matches(mol, sid: str, atoms: frozenset) -> list:
    """模板全覆盖环系的匹配：先直接匹配，失败退到氢化骨架匹配。"""
    from namepredict.layer2.ring_scaffold import _Q, _Q_H, _hydrogenated
    q = _Q.get(sid)
    matches = [m for m in mol.GetSubstructMatches(q, uniquify=False) if set(m) == atoms]
    if matches:
        return matches
    qh = _Q_H.get(sid)
    mol_h = memo.by_mol("hydrogenated", _hydrogenated, mol)
    if qh is None or mol_h is None:
        return []
    return [m for m in mol_h.GetSubstructMatches(qh, uniquify=False) if set(m) == atoms]


def _fixed_numbering(parent: dict, chain: list[int], substituents: list | None = None) -> list[int] | None:
    """P-14.4(a)：fused 环按模板映射固定编号；对称者取位次最小。"""
    sid = parent.get("scaffold_id")
    mol = parent.get("mol")
    if not sid or mol is None:
        return None
    from namepredict.layer2.ring_scaffold import _Q, standard_chain
    if _Q.get(sid) is None:
        return None
    atoms = frozenset(chain)
    chains = []
    for m in _template_matches(mol, sid, atoms):
        std = standard_chain(sid, tuple(m))
        if std is not None and set(std) == atoms:
            chains.append(std)
    if not chains:
        return None
    if len(chains) == 1:
        return chains[0]
    from namepredict.layer2.ring_scaffold import _STANDARD_LABELS
    suffixes = [a for a in _principal_atoms(parent) if a in chain]  # P-14.4(c)：principal 特征基团与自由价附着原子得最低位次。
    if parent.get("radical_c_idx") in chain:  # 自由基主基团：自由价连接点与 principal 同属 (c) 后缀类。
        suffixes.append(parent["radical_c_idx"])
    prefixes = [s["attach_idx"] for s in (substituents or []) if s.get("attach_idx") in chain]
    if not suffixes and not prefixes:
        return chains[0]
    labels = _STANDARD_LABELS.get(sid) or ()
    alpha_subs = [(alpha_order_key(s.get("en") or ""), s["attach_idx"])
                  for s in (substituents or []) if s.get("attach_idx") in chain]

    def _locant_key_of(std: list[int]) -> dict:
        """链上各原子的位次键，按 standard_path 标签而非链位置。"""
        if len(labels) == len(std):
            return {a: locant_key(labels[std.index(a)]) for a in std}
        return {a: (std.index(a) + 1, 0) for a in std}

    def _fixed_key(std: list[int]) -> tuple:
        """候选链 P-14.4 位次键：(c) 后缀、(f) 前缀、(g) 引用序。"""
        loc = _locant_key_of(std)
        return (tuple(sorted(loc[a] for a in suffixes if a in loc)),
                tuple(sorted(loc[a] for a in prefixes if a in loc)),
                tuple(sorted((k, loc[a]) for k, a in alpha_subs if a in loc)))
    return min(chains, key=_fixed_key)


def _fused_numbering(parent: dict, chain: list[int],
                     substituents: list | None = None) -> list[int] | None:
    """P-25.3.3 稠环编号；传统编号例外骨架保持固定编号，其余走优选取向。"""
    sid = parent.get("scaffold_id")
    if sid in TRADITIONAL_NUMBERING_IDS:
        return None
    if parent.get("bridged_node") is not None:  # 桥环走 P-23：chain_fused 判据对桥环恒真
        return None
    if parent.get("spiro_node") is not None:  # 螺环走 P-24：chain_fused 判据对螺环同样恒真
        return None
    mol = parent.get("mol")
    from namepredict.layer2.ring_scaffold import get_spec
    spec = get_spec(parent.get("scaffold_id") or "")
    registered_fused = bool(spec and spec.n_rings >= 2)  # 稠环保留模板：编号随骨架拓扑，与芳香性无关。
    chain_arom = bool(mol is not None and chain) and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in chain)
    chain_fused = bool(mol is not None and chain) and sum(  # 未注册稠环（fused_hetero）：放行走 P-25.3.3 外周编号。
        1 for r in sssr_rings(mol) if set(r) <= set(chain)) >= 2
    if mol is None or not chain or not (chain_arom or registered_fused or chain_fused):
        return None
    from namepredict.layer1.ring_systems import build_ring_systems
    from namepredict.layer4.fused_orientation import preferred_orientations
    from namepredict.layer4.fused_numbering import number_fused_system
    systems = [s for s in build_ring_systems(mol) if (s.get("atom_ids") or []) == sorted(set(chain))]
    if not systems:
        return None
    system = systems[0]
    sssr = list(system.get("sssr_indices") or [])
    if len(sssr) < 2:
        return None  # 单环走 P-14.4 通用枚举
    atom_rings = list(sssr_rings(mol))  # 仅取稠合系统自身环，避免无关环致 layout 失败。
    rings = [atom_rings[i] for i in sssr]
    idx_map = {orig: new for new, orig in enumerate(sssr)}
    fusion_edges = [(idx_map[i], idx_map[j], sh) for i, j, sh in (system.get("fusion_edges") or [])
                    if i in idx_map and j in idx_map]
    orients = preferred_orientations(mol, rings, fusion_edges)  # print("riings",rings) / print(orients)
    if not orients:
        return None
    chain_set = set(chain)  # 环外附着原子按优先级分层，逐层收窄镜像取向。
    layers: list[list[int]] = []
    radical = parent.get("radical_c_idx")
    if radical in chain_set:
        layers.append([radical])  # P-29: 取代基游离价连接点优先得最低位次
    principal_atoms = [a for a in _principal_atoms(parent) if a in chain_set]
    if principal_atoms:
        layers.append(principal_atoms)  # P-14.4(c): principal 特征基团优先于取代基
    from namepredict.layer4.fused_numbering import INDICATED_H
    layers.append(INDICATED_H)  # P-25.3.3.1.2(f): 指示氢位次插在 (c) 之后。
    hydro_atoms = sorted(a for a in (parent.get("hydro_atoms") or ()) if a in chain_set)
    if hydro_atoms:
        layers.append(hydro_atoms)  # P-14.4(e)(i): 加氢位次先于 (f) 可分离前缀；P-31.2.2 指示氢仍优先
    sub_atoms = sorted(s["attach_idx"] for s in (substituents or [])
                       if s.get("attach_idx") in chain_set)
    if sub_atoms:
        layers.append(sub_atoms)  # P-14.4(f): 取代基位次集合最小化
    alpha_subs = [(alpha_order_key(s.get("en") or ""), s["attach_idx"])
                  for s in (substituents or []) if s.get("attach_idx") in chain_set]
    result = number_fused_system(mol, rings, [o.coord_dict() for o in orients], layers, alpha_subs)
    if result is None:
        return None
    fused_chain, labels = result
    parent["numbering_scaffold"] = {
        "scaffold_id": "fused", "labels": tuple(labels),
    }
    return fused_chain #外环原子顺序


# ── 入口 ────────────────────────────

def orient_numbering(parent: dict, substituents: list, *, float_hetero: bool = False) -> list[int] | None:
    """返回 P-14.4 定向后的原子顺序，不适用返回 None。"""
    chain = parent.get("chain") or []
    if not chain:
        return None
    if parent.get("fbs_nodes") is not None:  # P-24.5+ 组分式螺环：逐组分裁决
        from namepredict.layer4.spiro_numbering import fbs_numbering
        return fbs_numbering(parent, substituents)
    if parent.get("spiro_nodes") is not None:  # P-24 螺环：编号只能由 L2 候选裁决
        from namepredict.layer4.spiro_numbering import spiro_numbering
        return spiro_numbering(parent, substituents)  # 无论成败都不再下落（_fused 会误吞螺环）
    if parent.get("bridged_nodes") is not None:  # P-23 桥环：编号只能由 L2 候选裁决
        from namepredict.layer4.bridged_numbering import bridged_numbering
        return bridged_numbering(parent, substituents)  # 无论成败都不再下落（_fused 会误吞桥环）
    fixed = _fixed_numbering(parent, chain, substituents)
    if fixed is not None:
        return fixed
    fused = _fused_numbering(parent, chain, substituents)
    if fused is not None:
        return fused
    mol = parent.get("mol")
    if parent.get("bridged_node") is not None:  # 桥环兜底：chain 非环序，旋转/翻转枚举无意义
        return None
    is_ring = _is_ring(parent)
    if mol is not None and is_ring and any(
            mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in chain):
        cands = _narrow_hetero_ring(_ring_cands(chain), mol, chain, float_hetero)  # 杂环：P-22.2.2.1.3 元素序窄化先于 principal。
    else:
        cands = _ring_cands(chain) if is_ring else [  # 碳环/链：P-14.4(a) 固定 locant 1 锚定后退化
            _numbered(chain), _numbered(list(reversed(chain)))]  # 链：正反两个方向的编号候选
    principal = _principal_atoms(parent)
    if principal:
        cands = narrow(cands, lambda c: _locant_set(c, principal), skip_none=True)
    hydro = _hydro_indicated_atoms(parent, chain)
    if hydro:  # P-14.4(b)(d)(e)(i)：加氢/指示氢位次先于多重键与取代基
        cands = narrow(cands, lambda c: _locant_set(c, hydro), skip_none=True)
    bonds, doubles = _unsat_bonds(parent)
    if bonds:
        cands = narrow(cands, lambda c: _bond_locant_pairs(c, bonds, doubles), skip_none=True)
    subs = [s["attach_idx"] for s in substituents if s.get("attach_idx") in chain]

    if subs:
        cands = narrow(cands, lambda c: _locant_set(c, subs), skip_none=True)
    if len(cands) > 1 and substituents:  # P-14.4(f) 平局：最低位次给字母序最前的取代基。
        cands = narrow(cands, lambda c: _stem_loc_pairs(_to_chain(c), substituents), skip_none=True)
    if len(cands) > 1:
        codes = _chain_rs_codes(mol, chain)  # P-14.4(j) 立体平局：按 CIP 描述符定方向，低位次给 R/M/r。
        if codes:
            cands = narrow(cands, lambda c: _rs_locant_key(codes, _to_chain(c)), skip_none=True)
    return _to_chain(cands[0])


def _component_labels(parent: dict, chain: list[int]) -> list[str]:
    """补并行 labels：优先标准 locant 标签，否则纯数字。"""
    scaffold = parent.get("numbering_scaffold")
    labels = scaffold.get("labels") if scaffold else None
    if labels and len(labels) == len(chain):
        return list(labels)
    from namepredict.layer2.ring_scaffold import _STANDARD_LABELS
    labels = _STANDARD_LABELS.get(parent.get("scaffold_id") or "")
    if labels and len(labels) == len(chain):
        return list(labels)
    return [str(i + 1) for i in range(len(chain))]


def fused_component_numbering(mol, scaffold_id, sub_rings, shared=None, sub_edges=None):
    """稠合组分自身编号 P-25.4/P-25.3.3，稠合点作取代基最小化。"""
    if not sub_rings:
        return None, None
    subs = [{"attach_idx": a} for a in (shared or ())] if shared else []
    chain0 = sorted(set().union(*sub_rings))
    from namepredict.layer2.ring_scaffold import _STANDARD_ORDERS
    if len(sub_rings) == 1 or (scaffold_id and scaffold_id in _STANDARD_ORDERS):
        if not scaffold_id:
            ring = list(sub_rings[0])
            return ring, [str(i + 1) for i in range(len(ring))]
        parent = {"mol": mol, "scaffold_id": scaffold_id, "chain": chain0}
        res = orient_numbering(parent, subs, float_hetero=bool(shared))  # float_hetero: 对称杂环留两镜像，locant 1 由稠合原子定。
        if not res:
            return None, None
        return res, _component_labels(parent, res)
    from namepredict.layer4.fused_orientation import preferred_orientations  # 多环无固定编号: P-25.3.3 通用编号(fused_numbering)
    from namepredict.layer4.fused_numbering import number_fused_system
    orients = preferred_orientations(mol, sub_rings, sub_edges or [])
    if not orients:
        return None, None
    result = number_fused_system(mol, sub_rings, [o.coord_dict() for o in orients],
                                 [sorted(shared)] if shared else None)  # 稠合点作 sub_layers 逐层最小化位次(P-25.3.1.3)。
    if result is None:  # print(result)
        return None, None
    return result[0], result[1]
