"""L4 定向编号引擎：按 P-14.4 规则筛出链/环原子顺序候选。"""
from __future__ import annotations



# ── 候选生成 ──────────────────────────────────────────────────

def _numbered(chain: list[int]) -> dict[int, int]:
    """把链原子顺序映射为 {原子: 位次}。"""
    return {a: i + 1 for i, a in enumerate(chain)}


def _to_chain(cand: dict[int, int]) -> list[int]:
    """按位次升序还原原子顺序列表。"""
    return sorted(cand, key=cand.get)


def _chain_cands(chain: list[int]) -> list[dict[int, int]]:
    """生成线性链正反两个方向的编号候选。"""
    return [_numbered(chain), _numbered(list(reversed(chain)))]


def _ring_cands(chain: list[int]) -> list[dict[int, int]]:
    """生成环的全部旋转/翻转编号候选。"""
    n = len(chain)
    out: list[dict[int, int]] = []
    for i in range(n):
        fwd = chain[i:] + chain[:i]
        out.append(_numbered(fwd))
        out.append(_numbered(list(reversed(fwd))))
    return out


# ── 位次集合键 ───────────────────────────────────────────────────────

def _locant_set(cand: dict[int, int], atoms: list[int]) -> tuple[int, ...] | None:
    """计算原子集合在候选编号下的排序位次元组。"""
    locs = sorted(cand[a] for a in atoms if a in cand)
    return tuple(locs) if locs else None


def _bond_locants(cand: dict[int, int], bonds) -> tuple[int, ...] | None:
    """计算各多重键两端点位次的全排序元组（P-14.4(e) 多重键低位比较）。

    不能只取每键较小端点：环候选可绕行让闭合键一端成为 locant 1，把真实双键对
    （如环己-1,3-二烯的 (1,2),(3,4)）误判成更"小"的 (1,2)。两端点全 locant 排序
    迫使编号把整条不饱和键落到最低连续位置（1-2 与 3-4 → 1,3-diene）。
    """
    if not bonds:
        return None
    locs = []
    for b in bonds:
        if b[0] in cand and b[1] in cand:
            locs.append(cand[b[0]])
            locs.append(cand[b[1]])
    return tuple(sorted(locs)) if locs else None


def _narrow(cands: list[dict], key_fn) -> list[dict]:
    """保留位次集合键最小的候选；得到一个即提前停止。"""
    if len(cands) <= 1:
        return cands
    keys = [key_fn(c) for c in cands]
    if any(k is None for k in keys):
        return cands  # 特征全部缺失 → 规则不适用
    best = min(keys)
    return [c for c, k in zip(cands, keys) if k == best]


# ── 从 parent dict 提取 P-14.4 特征 ────────────────────────

def _principal_atoms(parent: dict) -> list[int]:
    """P-14.4(c)：principal characteristic group 的附着原子。"""
    atoms = parent.get("principal_attachment_atoms")
    if atoms:
        return list(atoms)
    facts = parent.get("principal_expression_facts")
    if facts is not None:
        attach = getattr(facts, "attachment_atoms", None)
        if attach:
            return sorted(attach)
    return []


def _unsat_bonds(parent: dict) -> tuple[list, list]:
    """（全部多重键、双键）端点对：单数 double_bond 亦进 all_bonds（混合烯炔整体最小化）。"""
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


# P-14.4(a)：parent dict 中标定必须为 locant 1 的原子的字段（杂环起点、环外羰基连接、自由基中心）。
_FIXED_START_KEYS = (
    "ring_attach_idx", "n_idx", "nh_idx", "hetero_idx", "radical_c_idx",
)


def _ring_hetero_start(parent: dict, chain: list[int], float_hetero: bool = False) -> int | None:
    """杂原子环：最优先杂原子为 locant 1（P-14.4，吡啶/嘧啶等）。

    多杂环（嘧啶双 N、咪唑等）同样固定一个杂原子为起点，其余杂原子在
    P-14.4 枚举中自然得低位（1,3 / 1,2）。1,3-二唑（咪唑/吡唑）N1 优先取
    带取代基（价 3）的 N、其次带 H 的 N（P-58.2.1）；其余退化为 Z 最小。

    ``float_hetero=True``（仅稠合组分编号 P-25.3.2.5 使用）：对称等价杂环
    （嘧啶双 N 等）无唯一 NH/Nsub 起点时允许 locant 1 在两个等价杂原子间浮动，
    交由后续取代基（稠合原子）位次最小化决定——避免按原子序号随机钉死起点
    得到与规范稠合描述符（d 侧等）不一致的镜像。
    """
    mol = parent.get("mol")
    if mol is None:
        return None
    heteros = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    if not heteros:
        return None
    n_sub = [a for a in heteros if mol.GetAtomWithIdx(a).GetAtomicNum() == 7
             and mol.GetAtomWithIdx(a).GetDegree() == 3]
    if len(n_sub) == 1:
        return n_sub[0]
    n_nh = [a for a in heteros if mol.GetAtomWithIdx(a).GetAtomicNum() == 7
            and mol.GetAtomWithIdx(a).GetTotalNumHs() > 0]
    if len(n_nh) == 1:
        return n_nh[0]
    if float_hetero and len(heteros) > 1:
        return None
    return min(heteros, key=lambda a: (mol.GetAtomWithIdx(a).GetAtomicNum(), a))


def _fixed_start(parent: dict, float_hetero: bool = False) -> int | None:
    """取固定 locant 1 起点原子：杂原子环优先杂原子（P-14.4，吡啶甲酸 N=1
    而非羧酸锚点），否则 FG 锚点/自由基字段，最后退化处理。"""
    hetero = _ring_hetero_start(parent, parent.get("chain") or [], float_hetero)
    if hetero is not None:
        return hetero
    for key in _FIXED_START_KEYS:
        v = parent.get(key)
        if v is not None:
            return v
    return None


def _fixed_numbering(parent: dict, chain: list[int], substituents: list | None = None) -> list[int] | None:
    """P-14.4(a)：经 L2 保留骨架固定编号（standard_path + 模板匹配映射）。

    fused 芳香环（喹啉/吲哚等）的 IUPAC 编号固定（P-25.4）：起点/方向不随取代基
    变化，P-14.4 通用环枚举会算错（如喹啉 10-氯 vs 2-氯）。standard_path 定义
    模板原子按标准 locant 的顺序；对称 scaffold（phenanthrene）子图匹配方向随
    分子原子编号漂移，用模板全部自同构 match 生成等价链，按取代基位次最小化。
    """
    sid = parent.get("scaffold_id")
    mol = parent.get("mol")
    if not sid or mol is None:
        return None
    from namepredict.layer2.ring_scaffold import _Q, standard_chain
    q = _Q.get(sid)
    if q is None:
        return None
    atoms = frozenset(chain)
    chains = []
    for m in mol.GetSubstructMatches(q, uniquify=False):
        if set(m) == atoms:
            std = standard_chain(sid, tuple(m))
            if std is not None and set(std) == atoms:
                chains.append(std)
    if not chains:
        return None
    if len(chains) == 1:
        return chains[0]
    subs = [s["attach_idx"] for s in (substituents or []) if s.get("attach_idx") in chain]
    if not subs:
        return chains[0]
    return min(chains, key=lambda std: tuple(sorted(std.index(a) + 1 for a in subs)))


def _fused_numbering(parent: dict, chain: list[int]) -> list[int] | None:
    """P-25.3.3 稠环编号：仅未注册(非保留模板)芳香稠环系统。

    护栏: 仅 scaffold_id 为 None/"carbocycle" 且全芳香多环走优选取向+外周编号;
    registered 模板(含对称 naphthalene/anthracene)保持固定编号/P-14.4 不被接管。
    """
    sid = parent.get("scaffold_id")
    if sid not in (None, "carbocycle", "fused_hetero"):
        return None
    mol = parent.get("mol")
    if mol is None or not chain or not all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in chain):
        return None
    from namepredict.layer1.ring_systems import build_ring_systems
    from namepredict.layer4.fused_orientation import preferred_orientations
    from namepredict.layer4.fused_numbering import number_fused_system
    systems = [s for s in build_ring_systems(mol) if (s.get("atom_ids") or []) == sorted(set(chain))]
    if not systems:
        return None
    system = systems[0]
    if len(system.get("sssr_indices") or []) < 2:
        return None  # 单环走 P-14.4 通用枚举
    rings = list(mol.GetRingInfo().AtomRings())
    print("riings",rings)
    orients = preferred_orientations(mol, rings, system["fusion_edges"])
    # print(orients)
    if not orients:
        return None
    result = number_fused_system(mol, rings, [o.coord_dict() for o in orients])
    if result is None:
        return None
    fused_chain, labels = result
    parent["numbering_scaffold"] = {
        "scaffold_id": "fused", "labels": tuple(labels), "relative_stereo": None,
    }
    return fused_chain #外环原子顺序


# ── 入口 ─────────────────────────────────────────────────────────────────

def orient_numbering(parent: dict, substituents: list, *, float_hetero: bool = False) -> list[int] | None:
    """返回 P-14.4 定向后的原子顺序；不适用则返回 None。

    ``float_hetero`` 见 ``_ring_hetero_start``（仅稠合组分编号路径开启）。
    """
    chain = parent.get("chain") or []
    if not chain:
        return None
    fixed = _fixed_numbering(parent, chain, substituents)
    if fixed is not None:
        return fixed
    fused = _fused_numbering(parent, chain)
    if fused is not None:
        return fused
    if _is_ring(parent):
        cands = _ring_cands(chain)
    else:
        cands = _chain_cands(chain)
    start = _fixed_start(parent, float_hetero)
    if start is not None:
        cands = [c for c in cands if c.get(start) == 1]
    if not cands:
        # 固定起点原子在链候选中不可能为 locant 1（如线性链中部的杂环原子）：保留原顺序。
        return chain
    mol = parent.get("mol")
    if mol is not None and _is_ring(parent):
        # P-14.4：环系编号使杂原子得最低 locant（唑类 N 必须 1,3/1,2）；
        # 先于 principal 最小化，避免 principal 翻转破坏唑环固定编号。
        heteros = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
        if heteros:
            cands = _narrow(cands, lambda c: _locant_set(c, heteros))
    principal = _principal_atoms(parent)
    if principal:
        cands = _narrow(cands, lambda c: _locant_set(c, principal))
    bonds, doubles = _unsat_bonds(parent)
    if bonds:
        cands = _narrow(cands, lambda c: (_bond_locants(c, bonds), _bond_locants(c, doubles)))
    subs = [s["attach_idx"] for s in substituents if s.get("attach_idx") in chain]
    
    if subs:
        cands = _narrow(cands, lambda c: _locant_set(c, subs))
    if len(cands) > 1 and substituents:
        # P-14.4(f) 平局：最低位次集合已相同 → 把最低位次给字母序最前的取代基（stem-alpha 对，P-14.5）。
        from namepredict.layer4._chain_orient import _stem_loc_pairs
        cands = _narrow(cands, lambda c: _stem_loc_pairs(_to_chain(c), substituents))
    return _to_chain(cands[0])


def _component_labels(parent: dict, chain: list[int]) -> list[str]:
    """orient_numbering 只返回链, 补并行 labels: 优先标准 locant 标签, 否则纯数字。"""
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
    """稠合组分的自身编号(P-25.4/P-25.3.3), 稠合点 shared 作取代基最小化位次。
    返回 (chain, labels) 或 (None, None); 供 L5 fused_namer 组装稠合名。"""
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
        # float_hetero: 对称等价杂环(嘧啶双 N 等)的 locant 1 交给稠合原子位次
        # 最小化决定, 使碱环取向与规范稠合字母(d 侧)一致(见 _ring_hetero_start)。
        res = orient_numbering(parent, subs, float_hetero=bool(shared))
        if not res:
            return None, None
        return res, _component_labels(parent, res)
    # 多环无固定编号: P-25.3.3 通用编号(fused_numbering)。
    from namepredict.layer4.fused_orientation import preferred_orientations
    from namepredict.layer4.fused_numbering import number_fused_system
    orients = preferred_orientations(mol, sub_rings, sub_edges or [])
    if not orients:
        return None, None
    result = number_fused_system(mol, sub_rings, [o.coord_dict() for o in orients])
    # print(result)
    if result is None:
        return None, None
    return result[0], result[1]
