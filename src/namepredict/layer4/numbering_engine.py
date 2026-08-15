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
    """计算各键较小端点位次的排序元组。"""
    if not bonds:
        return None
    mins = []
    for b in bonds:
        if b[0] in cand and b[1] in cand:
            mins.append(min(cand[b[0]], cand[b[1]]))
    return tuple(sorted(mins)) if mins else None


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
    """（全部多重键、双键）端点对。"""
    all_bonds, doubles = [], []
    for key, target in (("double_bond", doubles), ("triple_bond", all_bonds)):
        v = parent.get(key)
        if v:
            pair = (v[0], v[1])
            (target if key == "double_bond" else all_bonds).append(pair)
   
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


def _ring_hetero_start(parent: dict, chain: list[int]) -> int | None:
    """单杂原子环：唯一的杂原子为 locant 1（P-14.4，如吡啶）。"""
    mol = parent.get("mol")
    if mol is None:
        return None
    heteros = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    return heteros[0] if len(heteros) == 1 else None


def _fixed_start(parent: dict) -> int | None:
    """取固定 locant 1 起点原子；无字段则尝试单杂原子环。"""
    for key in _FIXED_START_KEYS:
        v = parent.get(key)
        if v is not None:
            return v
    return _ring_hetero_start(parent, parent.get("chain") or [])


def _fixed_numbering(parent: dict, chain: list[int]) -> list[int] | None:
    """P-14.4(a)：经 L2 numbering_scaffold 的保留骨架固定编号。"""
    plan = None
    return list(plan.atom_order) if plan is not None else None


# ── 入口 ─────────────────────────────────────────────────────────────────

def orient_numbering(parent: dict, substituents: list) -> list[int] | None:
    """返回 P-14.4 定向后的原子顺序；不适用则返回 None。"""
    chain = parent.get("chain") or []
    if not chain:
        return None
    fixed = _fixed_numbering(parent, chain)
    if fixed is not None:
        return fixed
    if _is_ring(parent):
        cands = _ring_cands(chain)
    else:
        cands = _chain_cands(chain)
    start = _fixed_start(parent)
    if start is not None:
        cands = [c for c in cands if c.get(start) == 1]
    if not cands:
        # 固定起点原子在链候选中不可能为 locant 1（如线性链中部的杂环原子）：保留原顺序。
        return chain
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
