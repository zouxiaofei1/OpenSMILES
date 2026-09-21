"""P-24 螺环拆解：环组分划分 + von Baeyer 螺描述符 + 组分式螺环候选。

全为单环组分时走 P-24.2 螺描述符，含多环组分时走 P-24.5~24.7 组分式命名。
L2 只枚举候选，编号裁决交 L4；L5 不得 import L2，node 按鸭子类型读。
"""
from __future__ import annotations

import re
from dataclasses import dataclass, replace

from namepredict.constants import C
from namepredict.layer1.ring_systems import sssr_rings
from namepredict.layer2.ring_scaffold import match_retained

_BRACKET = re.compile(r"\[([^\]]*)\]")
_MAX_WALKS = 64  # 段回路候选上限
SPIRO_SCAFFOLDS = ("mono_spiro", "fused_bridged_spiro")


@dataclass(frozen=True)
class SpiroComponent:
    """一个环组分：螺环系内与他环仅共用螺原子的单环（P-24.2）。"""

    index: int
    atom_ids: tuple[int, ...]
    ring_index: int
    spiro_atoms: tuple[int, ...]


@dataclass(frozen=True)
class _Seg:
    """螺描述符中的一段：环内沿 a→b 的非螺原子序列（端环则 a=b）。"""

    ring_index: int
    a: int
    b: int
    atoms: tuple[int, ...]
    terminal: bool


@dataclass(frozen=True)
class SpiroNode:
    """P-24 螺环母体：环组分 + 螺描述符 + 原子→位次映射。"""

    scaffold_id: str
    atom_ids: tuple[int, ...]
    free_spiro_atoms: tuple[int, ...]
    components: tuple[SpiroComponent, ...]
    descriptor: tuple[int, ...]
    descriptor_superscripts: tuple[int, ...]
    numbering: dict[int, int]
    ring: str
    n_rings: int

    def scaffold_identity(self):
        """转为 L2 scaffold 身份（id 与命名类同为族名）。"""
        from namepredict.layer2.ring_scaffold import ScaffoldIdentity
        return ScaffoldIdentity(self.scaffold_id, self.scaffold_id, self.n_rings, self.ring)


def _cycle_from(ring, s: int) -> list[int]:
    """环序旋转为从 s 出发的单向序列，s 居首。"""
    i = list(ring).index(s)
    return list(ring[i:]) + list(ring[:i])


def _flip(seg: _Seg) -> _Seg:
    """掉转段的方向，使 atoms 沿 b→a 序。"""
    return replace(seg, a=seg.b, b=seg.a, atoms=tuple(reversed(seg.atoms)))


def _segments(rings, indices, spiros: frozenset[int]) -> list[_Seg] | None:
    """按环分解出全部段；环组分非单环（含稠合环）返回 None。"""
    out: list[_Seg] = []
    for i in indices:
        ring = list(rings[i])
        sp = [a for a in ring if a in spiros]
        if not sp:
            return None
        if len(sp) == 1:  # 端环：整环绕螺原子一周为一段
            rot = _cycle_from(ring, sp[0])
            out.append(_Seg(i, sp[0], sp[0], tuple(rot[1:]), True))
            continue
        rot = _cycle_from(ring, sp[0])  # 中心环：环序上相邻螺原子间各成一段
        pos = [j for j, a in enumerate(rot) if a in spiros]
        for k, j in enumerate(pos):
            j2 = pos[(k + 1) % len(pos)]
            atoms = tuple(rot[j + 1:j2]) if k + 1 < len(pos) else tuple(rot[j + 1:] + rot[:j2])
            out.append(_Seg(i, rot[j], rot[j2], atoms, False))
    return out


def _step(segs, cur: int, remaining: tuple[int, ...], order: tuple, out: list, cap: int) -> None:
    """深度优先枚举段回路：端环段优先、短段优先（P-24.2.2 最短路径）。"""
    if len(out) >= cap:
        return
    if not remaining:
        out.append(order)
        return
    picks = []
    for k in remaining:
        for rev in (False, True):
            f = _flip(segs[k]) if rev else segs[k]
            if f.a == cur:
                picks.append((0 if f.terminal else 1, len(f.atoms), k, rev))
    picks.sort()
    for _, _, k, rev in picks:
        f = _flip(segs[k]) if rev else segs[k]
        _step(segs, f.b, tuple(x for x in remaining if x != k), order + ((k, rev),), out, cap)


def _enumerate_walks(segs: list[_Seg], cap: int = _MAX_WALKS) -> list[tuple]:
    """从各端环出发枚举段回路（P-24.2.2 起自端环，绕行回首个螺原子）。"""
    out: list[tuple] = []
    starts = [k for k, s in enumerate(segs) if s.terminal] or list(range(len(segs)))
    for k in starts:
        for rev in (False, True):
            if len(out) >= cap:
                return out
            f = _flip(segs[k]) if rev else segs[k]
            _step(segs, f.b, tuple(x for x in range(len(segs)) if x != k), ((k, rev),), out, cap)
    return out


def _number_walk(segs: list[_Seg], walk: tuple) -> tuple[dict, tuple, tuple]:
    """按段引用顺序编号：段内原子先编，螺原子首次遇到时得位次（P-24.2.2）。"""
    num: dict[int, int] = {}
    desc: list[int] = []
    sups: list[int] = []
    nxt = 1
    for k, rev in walk:
        seg = _flip(segs[k]) if rev else segs[k]
        desc.append(len(seg.atoms))
        for a in seg.atoms:
            num[a] = nxt
            nxt += 1
        if seg.b in num:
            sups.append(num[seg.b])  # 重访螺原子：已定位次作上标
        else:
            num[seg.b] = nxt
            nxt += 1
            sups.append(0)
    return num, tuple(desc), tuple(sups)


def _ring_kind(mol, atom_ids) -> str:
    """环组分的元素类型：全碳为 carbo，否则 hetero。"""
    return "carbo" if all(mol.GetAtomWithIdx(a).GetAtomicNum() == C
                          for a in atom_ids) else "hetero"


def spiro_scaffold_identity(info: dict, system: dict):
    """螺环系的 scaffold 身份：全单环组分为 mono_spiro，否则 fused_bridged_spiro。"""
    from namepredict.layer2.ring_scaffold import ScaffoldIdentity
    mol = info["mol"]
    rings = list(sssr_rings(mol))
    indices = sorted(system.get("sssr_indices") or ())
    spiros = frozenset(system.get("free_spiro_atoms") or ())
    single = _segments(rings, indices, spiros) is not None  # 每个环都是组分单环
    sid = "mono_spiro" if single else "fused_bridged_spiro"
    atom_ids = tuple(system.get("atom_ids") or ())
    return ScaffoldIdentity(sid, sid, len(indices), _ring_kind(mol, atom_ids))


def decompose_spiro_system(info: dict, system: dict) -> list:
    """公共入口：环系拆解为并列螺环候选（单环组分→SpiroNode，多环→FbsNode）。"""
    free = tuple(system.get("free_spiro_atoms") or ())
    if not free:
        return []
    mol = info["mol"]
    rings = list(sssr_rings(mol))
    indices = sorted(system.get("sssr_indices") or ())
    atom_ids = tuple(system.get("atom_ids") or ())
    kind = _ring_kind(mol, atom_ids)
    spiros = frozenset(free)
    segs = _segments(rings, indices, spiros)
    if segs is None:  # 含多环组分：转 P-24.5~24.7 组分式命名
        return decompose_fbs_system(info, system)
    comps = tuple(SpiroComponent(k, tuple(sorted(rings[i])), i,
                                 tuple(a for a in rings[i] if a in spiros))
                  for k, i in enumerate(indices))
    walks = [(_number_walk(segs, w)) for w in _enumerate_walks(segs)]
    if len(free) == 1:  # P-24.2.1 单螺描述符无上标（上标是 P-24.2.2 多螺专用）
        walks = [(num, desc, tuple(0 for _ in sups)) for num, desc, sups in walks]
    out = [SpiroNode("mono_spiro", atom_ids, free, comps, desc, sups, num, kind, len(indices))
           for num, desc, sups in walks]
    return [nd for nd in out if _numbers_all(nd, atom_ids)]


def _numbers_all(node: SpiroNode, atom_ids) -> bool:
    """编号是否覆盖骨架全部原子且无重号。"""
    return len(node.numbering) == len(atom_ids) and set(node.numbering.values()) == set(
        range(1, len(atom_ids) + 1))


# ── P-24.5~24.7 组分式螺环（含多环组分） ────────────────

@dataclass(frozen=True)
class FbsNumbering:
    """单个组分的一个编号候选（本位次，不含撇号）。"""

    chain: tuple[int, ...]
    labels: tuple[str, ...]

    @property
    def locants(self) -> dict[int, str]:
        """原子 → 本位次标签。"""
        return dict(zip(self.chain, self.labels))


@dataclass(frozen=True)
class FbsComponent:
    """一个环组分：单环/稠环/桥环 + 自身基名与编号候选。"""

    index: int
    kind: str  # mono_ring / fused_ring / bridged_ring
    atom_ids: tuple[int, ...]
    ring_indices: tuple[int, ...]
    spiro_atoms: tuple[int, ...]
    base_en: str | None
    base_zh: str | None
    bare_en: str | None  # 裸词干（可接 a-1,3-diene）；不可接时为 None
    bare_zh: str | None
    sid: str | None
    node: object | None  # BridgedNode / FusedNode / None
    numberings: tuple[FbsNumbering, ...]
    order_extra: tuple = ()  # P-24.5.3 平局键（斜体稠合字母 / von Baeyer 描述符）


@dataclass(frozen=True)
class FbsNode:
    """P-24.5~24.7 组分式螺环母体：引用序已定的组分 + 螺接连边。"""

    atom_ids: tuple[int, ...]
    free_spiro_atoms: tuple[int, ...]
    components: tuple[FbsComponent, ...]
    links: tuple[tuple[int, int, int], ...]  # (螺原子, 父槽位, 本槽位)，按引用序
    cites: tuple[tuple[int, ...], ...]  # 每个引用项覆盖的槽位（>1 即 bis/tris 组）
    ring: str
    n_rings: int

    @property
    def scaffold_id(self) -> str:
        """L2 scaffold 身份 id。"""
        return "fused_bridged_spiro"

    def scaffold_identity(self):
        """转为 L2 scaffold 身份（id 与命名类同为 fused_bridged_spiro）。"""
        from namepredict.layer2.ring_scaffold import ScaffoldIdentity
        return ScaffoldIdentity("fused_bridged_spiro", "fused_bridged_spiro",
                                self.n_rings, self.ring)


# ── 环组分划分 ──────────────────────────────────────────

def component_indices(system: dict) -> list[list[int]]:
    """环系按自由螺原子切成环组分（非自由螺原子不切断），返回环下标组。"""
    idx = list(system.get("sssr_indices") or ())
    free = set(system.get("free_spiro_atoms") or ())
    parent = {r: r for r in idx}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i, j, _ in system.get("fusion_edges") or ():
        union(i, j)
    for i, j, s in system.get("spiro_edges") or ():
        if s not in free:  # 非自由螺原子不切断，两侧仍属同一多环组分
            union(i, j)
    out: dict[int, list[int]] = {}
    for r in idx:
        out.setdefault(find(r), []).append(r)
    return [v for _, v in sorted(out.items(), key=lambda kv: min(kv[1]))]


def _component_system(system: dict, comp: list[int], atoms: tuple[int, ...]) -> dict:
    """构造只含本组分的环系 dict（atom_ids 须为组分自身，否则环数会算错）。"""
    cs = set(comp)
    local = {r: k for k, r in enumerate(sorted(comp))}
    return {
        "atom_ids": list(atoms),
        "sssr_indices": sorted(comp),
        "fusion_edges": [(local[i], local[j], sh) for i, j, sh in
                         system.get("fusion_edges") or () if i in cs and j in cs],
        "spiro_edges": [], "free_spiro_atoms": [],
    }


# ── 组分编号候选 ────────────────────────────────────────

def _from_numbering(num: dict[int, int]) -> FbsNumbering:
    """整数编号映射 → 编号候选（位次写成字符串）。"""
    chain = tuple(sorted(num, key=num.get))
    return FbsNumbering(chain, tuple(str(num[a]) for a in chain))


def _mono_numberings(mol, ring) -> list[FbsNumbering]:
    """单环组分：全旋转/翻转候选，按杂原子优先性收窄（P-22.2.2.1.3）。"""
    from namepredict.layer4.numbering_engine import _narrow_hetero_ring, _ring_cands
    chain = list(ring)
    cands = _ring_cands(chain)
    if any(mol.GetAtomWithIdx(a).GetAtomicNum() != C for a in chain):
        cands = _narrow_hetero_ring(cands, mol, chain, False)
    return [_from_numbering(c) for c in cands]


def _fused_component_numbering(mol, sid, sub_rings, spiros, sub_edges):
    """委托 L4 给多环组分编号（稠合点作 sub_layer 逐层最小化，P-25.3.1.3）。"""
    from namepredict.layer4.numbering_engine import fused_component_numbering
    return fused_component_numbering(mol, sid, sub_rings, list(spiros) or None, sub_edges)


def _retained_numberings(mol, sid: str, matches: list, sub_rings: list,
                         spiros, sub_edges: list) -> list[FbsNumbering]:
    """保留名组分：固定编号视图优先，否则按组分自身取向编号（P-25.3.3）。

    对称保留母体的模板自同构各给一个候选（如 2-benzofuran 的 1/3 位互换），
    哪个取向入选由 L4 按 P-24.5.2 的螺位次最小定。
    """
    from namepredict.layer2.ring_scaffold import _STANDARD_LABELS, standard_chain
    labels = _STANDARD_LABELS.get(sid) or ()
    fixed = [FbsNumbering(tuple(c), tuple(labels)) for c in
             (standard_chain(sid, m) for m in matches) if c and len(labels) == len(c)]
    if fixed:
        return fixed
    if len(sub_rings) == 1:
        return _mono_numberings(mol, sub_rings[0])
    chain, labs = _fused_component_numbering(mol, sid, sub_rings, spiros, sub_edges)
    return [] if not chain or not labs else [FbsNumbering(tuple(chain), tuple(labs))]


def _retained_matches(info: dict, sid: str, atoms, primary: tuple[int, ...]) -> list:
    """保留模板在组分原子集上的全部映射（含氢化骨架回退，对称母体给出取向候选）。"""
    from namepredict.layer2.ring_scaffold import _Q, _Q_H, _hydrogenated
    from namepredict.tools import memo
    want = set(atoms)
    hits = [tuple(m) for m in info["mol"].GetSubstructMatches(_Q[sid], uniquify=False)
            if set(m) == want]
    if hits:
        return hits
    mol_h = memo.by_mol("hydrogenated", _hydrogenated, info["mol"])
    if mol_h is None:
        return [primary]
    hits = [tuple(m) for m in mol_h.GetSubstructMatches(_Q_H[sid], uniquify=False)
            if set(m) == want]
    return hits or [primary]


def _flat_numberings(info, system, comp, atoms, spiros, sub_rings, sub_edges):
    """多环非保留组分：先后试 von Baeyer（P-23）与稠合树（P-25）。"""
    from namepredict.layer2.bridged_system import decompose_bridged_system
    from namepredict.layer2.fused_system import decompose_fused_system
    sub = _component_system(system, comp, atoms)
    nodes = decompose_bridged_system(info, sub)
    if nodes:
        from namepredict.layer5.bridged_namer import bridged_body_names
        body = bridged_body_names(nodes[0])
        if body:
            nums = tuple(FbsNumbering(tuple(sorted(n.numbering, key=n.numbering.get)),
                                      tuple(str(n.numbering[a]) for a in
                                            sorted(n.numbering, key=n.numbering.get)))
                         for n in nodes)
            return "bridged_ring", nodes[0], body, nums, nodes[0].descriptor
    tree = decompose_fused_system(info, sub)
    if tree is None:
        return None
    from namepredict.layer5.fused_namer import fused_parent_names
    names = fused_parent_names(info["mol"], tree)
    if not names or not names[0] or not names[1]:
        return None  # 稠合树渲染不出基名（基不是稠合零件）：显式失败而非错名
    chain, labs = _fused_component_numbering(info["mol"], None, sub_rings, spiros, sub_edges)
    if not chain or not labs:
        return None
    extra = tuple("".join(_BRACKET.findall(names[0])))
    return ("fused_ring", tree, names, (FbsNumbering(tuple(chain), tuple(labs)),), extra)


def _retained_base(mol, sid: str, chain) -> tuple:
    """保留名组分基名（含 'a' 位次前缀，同 pack_parent_stem 口径）。"""
    from namepredict.layer2.kind_registry import pack_parent_stem
    packed = pack_parent_stem({"scaffold_id": sid, "mol": mol, "chain": list(chain or ())}, mol)
    return packed.get("stem_en"), packed.get("stem_zh")


def _cyclo_base(bare: bool):
    """单环烃组分基名生成器：全名 / 裸词干（可接 a-ene）。"""
    def make(n_ring: int):
        from namepredict.layer5.stems import alkane_en, alkane_zh
        en, zh = alkane_en(n_ring), alkane_zh(n_ring)
        if not en or not zh:
            return None
        if bare:
            return f"cyclo{en[:-3]}", f"环{zh[:-1]}"
        return f"cyclo{en}", f"环{zh}"
    return make


def _build_component(mol, info, system, comp: list[int], free, k: int) -> FbsComponent | None:
    """识别单个组分：基名 + 编号候选 + 螺原子（不可命名返回 None）。"""
    rings = list(sssr_rings(mol))
    atoms = tuple(sorted({a for r in comp for a in rings[r]}))
    sub_rings = [rings[r] for r in sorted(comp)]
    sub_edges = _component_system(system, comp, atoms)["fusion_edges"]
    spiros = tuple(sorted(a for a in atoms if a in free))
    sid = match_retained(info, atoms)
    bare_en = bare_zh = None
    if sid is not None:  # 保留名组分：杂原子与不饱和已含在名内
        from namepredict.layer2.ring_scaffold import _match_with_map
        hit = _match_with_map(info, atoms)
        match = tuple(hit[1]) if hit else ()
        matches = _retained_matches(info, sid, atoms, match)
        kind = "mono_ring" if len(comp) == 1 else "fused_ring"
        nums = _retained_numberings(mol, sid, matches, sub_rings, spiros, sub_edges)
        base_en, base_zh = _retained_base(mol, sid, nums[0].chain if nums else ())
        extra = () if kind == "mono_ring" else tuple("".join(_BRACKET.findall(base_en or "")))
        node = None
    elif len(comp) == 1:  # 未注册单环：碳环出 cyclo 基名，杂环本版不支持
        chain = tuple(sub_rings[0])
        if any(mol.GetAtomWithIdx(a).GetAtomicNum() != C for a in chain):
            return None
        kind, sid, node, extra = "mono_ring", None, None, ()
        nums = _mono_numberings(mol, chain)
        base_en, base_zh = _cyclo_base(False)(len(chain)) or (None, None)
        bare_en, bare_zh = _cyclo_base(True)(len(chain)) or (None, None)
    else:
        flat = _flat_numberings(info, system, comp, atoms, spiros, sub_rings, sub_edges)
        if flat is None:
            return None
        kind, node, (base_en, base_zh), nums, extra = flat
        if kind == "bridged_ring":  # 'a' 前缀按 P-24.5.2 提到 spiro 之前，故取裸词干
            from namepredict.layer5.bridged_namer import bridged_body_names
            body = bridged_body_names(node)
            if body is None:
                return None
            base_en, base_zh = body[0]
            bare_en, bare_zh = body[1]
    if not nums or base_en is None or base_zh is None:
        return None
    return FbsComponent(k, kind, atoms, tuple(sorted(comp)), spiros, base_en, base_zh,
                        bare_en, bare_zh, sid, node, tuple(nums), extra)


# ── 引用序（P-24.5.3 / P-24.6 / P-24.7.2） ──────────────

def _order_key(comp: FbsComponent) -> tuple:
    """引用序键：字母数字序 → 斜体稠合字母/von Baeyer 描述符（P-24.5.3）。"""
    from namepredict.tools.re import alpha_order_key
    return (alpha_order_key(comp.base_en or ""), comp.order_extra)


def citation_order(comps: list[FbsComponent], edges: list[list[tuple[int, int]]]):
    """P-24.5.1/24.6/24.7 定组分引用序，返回按引用项分组的组分下标（失败返回 None）。

    每个引用项是一个元组：长度 1 为单组分，>1 为同名端组分的倍增词组
    （P-24.7.1 的 tris(...) 与 P-24.7.2 的 bis(...)）。
    """
    n = len(comps)
    deg = [len(edges[k]) for k in range(n)]
    term = [k for k in range(n) if deg[k] == 1]
    if n < 2 or not term:
        return None
    if n == 2:  # P-24.5.1 两组分：字母数字序在前者列首
        return tuple((k,) for k in sorted(range(2), key=lambda k: _order_key(comps[k])))
    if all(d <= 2 for d in deg):  # P-24.6 直链：从字母序较小的端组分沿链展开
        start = min(term, key=lambda k: _order_key(comps[k]))
        out, prev, cur = [start], None, start
        while len(out) < n:
            nxt = [j for _, j in edges[cur] if j != prev]
            if len(nxt) != 1:
                return None
            prev, cur = cur, nxt[0]
            out.append(cur)
        return tuple((k,) for k in out)
    center = [k for k in range(n) if deg[k] >= 3]
    if len(center) != 1 or any(deg[k] != 1 for k in range(n) if k != center[0]):
        return None  # 更深的支链归 P-24.7.4，本版不做
    c = center[0]
    others = sorted((k for k in term if k != c), key=lambda k: _order_key(comps[k]))
    if len(others) < 2:
        return None
    if len({comps[k].base_en or str(k) for k in others}) == 1:  # P-24.7.1 端环全同：中心列首
        return ((c,), tuple(others))
    return _terminals_then_center(comps, c, others)  # P-24.7.2 端→心→端


def _terminals_then_center(comps: list[FbsComponent], center: int, others: list[int]):
    """P-24.7.2：字母序最早的端组分先列，随后中心，其余端组分按字母序（同名合并为组）。"""
    by_name: dict[str, list[int]] = {}
    for k in others:
        by_name.setdefault(comps[k].base_en or str(k), []).append(k)
    names = sorted(by_name, key=lambda nm: _order_key(comps[by_name[nm][0]]))
    return ((tuple(by_name[names[0]]), (center,))
            + tuple(tuple(by_name[nm]) for nm in names[1:]))


# ── 公共入口 ────────────────────────────────────────────

def decompose_fbs_system(info: dict, system: dict) -> list[FbsNode]:
    """公共入口：多环组分螺环系拆解为 FbsNode（不可命名返回空表）。"""
    free = tuple(system.get("free_spiro_atoms") or ())
    if not free:
        return []
    mol = info["mol"]
    rings = list(sssr_rings(mol))
    comps: list[FbsComponent] = []
    for k, comp in enumerate(component_indices(system)):
        built = _build_component(mol, info, system, comp, frozenset(free), k)
        if built is None:
            return []  # 任一组分不可命名即显式失败，不得静默丢组分
        comps.append(built)
    edges: list[list[tuple[int, int]]] = [[] for _ in comps]
    for i, j, s in system.get("spiro_edges") or ():
        if s not in free:
            continue
        a = next((c.index for c in comps if i in c.ring_indices), None)
        b = next((c.index for c in comps if j in c.ring_indices), None)
        if a is None or b is None or a == b:
            return []
        edges[a].append((s, b))
        edges[b].append((s, a))
    if sum(len(e) for e in edges) != 2 * (len(comps) - 1):
        return []  # 组分图须为树：边数 = 组分数 - 1
    groups = citation_order(comps, edges)
    if groups is None:
        return []
    flat = [k for grp in groups for k in grp]
    slot_of = {k: pos for pos, k in enumerate(flat)}
    ranked = tuple(replace(comps[k], index=slot_of[k]) for k in flat)  # index = 槽位，即撇号数
    placed = {flat[0]}
    links: list[tuple[int, int, int]] = []
    for pos, cur in enumerate(flat[1:], start=1):  # 每个后列组分须接在已列组分上
        hit = next(((s, j) for s, j in edges[cur] if j in placed), None)
        if hit is None:
            return []
        placed.add(cur)
        links.append((hit[0], slot_of[hit[1]], pos))
    cites = tuple(tuple(slot_of[k] for k in grp) for grp in groups)
    atom_ids = tuple(sorted({a for c in ranked for a in c.atom_ids}))
    ring = "carbo" if all(mol.GetAtomWithIdx(a).GetAtomicNum() == C for a in atom_ids) \
        else "hetero"
    return [FbsNode(atom_ids, free, ranked, tuple(links), cites, ring,
                    len(system.get("sssr_indices") or ()))]


