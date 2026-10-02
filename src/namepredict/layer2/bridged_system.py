"""P-23 扩展 von Baeyer 桥环拆解：主环/主桥/次级桥 + 全原子位次。"""
from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import permutations

from namepredict.layer1.ring_systems import sssr_rings
from namepredict.layer2.ring_scaffold import _ring_kind
from namepredict.layer4.numbering_engine import narrow

_MAX_CYCLES = 256  # 桥图简单环枚举上限
_MAX_ORDERS = 120  # 次级桥引用顺序枚举上限


@dataclass(frozen=True)
class BridgeSegment:
    """一段桥：两端桥头 + 桥内原子（沿 heads[0]→heads[1] 序）。"""

    heads: tuple[int, int]
    atoms: tuple[int, ...] = ()
    numbers: tuple[int, ...] = ()
    locants: tuple[int, int] = ()

    def __len__(self) -> int:
        """桥内原子数，即描述符中的数字（P-23.2.6.1.2）。"""
        return len(self.atoms)

    def reversed(self) -> "BridgeSegment":
        """掉转方向，使 atoms 沿 heads[1]→heads[0] 序。"""
        return replace(self, heads=(self.heads[1], self.heads[0]),
                       atoms=tuple(reversed(self.atoms)))


@dataclass(frozen=True)
class BridgedNode:
    """P-23 桥环母体：主环/主桥/次级桥划分 + 原子→位次映射。"""

    atom_ids: tuple[int, ...]
    numbering: dict[int, int]
    main_ring: tuple[int, ...]
    ring_segments: tuple[int, int]
    main_bridge: BridgeSegment
    independent_bridges: tuple[BridgeSegment, ...]
    dependent_bridges: tuple[BridgeSegment, ...]
    ring: str  # "carbo" / "hetero"
    n_rings: int

    @property
    def secondary_bridges(self) -> tuple[BridgeSegment, ...]:
        """按名称引用顺序：独立次级桥在前、依赖次级桥在后。"""
        return self.independent_bridges + self.dependent_bridges

    @property
    def descriptor(self) -> tuple[int, ...]:
        """von Baeyer 描述符：主环两段 + 主桥 + 各次级桥长度。"""
        return (self.ring_segments + (len(self.main_bridge),)
                + tuple(len(s) for s in self.secondary_bridges))

    @property
    def locant_pairs(self) -> tuple[tuple[int, int], ...]:
        """次级桥上标位次对（小者在前），按引用顺序。"""
        return tuple(s.locants for s in self.secondary_bridges)

    def scaffold_identity(self):
        """转为 L2 scaffold 身份（id 与命名类均为 bridged）。"""
        from namepredict.layer2.ring_scaffold import ScaffoldIdentity
        return ScaffoldIdentity("bridged", "bridged", self.n_rings, self.ring)


@dataclass(frozen=True)
class _Cand:
    """结构候选：主环环序 + 主桥头对 + 主桥边 + 次级桥边序。"""

    ring_v: tuple[int, ...]
    ring_e: tuple[int, ...]
    start: int
    end: int
    main_e: int
    sec_e: tuple[int, ...]


def _adj(mol, atoms: set[int]) -> dict[int, list[int]]:
    """环系内邻接表：只保留两端都在环系内的键。"""
    return {i: sorted(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()
                      if n.GetIdx() in atoms) for i in sorted(atoms)}


def _ring_count(adj: dict[int, list[int]]) -> int:
    """环数 = 键数 - 原子数 + 1（P-23.2.6.1.1）。"""
    return sum(len(v) for v in adj.values()) // 2 - len(adj) + 1


def _path_segment(adj, comp: set[int], hset: set[int]) -> BridgeSegment | None:
    """去桥头分量须为简单路径且两端各邻接一个桥头，返回该桥。"""
    if len(comp) == 1:
        v = next(iter(comp))
        hn = [n for n in adj[v] if n in hset]
        if len(hn) == 2 and hn[0] != hn[1]:
            return BridgeSegment((hn[0], hn[1]), (v,))
        return None
    inner = {v: [n for n in adj[v] if n in comp] for v in comp}
    if any(len(ns) > 2 for ns in inner.values()):
        return None
    ends = sorted(v for v, ns in inner.items() if len(ns) == 1)
    if len(ends) != 2:
        return None
    order, prev, cur = [ends[0]], None, ends[0]
    while cur != ends[1]:
        nxt = [n for n in inner[cur] if n != prev]
        if not nxt:
            return None
        prev, cur = cur, nxt[0]
        order.append(cur)
    if len(order) != len(comp):
        return None
    ha = [n for n in adj[ends[0]] if n in hset]
    hb = [n for n in adj[ends[1]] if n in hset]
    if len(ha) != 1 or len(hb) != 1 or ha[0] == hb[0]:
        return None
    return BridgeSegment((ha[0], hb[0]), tuple(order))


def _bridges(adj, heads: tuple[int, ...]) -> list[BridgeSegment] | None:
    """拆出全部桥：去桥头分量 + 桥头直键（0 原子桥，P-23.1.2）。"""
    hset, seen, out = set(heads), set(), []
    for v in sorted(adj):
        if v in hset or v in seen:
            continue
        comp, stack = set(), [v]
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            stack.extend(n for n in adj[cur] if n not in hset)
        seen |= comp
        seg = _path_segment(adj, comp, hset)
        if seg is None:
            return None
        out.append(seg)
    for i, a in enumerate(heads):
        for b in heads[i + 1:]:
            if b in adj[a]:
                out.append(BridgeSegment((a, b)))
    return out


def _h_graph(heads, bridges) -> dict[int, list[tuple[int, int]]]:
    """桥图 H：顶点=桥头，边=桥（可平行），返回 顶点→[(边序, 对端)]。"""
    g: dict[int, list[tuple[int, int]]] = {h: [] for h in heads}
    for i, seg in enumerate(bridges):
        g[seg.heads[0]].append((i, seg.heads[1]))
        g[seg.heads[1]].append((i, seg.heads[0]))
    return g


def _cycles(g, heads, limit: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    """H 的简单环：顶点序（首元素为环上最小顶点）+ 边序。"""
    out: list[tuple[tuple[int, ...], tuple[int, ...]]] = []

    def walk(start: int, cur: int, vs: list, es: list) -> None:
        if len(out) >= limit:
            return
        for eid, nxt in g[cur]:
            if eid in es:
                continue
            if nxt == start:
                if len(vs) >= 3 or len(es) == 1:  # 长环 / 两顶点平行边环
                    out.append((tuple(vs), tuple(es) + (eid,)))
                continue
            if nxt < start or nxt in vs:
                continue
            walk(start, nxt, vs + [nxt], es + [eid])

    for s in sorted(heads):
        walk(s, s, [s], [])
    return out


def _cycle_path(vs, es, u: int, v: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """主环上从 u 沿环序走到 v 的边序与顶点序。"""
    k, i, j = len(vs), vs.index(u), vs.index(v)
    edges, verts, cur = [], [u], i
    while cur != j:
        edges.append(es[cur])
        cur = (cur + 1) % k
        verts.append(vs[cur])
    return tuple(edges), tuple(verts)


def _path_atoms(vs, es, u, v, bridges) -> tuple[int, ...]:
    """主环上 u→v 链段内原子（不含两端桥头），按行进序。"""
    edges, verts = _cycle_path(vs, es, u, v)
    out: list[int] = []
    for e, (x, y) in zip(edges, zip(verts, verts[1:])):
        seg = bridges[e]
        out.extend(seg.atoms if seg.heads[0] == x else reversed(seg.atoms))
        if y != verts[-1]:
            out.append(y)
    return tuple(out)

def _secondary_ok(sec, bridges, core: set[int]) -> bool:
    """剩余桥须构成以 core 桥头为端点的简单链（P-23.1.7 / P-23.1.8）。"""
    deg: dict[int, int] = {}
    for i in sec:
        a, b = bridges[i].heads
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    if any(d != 2 for h, d in deg.items() if h not in core):
        return False
    return True


def _candidates(bridges, heads) -> list[_Cand]:
    """枚举全部 (主环, 主桥头对, 主桥) 结构候选（P-23.2.1 / P-23.2.6.2.1）。"""
    g, out = _h_graph(heads, bridges), []
    for vs, es in _cycles(g, heads, _MAX_CYCLES):
        rset, eset = set(vs), set(es)
        for m, seg in enumerate(bridges):
            u, v = seg.heads
            if m in eset or u not in rset or v not in rset:
                continue
            sec = tuple(i for i in range(len(bridges)) if i not in eset and i != m)
            if not _secondary_ok(sec, bridges, rset):
                continue
            out.append(_Cand(vs, es, u, v, m, sec))
            out.append(_Cand(vs, es, v, u, m, sec))
    return out


def _number_one(seg: BridgeSegment, cur: dict[int, int]) -> BridgeSegment:
    """从紧邻较高编号桥头的一端给桥原子编号（P-23.2.6.3）。"""
    a, b = seg.heads
    if cur.get(b, 0) > cur.get(a, 0):
        seg = seg.reversed()
    nxt = max(cur.values(), default=0)
    nums = []
    for atom in seg.atoms:
        nxt += 1
        cur[atom] = nxt
        nums.append(nxt)
    return replace(seg, numbers=tuple(nums),
                   locants=tuple(sorted((cur[seg.heads[0]], cur[seg.heads[1]]))))


def _number_secondary(core_num: dict[int, int], segs):
    """按给定引用顺序编号次级桥：独立在前、依赖在后（P-23.2.6.3）。"""
    cur, core = core_num, set(core_num)  # 就地续编号：新原子须写回调用方的映射
    ind = [s for s in segs if s.heads[0] in core and s.heads[1] in core]
    dep, out_dep = [s for s in segs if s not in ind], []
    out_ind = [_number_one(s, cur) for s in ind]
    while dep:
        ready = [s for s in dep if s.heads[0] in cur or s.heads[1] in cur]
        if not ready:
            return None
        s = ready[0]  # 保持给定顺序，仅跳过尚不可编号者
        dep.remove(s)
        out_dep.append(_number_one(s, cur))
    return tuple(out_ind), tuple(out_dep)


def _to_node(cand: _Cand, bridges, atom_ids, ring: str, r: int) -> BridgedNode | None:
    """给候选编号并组装节点；位次冲突或覆盖不全返回 None。"""
    a_at = _path_atoms(cand.ring_v, cand.ring_e, cand.start, cand.end, bridges)
    b_at = _path_atoms(cand.ring_v, cand.ring_e, cand.end, cand.start, bridges)
    if len(a_at) >= len(b_at):  # 长边在前（P-23.2.3）
        seg1, seg2 = a_at, b_at
    else:  # 改从另一端绕行：两条边的行进方向都要掉转，否则 seq 不是合法路径
        seg1, seg2 = tuple(reversed(b_at)), tuple(reversed(a_at))
    seq = (cand.start,) + seg1 + (cand.end,) + seg2
    if len(set(seq)) != len(seq):
        return None
    core = {a: i + 1 for i, a in enumerate(seq)}
    mseg, nxt, mnums = bridges[cand.main_e], len(seq), []
    main_atoms = mseg.atoms if mseg.heads[0] == cand.start else tuple(reversed(mseg.atoms))
    for a in main_atoms:
        nxt += 1
        core[a] = nxt
        mnums.append(nxt)
    main = BridgeSegment((cand.start, cand.end), main_atoms, tuple(mnums),
                         tuple(sorted((core[cand.start], core[cand.end]))))
    sec_segs = [bridges[i] for i in cand.sec_e]
    best = None
    for order in _orders(sec_segs):  # 引用顺序由 P-23.2.6.2.2/2.5 裁决，逐个试
        num = dict(core)
        sec = _number_secondary(num, list(order))
        if sec is None or set(num) != set(atom_ids):
            continue
        ind, dep = sec
        node = BridgedNode(atom_ids=tuple(sorted(atom_ids)), numbering=num, main_ring=seq,
                           ring_segments=(len(seg1), len(seg2)), main_bridge=main,
                           independent_bridges=ind, dependent_bridges=dep, ring=ring, n_rings=r)
        key = (tuple(-len(s) for s in node.secondary_bridges), _locants(node))
        if best is None or key < best[0]:
            best = (key, node)
    return best[1] if best else None


def _orders(segs: list):
    """次级桥引用顺序候选：长桥在前由 P-23.2.6.2.2 定，同长度组内交 P-23.2.6.2.5 裁决。"""
    if len(segs) < 2:
        return [tuple(segs)]
    groups: dict[int, list] = {}
    for seg in segs:
        groups.setdefault(len(seg), []).append(seg)
    out: list[tuple] = [()]
    for size in sorted(groups, reverse=True):
        out = [prefix + perm for prefix in out for perm in permutations(groups[size])]
        if len(out) > _MAX_ORDERS:  # 并列过多时退化为单一顺序
            return [tuple(sorted(segs, key=lambda s: -len(s)))]
    return out


def _locants(node: BridgedNode) -> tuple[int, ...]:
    """次级桥上标位次按引用顺序展平（P-23.2.6.2.5）。"""
    return tuple(n for pair in node.locant_pairs for n in pair)


def narrow_candidates(nodes: list[BridgedNode]) -> list[BridgedNode]:
    """按 P-23.2.6.2.1~2.5 依次收窄；末步不做确定性 tie-break。"""
    steps = [
        (lambda n: sum(n.ring_segments) + 2, True),          # P-23.2.1 主环原子数最大
        (lambda n: len(n.main_bridge), True),                # P-23.2.4 主桥最长（先于对称）
        (lambda n: abs(n.ring_segments[0] - n.ring_segments[1]), False),  # P-23.2.6.2.1 分割最对称
        (lambda n: tuple(len(s) for s in n.independent_bridges), True),   # P-23.2.6.2.2
        (lambda n: len(n.dependent_bridges), False),         # P-23.2.6.2.3 依赖桥最少
        (lambda n: tuple(sorted(_locants(n))), False),       # P-23.2.6.2.4 位次集合最低
        (lambda n: _locants(n), False),                      # P-23.2.6.2.5 引用序位次最低
    ]
    for key_fn, reverse in steps:
        nodes = narrow(nodes, key_fn, reverse=reverse)
    return nodes


def decompose_bridged_system(info: dict, system: dict) -> list[BridgedNode]:
    """公共入口：环系按 P-23 拆解为并列最优 BridgedNode 候选，非桥环返回空表。"""
    mol = info["mol"]
    atom_ids = tuple(sorted(system.get("atom_ids") or ()))
    if len(atom_ids) < 4:
        return []
    adj = _adj(mol, set(atom_ids))
    r = _ring_count(adj)
    if r < 2:
        return []
    heads = tuple(v for v in atom_ids if len(adj[v]) >= 3)  # P-23.1.1 桥头
    if len(heads) < 2:
        return []
    bridges = _bridges(adj, heads)
    if bridges is None or len(bridges) - len(heads) + 1 != r:
        return []
    ring = _ring_kind(mol, atom_ids)
    nodes = [n for c in _candidates(bridges, heads)
             if (n := _to_node(c, bridges, atom_ids, ring, r)) is not None]
    nodes = [n for n in nodes if sum(n.descriptor) + 2 == len(n.atom_ids)]  # P-23.2.6.1.4
    return narrow_candidates(nodes)


def _p25_names_ok(mol, fused_tree) -> bool:
    """P25 能否拼出稠合词干：直接问 L5 组装器，避免在 L2 复刻判据而漂移。"""
    from namepredict.layer5.fused_namer import fused_parent_names
    names = fused_parent_names(mol, fused_tree)
    return bool(names and names[0] and names[1])


def _cannot_be_mancude(mol, rings, atom_ids) -> bool:
    """组分环能否写出最大非累积双键（P-25.3.1.2）。

    mancude 性是环骨架的性质，与取代基无关：只看环系内 σ 键数。
    环系内已占满 4 根 σ 键的碳（4x 单键的环交界原子）腾不出 π 键。
    """
    atoms = set(atom_ids)
    for ring in rings:
        for i in ring:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetAtomicNum() != 6:
                continue
            if sum(1 for n in atom.GetNeighbors() if n.GetIdx() in atoms) >= 4:
                return True
    return False


def _fusion_naming_applies(mol, rings, atom_ids) -> bool:
    """稠合命名法是否适用（P-52.2.4.1 五元环要求 + P-25.3.1.2 mancude 要求）。"""
    if sum(1 for r in rings if len(r) >= 5) < 2:
        return False  # 不足两个五元或更多元环：von Baeyer 才是 PIN
    return not _cannot_be_mancude(mol, rings, atom_ids)


def _tree_ring_indices(node) -> set[int]:
    """稠环树覆盖的 SSSR 环下标（含全部附加组分）。"""
    out = set(node.ring_indices)
    for child in node.attached:
        out |= _tree_ring_indices(child)
    return out


def _tree_covers_rings(fused_tree, system: dict) -> bool:
    """稠环树是否覆盖环系的全部环（少一个环的稠合名会静默丢环）。"""
    if fused_tree is None:
        return False
    return _tree_ring_indices(fused_tree) == set(system.get("sssr_indices") or ())


def try_bridged_scaffold(info: dict, scaffold, fused_tree, system: dict) -> list[BridgedNode]:
    """路由到桥环：稠合命名法不适用，或 P25 命名失败。"""
    if system.get("free_spiro_atoms"):
        return []  # 含自由螺连接的环系归 P-24，不按 P-23 桥环命名
    nodes = decompose_bridged_system(info, system)
    if not nodes:
        return []
    mol = info["mol"]
    atom_ids = system.get("atom_ids") or ()
    sssr = system.get("sssr_indices") or ()
    all_rings = sssr_rings(mol)
    rings = [all_rings[i] for i in sssr]
    if not _fusion_naming_applies(mol, rings, atom_ids):  # 条件A：稠合不适用则直接走桥环
        return nodes
    if scaffold is None or scaffold.id != "fused_hetero":
        return []  # 保留模板身份与环系同构（P-25.2）→ 保持既有稠环路径
    if not _tree_covers_rings(fused_tree, system):
        return nodes  # 稠环树盖不住环系的全部环：稠合名必丢环
    if _p25_names_ok(mol, fused_tree):
        return []  # 条件B：P25 能产名则保持既有稠环路径
    return nodes
