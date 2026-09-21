"""P-24 螺环拆解：环组分划分 + von Baeyer 螺描述符 + 并列编号候选。

螺描述符本质是螺图上的一条欧拉回路：顶点为螺原子，边为环内两螺原子之间的
原子段。L2 只枚举候选（起始端环 × 各段取向），编号裁决交 L4 按 P-24.2.2.1/2
与 P-24.2.4.1.2 收窄。L5 不得 import L2，故 node 只按鸭子类型读
descriptor / descriptor_superscripts / free_spiro_atoms / numbering。
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from namepredict.layer1.ring_systems import sssr_rings

CARBON = 6
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
    return "carbo" if all(mol.GetAtomWithIdx(a).GetAtomicNum() == CARBON
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
        from namepredict.layer2.fbs_system import decompose_fbs_system
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


def has_spiro_system(system: dict | None) -> bool:
    """环系是否含自由螺连接（P-24.1）。"""
    return bool(system is not None and system.get("free_spiro_atoms"))


def _numbers_all(node: SpiroNode, atom_ids) -> bool:
    """编号是否覆盖骨架全部原子且无重号。"""
    return len(node.numbering) == len(atom_ids) and set(node.numbering.values()) == set(
        range(1, len(atom_ids) + 1))
