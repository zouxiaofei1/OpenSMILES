"""P-25.3.2.3 优选取向: 水平行摆放 + 四象限加权计数 + 取向选择。"""
from __future__ import annotations

import math
from dataclasses import dataclass

from namepredict.layer4.ring_geometry import (
    DEFORM_MAX,
    OVERLAP_FRAC,
    RING_TEMPLATES,
    apply_rigid,
    centroid,
    overlap_area,
    polygon_area,
    regular_polygon,
    rigid_fit,
    ring_cyclic,
)


@dataclass(frozen=True)
class Orientation:
    row: tuple[int, ...]                       # 水平行环索引(左→右)
    coords: tuple[tuple[int, float, float], ...]  # (原子, x, y) 元组(可哈希)
    quad: tuple[float, float, float, float]    # (Q1右上, Q2左上, Q3左下, Q4右下) 面积分数
    above: float                               # 水平轴上方环面积分数

    def coord_dict(self) -> dict[int, tuple[float, float]]:
        """转 {原子: (x,y)} dict。"""
        return {a: (x, y) for a, x, y in self.coords}


def _edge_map(fusion_edges) -> dict[frozenset[int], tuple[int, int]]:
    """共享 ≥2 原子融合边 → {环对: 共享原子对(排序)}。"""
    out = {}
    for i, j, sh in fusion_edges:
        if len(sh) >= 2:
            out[frozenset((i, j))] = tuple(sorted(sh)[:2])
    return out


def opposite_bonds(ring: tuple[int, ...], bond: tuple[int, int]) -> list[tuple[int, int]]:
    """偶环严格对面 1 条边; 奇环对面两条(P-25.3.2.3.1 奇数环 2 模板)。

    先由端点 index 定边位置 k(边为 (k, k+1) 或 (k+1, k)), 再取对面边。
    """
    n = len(ring)
    i, j = ring.index(bond[0]), ring.index(bond[1])
    k = i if (i + 1) % n == j else j
    if n % 2 == 0:
        kk = (k + n // 2) % n
        return [(ring[kk], ring[(kk + 1) % n])]
    k1 = (k + (n - 1) // 2) % n
    k2 = (k + (n + 1) // 2) % n
    return [(ring[k1], ring[(k1 + 1) % n]), (ring[k2], ring[(k2 + 1) % n])]


def _extend_row(rings, adj, left: int, mid: int, shared) -> tuple[int, ...]:
    """从 left→mid 方向递归向右扩展水平行。"""
    row = [left, mid]
    cur, entry = mid, frozenset(shared)
    while True:
        found = False
        for exit_bond in opposite_bonds(rings[cur], tuple(entry)):
            for nb, sh in adj[cur]:
                if nb not in row and frozenset(exit_bond) == sh:
                    row.append(nb)
                    cur, entry = nb, sh
                    found = True
                    break
            if found:
                break
        if not found:
            break
    return tuple(row)


def horizontal_rows(rings, fusion_edges) -> list[tuple[int, ...]]:
    """全部水平行候选(左→右), 按长度降序。"""
    adj = {r: [] for r in range(len(rings))}
    for i, j, sh in fusion_edges:
        if len(sh) >= 2:
            adj[i].append((j, frozenset(sh)))
            adj[j].append((i, frozenset(sh)))
    rows = set()
    for r in range(len(rings)):
        for nb, sh in adj[r]:
            rows.add(_extend_row(rings, adj, r, nb, sh))
    return sorted(rows, key=len, reverse=True)


def _place_ring(order: list[int], coords: dict, a: int, b: int, side: int) -> dict[int, tuple[float, float]]:
    """以共享边(a,b)摆放正 n 边形, 顶点0→a 顶点1→b, 心在边侧(side=+1 左法向, -1 右法向)。"""
    n = len(order)
    ax, ay = coords[a]
    bx, by = coords[b]
    length = math.hypot(bx - ax, by - ay)
    scale = length / (2 * math.sin(math.pi / n))
    ux, uy = (bx - ax) / length, (by - ay) / length   # 共享边方向
    vx, vy = -uy, ux                                  # 左法向
    out = {}
    for k, atom in enumerate(order):
        tx = math.cos(-math.pi / n + 2 * math.pi * k / n)
        ty = math.sin(-math.pi / n + 2 * math.pi * k / n)
        t = (ty + math.sin(math.pi / n)) * scale        # 沿共享边分量(顶点0→1 递增)
        nn = (math.cos(math.pi / n) - tx) * scale       # 垂直分量(左法向)
        if side < 0:
            nn = -nn
        out[atom] = (ax + t * ux + nn * vx, ay + t * uy + nn * vy)
    return out


def _opposite_side(coords: dict, prev_ring: tuple, a: int, b: int) -> int:
    """新环相对共享边(a,b)应取的 side: 使新环质心落在 prev 环质心的对面。

    共享边端点(a,b)顺序来自原子索引排序(_edge_map), 不保证几何方向, 故不能
    固定 side; 由 prev 环质心相对共享边左法向(-uy,ux)的符号决定(side=+1 时
    新环质心在左法向侧, -1 在右法向侧)。
    """
    ax, ay = coords[a]
    bx, by = coords[b]
    mx, my = (ax + bx) / 2, (ay + by) / 2
    cx = sum(coords[p][0] for p in prev_ring) / len(prev_ring)
    cy = sum(coords[p][1] for p in prev_ring) / len(prev_ring)
    d = (cx - mx) * (-(by - ay)) + (cy - my) * (bx - ax)
    return -1 if d > 0 else 1


def _layout(row: tuple[int, ...], rings, fusion_edges) -> dict | None:
    """刚性摆放水平行及其邻接环, 返回 {原子: 坐标}; 失败(奇环行内双侧)返回 None。"""
    edge = _edge_map(fusion_edges)
    coords: dict[int, tuple[float, float]] = {}
    placed: set[int] = set()
    for idx, r in enumerate(row):
        n = len(rings[r])
        if idx == 0:
            pair = edge.get(frozenset((r, row[1])))
            order = ring_cyclic(rings[r], pair[0], pair[1])
            tmpl = RING_TEMPLATES[n]
            for k, atom in enumerate(order):
                coords[atom] = tmpl[k]
        else:
            prev = row[idx - 1]
            pair = edge.get(frozenset((r, prev)))
            exit_pair = edge.get(frozenset((r, row[idx + 1]))) if idx < len(row) - 1 else None
            if n % 2 == 1 and exit_pair is not None:
                return None  # 行内奇环双侧融合需松弛, 本阶段刚性不可行
            order = ring_cyclic(rings[r], pair[0], pair[1])
            side = _opposite_side(coords, rings[prev], pair[0], pair[1])
            coords.update(_place_ring(order, coords, pair[0], pair[1], side))
        placed.add(r)
    for r in range(len(rings)):
        if r not in placed and not _place_neighbor(r, rings, coords, placed, edge):
            return None
    return coords


def _overlaps_any(cand: dict, rings, r: int, coords: dict) -> bool:
    """新环与任意环(含共享环)的重叠面积是否超过阈值。

    共享环正确摆放时只在共享边处重合(面积≈0), 若新环被摆到共享环同侧则会
    完全重合; 故不能跳过共享环——否则重合坐标无法被检测。"""
    cand_pts = [cand[a] for a in rings[r]]
    cand_area = polygon_area(cand_pts)
    for other, o_pts in _ring_polys(rings, coords):
        if other == r:
            continue
        inter = overlap_area(cand_pts, o_pts)
        if inter > OVERLAP_FRAC * min(cand_area, polygon_area(o_pts)):
            return True
    return False


def _ring_polys(rings, coords):
    """已摆环的多边形点集 [(环索引, [顶点])]。"""
    return [(r, [coords[a] for a in ring]) for r, ring in enumerate(rings)
            if all(a in coords for a in ring)]


def _place_neighbor(r: int, rings, coords: dict, placed: set[int], edge: dict) -> bool:
    """递归摆放环 r(至少一个共享边已定坐标); 返回是否成功。"""
    # edge 的 key 是 frozenset, 解包 (i, j) 顺序不可靠, 故用集合成员关系确定
    # 邻环索引(原实现解包顺序若翻转为 (j=r, i=邻环) 时 j 即 r 自身)。
    cands = []
    for key, pair in edge.items():
        if r in key:
            nb = next(x for x in key if x != r)
            if nb in placed:
                cands.append((nb, pair))
    if not cands:
        return False
    nb, (a, b) = cands[0]
    order = ring_cyclic(rings[r], a, b)
    # 首选共享环对侧(与水平行摆放一致的几何判定), 避免默认左法向把新环摆到
    # 与共享环同侧导致完全重合; 失败再试对侧(_overlaps_any 现含共享环检查)。
    side = _opposite_side(coords, rings[nb], a, b)
    best = None
    for cand_side in (side, -side):
        cand = _place_ring(order, coords, a, b, cand_side)
        if not _overlaps_any(cand, rings, r, coords):
            best = cand
            break
    if best is None:
        return False
    coords.update(best)
    placed.add(r)
    return True


def _ring_deform(pts: list, n: int) -> float:
    """环坐标相对正 n 边形模板的最大偏差(循环移位+镜像最优对齐)。"""
    tmpl = RING_TEMPLATES[n]
    best = float("inf")
    for shift in range(n):
        fwd = pts[shift:] + pts[:shift]
        for rev in (False, True):
            order = list(reversed(fwd)) if rev else fwd
            params = rigid_fit(order, tmpl)
            fitted = [apply_rigid(p, params) for p in order]
            best = min(best, max(math.dist(fitted[k], tmpl[k]) for k in range(n)))
    return best


def _valid_deform_overlap(coords: dict, rings, fusion_edges) -> bool:
    """每环刚体拟合偏差与环间重叠检查。"""
    for ring in rings:
        pts = [coords[a] for a in ring]
        dev = _ring_deform(pts, len(ring))
        length = max(math.dist(pts[k], pts[(k + 1) % len(pts)]) for k in range(len(pts)))
        if length > 1e-9 and dev / length > DEFORM_MAX:
            return False
    polys = _ring_polys(rings, coords)
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            ri, pi = polys[i]
            rj, pj = polys[j]
            if frozenset(rings[ri]) & frozenset(rings[rj]):
                continue
            inter = overlap_area(pi, pj)
            if inter > OVERLAP_FRAC * min(polygon_area(pi), polygon_area(pj)):
                return False
    return True


def _quadrant_fractions(coords: dict, rings) -> tuple[tuple[float, float, float, float], float]:
    """四象限(Q1右上,Q2左上,Q3左下,Q4右下)与水平轴上方环面积分数。"""
    cx, cy = centroid(list(coords.values()))
    big = 1e6
    q = [0.0, 0.0, 0.0, 0.0]
    above = 0.0
    for ring in rings:
        pts = [coords[a] for a in ring]
        area = polygon_area(pts)
        if area < 1e-12:
            continue
        from namepredict.layer4.ring_geometry import clip_polygon
        # 各象限裁剪矩形(以 cx,cy 为角, 逆时针; 足够大覆盖整个环系)。
        corners = [
            [(cx, cy), (cx + big, cy), (cx + big, cy + big), (cx, cy + big)],  # Q1 右上
            [(cx, cy), (cx, cy + big), (cx - big, cy + big), (cx - big, cy)],  # Q2 左上
            [(cx, cy), (cx - big, cy), (cx - big, cy - big), (cx, cy - big)],  # Q3 左下
            [(cx, cy), (cx, cy - big), (cx + big, cy - big), (cx + big, cy)],  # Q4 右下
        ]
        for k, clip in enumerate(corners):
            inter = clip_polygon(pts, clip)
            q[k] += polygon_area(inter) / area
        top = clip_polygon(pts, [(cx - big, cy), (cx + big, cy), (cx + big, cy + big), (cx - big, cy + big)])
        above += polygon_area(top) / area
    return tuple(q), above


def preferred_orientations(mol, rings, fusion_edges) -> list[Orientation]:
    """全部优选取向平局候选(水平行环数最多→右上最多→左下最少→上方最多)。

    对称环系(直线 acene 等)的左右/上下镜像 key 相同, 全部返回, 由编号阶段
    P-25.3.3.1.2 准则(a)-(d) 跨候选收窄; 否则依赖遍历顺序, 编号不稳定。
    """
    rows = horizontal_rows(rings, fusion_edges)
    print(rows,rings,fusion_edges,"\n")
    
    if not rows:
        return []
    max_len = len(rows[0])
    best_key: tuple | None = None
    bests: list[Orientation] = []
    for row in rows:
        if len(row) < max_len:
            continue
        for flip in (False, True):
            coords = _layout(row, rings, fusion_edges)
            # print("coords:",coords)
            if coords is None:
                continue
            if flip:
                coords = {a: (x, -y) for a, (x, y) in coords.items()}
            if not _valid_deform_overlap(coords, rings, fusion_edges):
                continue
            (q1, q2, q3, q4), above = _quadrant_fractions(coords, rings)
            # 面积分数有 ~1e-15 浮点尾差, round 消除后镜像才算平局
            key = (len(row), round(q1, 9), round(-q3, 9), round(above, 9))
            orient = Orientation(row, tuple((a, x, y) for a, (x, y) in sorted(coords.items())),
                                 (q1, q2, q3, q4), above)
            if best_key is None or key > best_key:
                best_key = key
                bests = [orient]
            elif key == best_key:
                bests.append(orient)
    print(bests)
    return bests


def preferred_orientation(mol, rings, fusion_edges) -> Orientation | None:
    """优选取向(取首个平局候选); 无候选返回 None。"""
    bests = preferred_orientations(mol, rings, fusion_edges)
    return bests[0] if bests else None
