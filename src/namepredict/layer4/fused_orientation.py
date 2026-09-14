"""P-25.3.2.3 优选取向: 水平行摆放 + 四象限加权计数 + 取向选择。"""
from __future__ import annotations

import math
from dataclasses import dataclass

from namepredict.layer4.ring_geometry import (
    DEFORM_MAX,
    OVERLAP_FRAC,
    RING_TEMPLATES,
    apply_rigid,
    overlap_area,
    polygon_area,
    rigid_fit,
    ring_cyclic,
    ring_shape_template,
)


@dataclass(frozen=True)
class Orientation:
    """P-25.3.2.3 取向候选：水平行环序、原子坐标、四象限与水平轴上方环计数。"""
    row: tuple[int, ...]                       # 水平行环索引(左→右)
    coords: tuple[tuple[int, float, float], ...]  # (原子, x, y) 元组(可哈希)
    quad: tuple[float, float, float, float]    # (Q1右上, Q2左上, Q3左下, Q4右下) 面积分数

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
    """偶环对面 1 条边、奇环对面两条(P-25.3.2.3.1)。"""
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


def _place_ring(order: list[int], coords: dict, a: int, b: int, side: int,
                template: list[tuple[float, float]] | None = None,
                used_tmpl: list | None = None) -> dict[int, tuple[float, float]]:
    """以共享边(a,b)摆放环(顶点0→a、1→b)；template 为环模板。"""
    n = len(order)
    ax, ay = coords[a]
    bx, by = coords[b]
    length = math.hypot(bx - ax, by - ay)
    ux, uy = (bx - ax) / length, (by - ay) / length   # 共享边方向
    vx, vy = -uy, ux                                  # 左法向
    if template is None:
        scale = length / (2 * math.sin(math.pi / n))
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
    base = [(tx, -ty) for tx, ty in template] if side < 0 else template  # side<0 用模板镜像使心朝右法向；共享边端点不动。
    p0, p1 = base[0], base[1]
    t0len = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
    scale = length / t0len
    theta = math.atan2(by - ay, bx - ax) - math.atan2(p1[1] - p0[1], p1[0] - p0[0])
    c, s = math.cos(theta), math.sin(theta)
    out = {}
    for k, atom in enumerate(order):
        dx, dy = base[k][0] - p0[0], base[k][1] - p0[1]
        out[atom] = (ax + scale * (c * dx - s * dy), ay + scale * (s * dx + c * dy))
    if used_tmpl is not None:
        used_tmpl[:] = [base]
    return out


def _opposite_side(coords: dict, prev_ring: tuple, a: int, b: int) -> int:
    """按 prev 环质心在左法向的符号定 side，使新环落在 prev 对面。"""
    ax, ay = coords[a]
    bx, by = coords[b]
    mx, my = (ax + bx) / 2, (ay + by) / 2
    cx = sum(coords[p][0] for p in prev_ring) / len(prev_ring)
    cy = sum(coords[p][1] for p in prev_ring) / len(prev_ring)
    d = (cx - mx) * (-(by - ay)) + (cy - my) * (bx - ax)
    return -1 if d > 0 else 1


def _layout(row: tuple[int, ...], rings, fusion_edges) -> tuple[dict | None, dict] | None:
    """摆放水平行及邻接环，返回 (坐标, 每环模板)；失败返回 None。"""
    edge = _edge_map(fusion_edges)
    coords: dict[int, tuple[float, float]] = {}
    ring_templates: dict[int, list] = {}
    placed: set[int] = set()
    for idx, r in enumerate(row):#遍历水平行每个环
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
            order = ring_cyclic(rings[r], pair[0], pair[1])
            tpl = None
            if n % 2 == 1 and exit_pair is not None:  # P-25.3.2.3.2：奇环行内双侧融合 → 变形环模板，仅当可构造时才放行
                exit_idx = next(
                    (k for k in range(n) if frozenset((order[k], order[(k + 1) % n])) == frozenset(exit_pair)),
                    None,
                )
                tpl = ring_shape_template(order, exit_idx) if exit_idx is not None else None
                if tpl is None:
                    return None  # 无法构造变形环 → 弃行(原逻辑, 不回归)
                ring_templates[r] = tpl
            side = _opposite_side(coords, rings[prev], pair[0], pair[1])
            used_tmpl: list = []
            coords.update(_place_ring(order, coords, pair[0], pair[1], side, tpl, used_tmpl))
            if used_tmpl:
                ring_templates[r] = used_tmpl[0]
        placed.add(r)
    for r in range(len(rings)):
        if r not in placed and not _place_neighbor(r, rings, coords, placed, edge):
            return None
    return coords, ring_templates


def _overlaps_any(cand: dict, rings, r: int, coords: dict) -> bool:
    """新环与任意环(含共享环)重叠面积是否超阈值。"""
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
    cands = []  # edge key 是 frozenset, 解包顺序不可靠, 用成员关系定。
    for key, pair in edge.items():
        if r in key:
            nb = next(x for x in key if x != r)
            if nb in placed:
                cands.append((nb, pair))
    if not cands:
        return False
    nb, (a, b) = cands[0]
    order = ring_cyclic(rings[r], a, b)
    side = _opposite_side(coords, rings[nb], a, b)  # 首选共享环对侧(与水平行一致), 失败再试对侧。
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


def _ring_deform(pts: list, n: int, tmpl: list | None = None) -> float:
    """环坐标相对模板(默认正 n 边形)的最大偏差，循环移位+镜像最优对齐。"""
    tmpl = tmpl if tmpl is not None else RING_TEMPLATES[n]
    is_distorted = tmpl is not None and tmpl is not RING_TEMPLATES[n]
    best = float("inf")
    for shift in range(n):
        fwd = pts[shift:] + pts[:shift]
        for rev in (False, True):
            order = list(reversed(fwd)) if rev else fwd
            variants = [order]
            if is_distorted:  # 变形环对 Kabsch 旋转方向敏感，补试镜像 x 的对齐。
                variants.append([(-x, y) for x, y in order])
            for cand in variants:
                params = rigid_fit(cand, tmpl)
                fitted = [apply_rigid(p, params) for p in cand]
                best = min(best, max(math.dist(fitted[k], tmpl[k]) for k in range(n)))
    return best


def _valid_deform_overlap(coords: dict, rings, ring_templates: dict | None = None) -> bool:
    """每环刚体拟合偏差与环间重叠检查。"""
    ring_templates = ring_templates or {}
    for r, ring in enumerate(rings):
        pts = [coords[a] for a in ring]
        dev = _ring_deform(pts, len(ring), ring_templates.get(r))
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


def _row_center(coords: dict, rings, row: tuple[int, ...]) -> tuple[float, float]:
    """水平行中心(P-25.3.2.3.3(b))：偶数环取共同键，奇数环取中心环。"""
    n = len(row)
    if n % 2 == 0:
        a, b = row[n // 2 - 1], row[n // 2]      # 中间一对相邻环
        shared = set(rings[a]) & set(rings[b])   # 中心共同键(竖直共用边)
        pts = [coords[x] for x in shared]
        return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
    c = row[n // 2]                              # 中心环
    pts = [coords[a] for a in rings[c]]
    return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)


def _ring_quadrant_contrib(pts, cx: float, cy: float) -> tuple[float, float, float, float]:
    """单环四象限贡献：两轴平分各 1/4，一轴平分两侧各 1/2，否则整环归一象限。"""
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    cross_v = min(xs) < cx and max(xs) > cx      # 顶点严格跨竖直轴(仅触轴不算)
    cross_h = min(ys) < cy and max(ys) > cy      # 顶点严格跨水平轴
    gx = sum(xs) / len(xs)
    gy = sum(ys) / len(ys)
    q = [0.0, 0.0, 0.0, 0.0]
    if cross_v and cross_h:                      # 被两条轴平分 → 四象限各 1/4
        return (0.25, 0.25, 0.25, 0.25)
    if cross_h:                                  # 只被水平轴平分 → 上下各 1/2，所在左右由质心定
        if gx >= cx:                             # 右：上 Q1 / 下 Q4
            q[0] = q[3] = 0.5
        else:                                    # 左：上 Q2 / 下 Q3
            q[1] = q[2] = 0.5
        return tuple(q)
    if cross_v:                                  # 只被竖直轴平分 → 左右各 1/2，所在上下由质心定
        if gy >= cy:                             # 上：左 Q2 / 右 Q1
            q[1] = q[0] = 0.5
        else:                                    # 下：左 Q3 / 右 Q4
            q[2] = q[3] = 0.5
        return tuple(q)
    k = 0 if (gx >= cx and gy >= cy) else \
        1 if (gx < cx and gy >= cy) else \
        3 if (gx >= cx and gy < cy) else 2       # Q1/Q2/Q4/Q3
    q[k] = 1.0
    return tuple(q)


def _above_contrib(pts, cy: float) -> float:
    """单个环在水平轴上方环数：被水平轴平分计 1/2，完全在上方计 1，否则 0。"""
    ys = [p[1] for p in pts]
    if min(ys) < cy and max(ys) > cy:
        return 0.5
    return 1.0 if sum(ys) / len(ys) >= cy else 0.0


def _quadrant_fractions(coords: dict, rings, row) -> tuple[tuple[float, float, float, float], float]:
    """四象限与上方环数：环计数法，原点取水平行中心 (P-25.3.2.3.3)。"""
    cx, cy = _row_center(coords, rings, row)
    q = [0.0, 0.0, 0.0, 0.0]
    above = 0.0
    for ring in rings:
        pts = [coords[a] for a in ring]
        cr = _ring_quadrant_contrib(pts, cx, cy)
        for k in range(4):
            q[k] += cr[k]
        above += _above_contrib(pts, cy)
    return tuple(q), above


def preferred_orientations(mol, rings, fusion_edges) -> list[Orientation]:
    """全部优选取向平局候选(水平行环数→右上→左下→上方)，镜像一并返回。"""
    rows = horizontal_rows(rings, fusion_edges)  # print(rows,rings,fusion_edges,"\n")

    if not rows:
        return []
    max_len = len(rows[0])
    best_key: tuple | None = None
    bests: list[Orientation] = []
    for row in rows:
        if len(row) < max_len:
            continue
        laid = _layout(row, rings, fusion_edges)  # print("coords:",coords)  # 与 flip 无关
        if laid is None:
            continue
        for flip in (False, True):
            coords, ring_templates = laid
            if flip:
                coords = {a: (x, -y) for a, (x, y) in coords.items()}
            if not _valid_deform_overlap(coords, rings, ring_templates):
                continue
            (q1, q2, q3, q4), above = _quadrant_fractions(coords, rings, row)
            key = (len(row), q1, -q3, above)  # 象限/环数为离散值 0/0.25/0.5/1，无浮点累计尾差
            orient = Orientation(row, tuple((a, x, y) for a, (x, y) in sorted(coords.items())),
                                 (q1, q2, q3, q4))
            if best_key is None or key > best_key:
                best_key = key
                bests = [orient]
            elif key == best_key:
                bests.append(orient)
    return bests  # print(bests)


def preferred_orientation(mol, rings, fusion_edges) -> Orientation | None:
    """似乎和前端有关"""
    bests = preferred_orientations(mol, rings, fusion_edges)
    return bests[0] if bests else None
