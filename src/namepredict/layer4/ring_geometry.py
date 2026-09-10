"""3-8 元环标准形状坐标模板与平面几何原语(P-25.3.2.3.1 自建模板坐标)。"""
from __future__ import annotations

import math

# 变形/重叠阈值：单原子最大偏差/环边长上限、非共享环重叠面积/较小环面积上限。
DEFORM_MAX = 0.30
OVERLAP_FRAC = 0.05


def regular_polygon(n: int, *, right_edge_vertical: bool = True) -> list[tuple[float, float]]:
    """单位圆内接正 n 边形坐标，一条边竖直。 """
    pts = [(math.cos(-math.pi / n + 2 * math.pi * k / n),
            math.sin(-math.pi / n + 2 * math.pi * k / n)) for k in range(n)]
    if not right_edge_vertical:
        pts = [(-x, y) for x, y in pts]
    return pts


RING_TEMPLATES: dict[int, list[tuple[float, float]]] = {
    n: regular_polygon(n) for n in range(3, 20)
}

# 行内中间奇环的"行宽"：铺成水平行时相邻环共享边间的水平间距。
# 行首环为正六边形(边长 1)时对边距 = sqrt(3)。变形五元/七元模板用该行宽，
# 使其左右两条共享竖边与邻环衔接。P-25.3.2.3.2：变形环应尽可能小。
_ROW_WIDTH = 3 ** 0.5


def ring_shape_template(order: list[int], exit_idx: int,
                        row_width: float = _ROW_WIDTH) -> list[tuple[float, float]] | None:
    """P-25.3.2.3.2 变形环模板：让奇数环在水平行中间也能"两侧竖直边"稠合。 """
    n = len(order)
    if n not in (5, 7):
        return None
    e = exit_idx % n
    if n == 5:
        # 顶点落在哪个索引由行方向决定：出口边为 (2,3) 时顶点是末尾索引 4，出口边为 (3,4) 时顶点是索引 2(镜像形)。
        # 两个方向都要给模板，否则整行因摆不出被丢弃，行内最右环也一并丢失(6-5-6 会退化成只能从另一侧起编)。
        if e == 2:
            return [(0.0, 0.0), (1.0, 0.0), (1.0, row_width), (0.0, row_width), (-0.4, 0.0)]
        if e == 3:
            return [(0.0, 0.0), (1.0, 0.0), (1.4, 0.0), (1.0, row_width), (0.0, row_width)]
        return None
    if n == 7 and e not in (2, 3):
        return None
    positions: list[tuple[float, float]] = []  # 逐顶点分配：索引 0,1 为左共享边(单位水平, y=0)；索引 e,e+1 为右共享边(平行水平, y=row_width)
    mid_slot = 0.0  # 已用过的中间剩余顶点计数（n=7 时 2..e-1 会有前向连接顶点）
    for i in range(n):
        if i == 0:
            positions.append((0.0, 0.0))
        elif i == 1:
            positions.append((1.0, 0.0))
        elif i == e:
            positions.append((1.0, row_width))   # 右共享边起点：正对 order[1] 上方
        elif i == e + 1:
            positions.append((0.0, row_width))   # 右共享边终点：正对 order[0] 上方
        elif i == n - 1:
            positions.append((-0.4, 0.0))        # 末尾顶点：左侧折线连接 order[n-1]→order[0]
        else:
            positions.append((-0.4, row_width * 0.4))  # 其余顶点：左侧中部
    return positions


def ring_cyclic(ring_tuple: tuple[int, ...], a: int, b: int) -> list[int]:
    """由 RDKit AtomRings 环元组求以 (a,b) 为首边（a→b 正序/反序）的环序，返回 a 开头、b 第二。"""
    n = len(ring_tuple)
    i = ring_tuple.index(a)
    if ring_tuple[(i + 1) % n] == b:
        return [ring_tuple[(i + k) % n] for k in range(n)]
    return [ring_tuple[(i - k) % n] for k in range(n)]


def centroid(pts: list[tuple[float, float]]) -> tuple[float, float]:
    """点集质心。"""
    n = len(pts)
    return (sum(x for x, _ in pts) / n, sum(y for _, y in pts) / n)


def apply_rigid(p: tuple[float, float], params) -> tuple[float, float]:
    """应用 rigid_fit 返回的 (scale, cos, sin, tx, ty) 变换。"""
    scale, cos_t, sin_t, tx, ty = params
    return (scale * (cos_t * p[0] - sin_t * p[1]) + tx,
            scale * (sin_t * p[0] + cos_t * p[1]) + ty)


def rigid_fit(src: list[tuple[float, float]], dst: list[tuple[float, float]]):
    """Procrustes 最优刚体+缩放 src→dst (Kabsch 2D)，返回 (scale, cos, sin, tx, ty) 供贴合模板与偏差打分。"""
    n = len(src)
    cx_s, cy_s = centroid(src)
    cx_d, cy_d = centroid(dst)
    sx = [(x - cx_s, y - cy_s) for x, y in src]
    dx = [(x - cx_d, y - cy_d) for x, y in dst]
    h11 = sum(dx[i][0] * sx[i][0] + dx[i][1] * sx[i][1] for i in range(n))
    h12 = sum(dx[i][1] * sx[i][0] - dx[i][0] * sx[i][1] for i in range(n))  # Kabsch 最优旋转角 θ = atan2(Σ d_y s_x - d_x s_y, Σ d·s)，让 Rθ·s ≈ d
    theta = math.atan2(h12, h11)
    s_sq = sum(x * x + y * y for x, y in sx)
    d_sq = sum(x * x + y * y for x, y in dx)
    scale = math.sqrt(d_sq / s_sq) if s_sq > 1e-12 else 1.0
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    tx = cx_d - scale * (cos_t * cx_s - sin_t * cy_s)
    ty = cy_d - scale * (sin_t * cx_s + cos_t * cy_s)
    return scale, cos_t, sin_t, tx, ty


def _inside(p: tuple[float, float], a: tuple[float, float], b: tuple[float, float]) -> bool:
    """点 p 是否在定向边 a→b 的左侧（内侧）。"""
    return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) >= 0


def _line_inter(p: tuple[float, float], q: tuple[float, float],
                a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
    """线段 pq 与直线 ab 的交点。"""
    d_pq = (q[0] - p[0], q[1] - p[1])
    d_ab = (b[0] - a[0], b[1] - a[1])
    denom = d_pq[0] * d_ab[1] - d_pq[1] * d_ab[0]
    if abs(denom) < 1e-12:
        return q
    t = ((a[0] - p[0]) * d_ab[1] - (a[1] - p[1]) * d_ab[0]) / denom
    return (p[0] + t * d_pq[0], p[1] + t * d_pq[1])


def clip_polygon(subject: list[tuple[float, float]],
                 clip: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Sutherland–Hodgman 多边形裁剪，返回 subject∩clip 顶点（可能空）。"""
    out = list(subject)
    for i in range(len(clip)):
        a, b = clip[i], clip[(i + 1) % len(clip)]
        if not out:
            return []
        inp, out = out, []
        s = inp[-1]
        for e in inp:
            if _inside(e, a, b):
                if not _inside(s, a, b):
                    out.append(_line_inter(s, e, a, b))
                out.append(e)
            elif _inside(s, a, b):
                out.append(_line_inter(s, e, a, b))
            s = e
    return out


def polygon_area(poly: list[tuple[float, float]]) -> float:
    """鞋带公式求多边形面积（绝对值）。"""
    if len(poly) < 3:
        return 0.0
    return abs(sum(poly[i][0] * poly[(i + 1) % len(poly)][1]
                   - poly[(i + 1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))) / 2.0


def overlap_area(p1: list[tuple[float, float]], p2: list[tuple[float, float]]) -> float:
    """两多边形重叠面积（用 p1 裁剪 p2 后求面积，方向对称）。"""
    return polygon_area(clip_polygon(p1, p2))
