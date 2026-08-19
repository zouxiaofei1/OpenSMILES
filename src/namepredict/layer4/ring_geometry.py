"""3-8 元环标准形状坐标模板与平面几何原语(P-25.3.2.3.1 自建模板坐标)。"""
from __future__ import annotations

import math

# 变形/重叠阈值：单原子最大偏差/环边长上限、非共享环重叠面积/较小环面积上限。
DEFORM_MAX = 0.30
OVERLAP_FRAC = 0.05


def regular_polygon(n: int, *, right_edge_vertical: bool = True) -> list[tuple[float, float]]:
    """单位圆内接正 n 边形坐标，一条边竖直。

    顶点 k 取 θ_k = -π/n + 2πk/n；边 (0,1) 竖直在右侧 (x=cos(π/n))，
    ``right_edge_vertical=False`` 时 x→-x 镜像。偶 n 左、右各一条竖直边，
    奇 n 仅一侧竖直边另一侧为尖点（P-25.3.2.3.1 奇数环 2 模板）。
    """
    pts = [(math.cos(-math.pi / n + 2 * math.pi * k / n),
            math.sin(-math.pi / n + 2 * math.pi * k / n)) for k in range(n)]
    if not right_edge_vertical:
        pts = [(-x, y) for x, y in pts]
    return pts


RING_TEMPLATES: dict[int, list[tuple[float, float]]] = {
    n: regular_polygon(n) for n in range(3, 9)
}


def ring_cyclic(ring_tuple: tuple[int, ...], a: int, b: int) -> list[int]:
    """由 RDKit AtomRings 环元组求以 (a,b) 为第一条边（a→b）的环序。

    b 为 a 在元组中的后继则正序，否则反序；返回以 a 开头、b 第二的顺序。
    """
    n = len(ring_tuple)
    i = ring_tuple.index(a)
    if ring_tuple[(i + 1) % n] == b:
        return [ring_tuple[(i + k) % n] for k in range(n)]
    return [ring_tuple[(i - k) % n] for k in range(n)]


def translate(pts: list[tuple[float, float]], dx: float, dy: float) -> list[tuple[float, float]]:
    """整体平移坐标。"""
    return [(x + dx, y + dy) for x, y in pts]


def reflect_y(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """y 轴对称翻转（x → -x）。"""
    return [(-x, y) for x, y in pts]


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
    """Procrustes 最优刚体+缩放: src→dst，返回 (scale, cos, sin, tx, ty)。

    用于把环实际坐标贴合到正则模板并计算偏差打分（Kabsch 2D）。
    """
    n = len(src)
    cx_s, cy_s = centroid(src)
    cx_d, cy_d = centroid(dst)
    sx = [(x - cx_s, y - cy_s) for x, y in src]
    dx = [(x - cx_d, y - cy_d) for x, y in dst]
    h11 = sum(dx[i][0] * sx[i][0] + dx[i][1] * sx[i][1] for i in range(n))
    h12 = sum(dx[i][0] * sx[i][1] - dx[i][1] * sx[i][0] for i in range(n))
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
