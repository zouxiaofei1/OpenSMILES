# IUPAC: P-25.3.2.3.1
# Layer: L4
"""ring_geometry: 3-8 元环模板坐标与平面几何原语。"""
from __future__ import annotations

import math

from namepredict.layer4.ring_geometry import (
    RING_TEMPLATES,
    clip_polygon,
    overlap_area,
    polygon_area,
    regular_polygon,
    ring_cyclic,
    ring_shape_template,
)


def _edge_len(p1, p2):
    return math.hypot(p2[0] - p1[0], p2[1] - p1[1])


def test_regular_polygon_edges_equal_and_vertical():
    for n in range(3, 9):
        pts = RING_TEMPLATES[n]
        assert len(pts) == n
        lens = [_edge_len(pts[i], pts[(i + 1) % n]) for i in range(n)]
        assert max(lens) - min(lens) < 1e-9, n
        # 边(0,1) 竖直（右侧，x 相等）
        assert abs(pts[0][0] - pts[1][0]) < 1e-9
        assert pts[0][0] > 0


def test_regular_polygon_odd_mirror_has_two_variants():
    """奇环模板镜像后仍为正多边形且竖直边在左侧。"""
    for n in (3, 5, 7):
        pts = regular_polygon(n, right_edge_vertical=False)
        assert abs(pts[0][0] - pts[1][0]) < 1e-9
        assert pts[0][0] < 0  # 竖直边在左侧


def test_ring_shape_template_distorted_pentagon():
    """5 元行中双侧融合的变形五元环: 两共享边平行(铺行后竖直)、5 顶点、边01 单位对齐。"""
    order = [0, 1, 2, 3, 4]  # 左共享边(0,1), 右共享边(2,3)
    tpl = ring_shape_template(order, 2)
    assert tpl is not None
    assert len(tpl) == 5
    assert tpl[0] == (0.0, 0.0) and tpl[1] == (1.0, 0.0)  # 边01 为单位对齐边
    # 右共享边(2,3) 与左共享边(0,1) 平行(都沿 x) → 铺行后都竖直
    def _parallel(p0, p1, q0, q1):
        return abs((p1[1] - p0[1]) * (q1[0] - q0[0]) - (p1[0] - p0[0]) * (q1[1] - q0[1])) < 1e-9
    assert _parallel(tpl[0], tpl[1], tpl[2], tpl[3])
    # 非双边(偶环/奇数环端部)返回 None, 调用方保持原弃行逻辑
    assert ring_shape_template(list(range(6)), 2) is None


def test_ring_cyclic_forward_and_reverse():
    ring = (0, 1, 2, 3, 4, 5)  # 六元环
    # b 是 a 的后继 → 正序
    assert ring_cyclic(ring, 0, 1) == [0, 1, 2, 3, 4, 5]
    # b 是 a 的前驱 → 反序
    assert ring_cyclic(ring, 0, 5) == [0, 5, 4, 3, 2, 1]


def test_clip_polygon_disjoint_returns_empty():
    a = [(0, 0), (1, 0), (1, 1), (0, 1)]
    b = [(3, 3), (4, 3), (4, 4), (3, 4)]
    assert clip_polygon(a, b) == []


def test_clip_polygon_intersect_positive_area():
    a = [(0, 0), (2, 0), (2, 2), (0, 2)]
    b = [(1, 1), (3, 1), (3, 3), (1, 3)]
    inter = clip_polygon(a, b)
    assert polygon_area(inter) > 0


def test_polygon_area_rectangle():
    assert abs(polygon_area([(0, 0), (2, 0), (2, 1), (0, 1)]) - 2.0) < 1e-9


def test_overlap_area_known():
    a = [(0, 0), (2, 0), (2, 2), (0, 2)]
    b = [(1, 1), (3, 1), (3, 3), (1, 3)]
    assert abs(overlap_area(a, b) - 1.0) < 1e-9
    # 重叠面积方向对称
    assert abs(overlap_area(b, a) - 1.0) < 1e-9
