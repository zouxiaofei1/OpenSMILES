# IUPAC: P-25.3.2.3
# Layer: L4
"""fused_orientation: 优选取向(水平行/坐标/象限加权)。"""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.ring_systems import build_ring_systems
from namepredict.layer4.fused_orientation import Orientation, preferred_orientation


def _orient(smiles: str) -> Orientation | None:
    mol = preprocess(smiles)
    assert mol is not None
    systems = build_ring_systems(mol)
    if not systems:
        return None
    rings = list(mol.GetRingInfo().AtomRings())
    return preferred_orientation(mol, rings, systems[0]["fusion_edges"])


def _shared_vertical(orient: Orientation, rings, row) -> bool:
    """相邻行环共享边是否竖直。"""
    cd = orient.coord_dict()
    for k in range(len(row) - 1):
        shared = sorted(set(rings[row[k]]) & set(rings[row[k + 1]]))
        x0, _ = cd[shared[0]]
        x1, _ = cd[shared[1]]
        if abs(x0 - x1) > 1e-6:
            return False
    return True


def test_naphthalene_two_ring_row_symmetric():
    o = _orient("c1ccc2ccccc2c1")
    assert o is not None
    assert len(o.row) == 2
    q1, _, q3, _ = o.quad
    assert abs(q1 - q3) < 0.01  # 对称碳环上下镜像等价


def test_anthracene_three_ring_row():
    """anthracene 线性 3 苯环 → 水平行 3 环。"""
    mol = preprocess("c1ccc2cc3ccccc3cc2c1")
    rings = list(mol.GetRingInfo().AtomRings())
    o = _orient("c1ccc2cc3ccccc3cc2c1")
    assert o is not None and len(o.row) == 3
    assert _shared_vertical(o, rings, o.row)


def test_phenanthrene_angular_two_ring_row():
    """phenanthrene 角型 3 环 → 水平行最多 2 环，右上象限环数最多。"""
    mol = preprocess("c1ccc2c(c1)ccc1ccccc12")
    rings = list(mol.GetRingInfo().AtomRings())
    o = _orient("c1ccc2c(c1)ccc1ccccc12")
    assert o is not None and len(o.row) == 2
    q1, _, q3, _ = o.quad
    assert q1 > q3  # 右上环数优先于左下（phenanthrene 取向）


def test_indole_and_quinoline_two_ring_row():
    for smi in ("c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1"):
        o = _orient(smi)
        assert o is not None and len(o.row) == 2


def test_azulene_odd_ring_pair_row():
    """azulene 5+7 奇环对：共享边竖直（奇环端部模板）。"""
    mol = preprocess("C1=CC=C2C=CC=CC2=C1")
    rings = list(mol.GetRingInfo().AtomRings())
    o = _orient("C1=CC=C2C=CC=CC2=C1")
    assert o is not None and len(o.row) == 2
    assert _shared_vertical(o, rings, o.row)


def test_pyrene_max_row_two():
    o = _orient("c1cc2ccc3cccc4ccc(c1)c2c34")
    assert o is not None and len(o.row) == 2


def test_6_5_6_linear_three_ring_row():
    """6-5-6 线性稠环(中间 5 元奇环两侧稠合) → P-25.3.2.3.2 变形五元环使三环成水平行。"""
    smi = "ClC1C2NC3C=CC=CC=3C=2N=CN=1"
    mol = preprocess(smi)
    rings = list(mol.GetRingInfo().AtomRings())
    o = _orient(smi)
    assert o is not None
    assert len(o.row) == 3
    assert _shared_vertical(o, rings, o.row)
