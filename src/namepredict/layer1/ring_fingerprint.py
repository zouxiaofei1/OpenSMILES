"""布局级环系指纹（环大小 + 稠合 + 杂原子 + 芳香性）。"""
from __future__ import annotations

from namepredict.layer1.ring_ir import RingComponent, RingSystemIR


def _comp_key(c: RingComponent) -> tuple:
    """返回环组件的排序键（大小、杂原子序数、芳香性）。"""
    zs = tuple(sorted(z for _, z in c.hetero))
    return (c.size, zs, int(c.aromatic))


def _sorted_comps(ir: RingSystemIR) -> tuple[RingComponent, ...]:
    """按排序键返回环系中排好序的组件元组。"""
    return tuple(sorted(ir.components, key=_comp_key))


def _sizes_token(comps: tuple[RingComponent, ...]) -> str:
    """生成环大小序列的指纹片段。"""
    return ",".join(str(c.size) for c in comps)


def _fusion_token(ir: RingSystemIR) -> str:
    """生成稠合信息的指纹片段。"""
    if not ir.fusions:
        return "0"
    shared = sorted(len(f.shared) for f in ir.fusions)
    bonds = sum(1 for f in ir.fusions if f.bond)
    return f"{len(ir.fusions)}:{','.join(map(str, shared))}:b{bonds}"


def _hetero_part(c: RingComponent) -> str:
    """生成单个组件杂原子序数串。"""
    zs = ",".join(str(z) for z in sorted(z for _, z in c.hetero))
    return f"{c.size}:{zs}"


def _hetero_token(comps: tuple[RingComponent, ...]) -> str:
    """生成全部组件杂原子部分的指纹片段。"""
    return ";".join(_hetero_part(c) for c in comps)


def _arom_token(comps: tuple[RingComponent, ...]) -> str:
    """生成全部组件芳香性标志串。"""
    return "".join("1" if c.aromatic else "0" for c in comps)


def ring_fingerprint(ir: RingSystemIR) -> str:
    """布局指纹：topology|sizes|fusion|hetero|aromatic。"""
    comps = _sorted_comps(ir)
    parts = (
        ir.topology,
        _sizes_token(comps),
        _fusion_token(ir),
        _hetero_token(comps),
        _arom_token(comps),
    )
    return "|".join(parts)
