"""Layout-level ring-system fingerprint (sizes + fusion + hetero + aromatic)."""
from __future__ import annotations

from namepredict.layer1.ring_ir import RingComponent, RingSystemIR


def _comp_key(c: RingComponent) -> tuple:
    zs = tuple(sorted(z for _, z in c.hetero))
    return (c.size, zs, int(c.aromatic))


def _sorted_comps(ir: RingSystemIR) -> tuple[RingComponent, ...]:
    return tuple(sorted(ir.components, key=_comp_key))


def _sizes_token(comps: tuple[RingComponent, ...]) -> str:
    return ",".join(str(c.size) for c in comps)


def _fusion_token(ir: RingSystemIR) -> str:
    if not ir.fusions:
        return "0"
    shared = sorted(len(f.shared) for f in ir.fusions)
    bonds = sum(1 for f in ir.fusions if f.bond)
    return f"{len(ir.fusions)}:{','.join(map(str, shared))}:b{bonds}"


def _hetero_part(c: RingComponent) -> str:
    zs = ",".join(str(z) for z in sorted(z for _, z in c.hetero))
    return f"{c.size}:{zs}"


def _hetero_token(comps: tuple[RingComponent, ...]) -> str:
    return ";".join(_hetero_part(c) for c in comps)


def _arom_token(comps: tuple[RingComponent, ...]) -> str:
    return "".join("1" if c.aromatic else "0" for c in comps)


def ring_fingerprint(ir: RingSystemIR) -> str:
    """Layout fingerprint: topology|sizes|fusion|hetero|aromatic."""
    comps = _sorted_comps(ir)
    parts = (
        ir.topology,
        _sizes_token(comps),
        _fusion_token(ir),
        _hetero_token(comps),
        _arom_token(comps),
    )
    return "|".join(parts)
