"""由 SSSR 稠合拓扑构建的带类型 RingSystemIR（仅 L1 事实）。"""
from __future__ import annotations

from dataclasses import dataclass, replace

from rdkit.Chem import Mol

from namepredict.layer1.ring_systems import (
    _is_arom_mancude,
    _sssr,
    build_ring_systems,
)


@dataclass(frozen=True)
class RingComponent:
    sssr_idx: int
    size: int
    atom_ids: tuple[int, ...]
    hetero: tuple[tuple[int, int], ...]
    aromatic: bool


@dataclass(frozen=True)
class FusionEdge:
    a: int
    b: int
    shared: tuple[int, ...]
    bond: bool


@dataclass(frozen=True)
class RingSystemIR:
    atom_ids: tuple[int, ...]
    components: tuple[RingComponent, ...]
    fusions: tuple[FusionEdge, ...]
    topology: str
    fingerprint: str | None = ""


def _ring_hetero(mol: Mol, atom_ids: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    out: list[tuple[int, int]] = []
    for i in atom_ids:
        z = mol.GetAtomWithIdx(i).GetAtomicNum()
        if z != 6:
            out.append((i, z))
    return tuple(out)


def _component(mol: Mol, sssr_idx: int, ring: tuple[int, ...]) -> RingComponent:
    atoms = tuple(sorted(ring))
    return RingComponent(
        sssr_idx=sssr_idx,
        size=len(atoms),
        atom_ids=atoms,
        hetero=_ring_hetero(mol, atoms),
        aromatic=_is_arom_mancude(mol, atoms),
    )


def _shared_bond(mol: Mol, shared: tuple[int, ...]) -> bool:
    if len(shared) < 2:
        return False
    a, b = shared[0], shared[1]
    return mol.GetBondBetweenAtoms(a, b) is not None


def _fusion(mol: Mol, edge: tuple) -> FusionEdge:
    a, b, shared_list = edge
    shared = tuple(shared_list)
    return FusionEdge(a=a, b=b, shared=shared, bond=_shared_bond(mol, shared))


def _components(
    mol: Mol, rings: list[tuple[int, ...]], indices: list[int],
) -> tuple[RingComponent, ...]:
    return tuple(_component(mol, i, rings[i]) for i in sorted(indices))


def _fusions(mol: Mol, edges: list) -> tuple[FusionEdge, ...]:
    return tuple(_fusion(mol, e) for e in edges)


def _system_ir(
    mol: Mol, rings: list[tuple[int, ...]], system: dict,
) -> RingSystemIR:
    from namepredict.layer1.ring_fingerprint import ring_fingerprint

    ir = RingSystemIR(
        atom_ids=tuple(system["atom_ids"]),
        components=_components(mol, rings, system["sssr_indices"]),
        fusions=_fusions(mol, system["fusion_edges"]),
        topology=system["topology"],
        fingerprint="",
    )
    return replace(ir, fingerprint=ring_fingerprint(ir))


def build_ring_ir(mol: Mol) -> list[RingSystemIR]:
    """由 SSSR 稠合图构建带类型的环系 IR。"""
    rings = _sssr(mol)
    systems = build_ring_systems(mol)
    return [_system_ir(mol, rings, s) for s in systems]
