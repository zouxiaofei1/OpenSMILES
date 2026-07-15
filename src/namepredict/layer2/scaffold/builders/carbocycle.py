"""Mono carbocycle scaffold builder: alkane / ene / polyene by endocyclic C=C."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem import Mol

from namepredict.layer1.ring_ir import RingSystemIR, build_ring_ir
from namepredict.layer2.ring_parent import _ring_double_count


@dataclass(frozen=True)
class ScaffoldHit:
    spec_id: str
    ring_size: int
    n_double: int
    atom_ids: tuple[int, ...]


def _as_mol(src) -> Mol | None:
    if isinstance(src, Mol):
        return src
    if isinstance(src, dict) and isinstance(src.get("mol"), Mol):
        return src["mol"]
    return None


def _mono_systems(mol: Mol) -> list[RingSystemIR]:
    return [s for s in build_ring_ir(mol) if s.topology == "mono"]


def _all_carbon(mol: Mol, atom_ids: tuple[int, ...]) -> bool:
    return all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in atom_ids)


def _sssr_order(mol: Mol, atom_ids: tuple[int, ...]) -> tuple[int, ...] | None:
    want = frozenset(atom_ids)
    for ring in mol.GetRingInfo().AtomRings():
        if frozenset(ring) == want:
            return tuple(ring)
    return None


def _is_arom_c6(mol: Mol, atom_ids: tuple[int, ...]) -> bool:
    if len(atom_ids) != 6:
        return False
    return all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids)


def _spec_for(n_double: int) -> str | None:
    if n_double == 0:
        return "cycloalkane"
    if n_double == 1:
        return "cycloalkene"
    if n_double >= 2:
        return "cyclopolyene"
    return None


def _carbo_ok(mol: Mol, ir: RingSystemIR) -> bool:
    if len(ir.components) != 1 or ir.components[0].hetero:
        return False
    return _all_carbon(mol, ir.atom_ids) and not _is_arom_c6(mol, ir.atom_ids)


def _hit_from(mol: Mol, ir: RingSystemIR) -> ScaffoldHit | None:
    if not _carbo_ok(mol, ir):
        return None
    order = _sssr_order(mol, ir.atom_ids)
    if order is None:
        return None
    n = _ring_double_count(mol, order)
    sid = _spec_for(n)
    return None if sid is None else ScaffoldHit(sid, len(order), n, order)


def try_carbocycle(src) -> ScaffoldHit | None:
    """Identify mono all-C carbocycle; None for arene / open chain / multi-ring."""
    mol = _as_mol(src)
    if mol is None:
        return None
    mono = _mono_systems(mol)
    if len(mono) != 1:
        return None
    return _hit_from(mol, mono[0])
