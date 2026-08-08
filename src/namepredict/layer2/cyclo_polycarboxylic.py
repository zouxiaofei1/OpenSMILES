"""Supported neutral cycloalkane di- and tricarboxylic acid parents."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import _cooh_oxygen_idxs
from namepredict.layer2.scaffold.cyclo_pick import _cyclo_parent_candidates
from namepredict.layer2.scaffold.ring_parent import _hetero_or_ring_halo, _outside_carbons, _ring_halo_n, _ring_side_starts
from namepredict.layer3.side_alkyl import _walk_linear
from namepredict.layer2.scaffold.specs import (
    cycloalkane_polyacid_stem, get_spec, numbering_scaffold_facts,
)

def _attach(mol: Mol, ring: set[int], acid: dict) -> int | None:
    hits = [n.GetIdx() for n in mol.GetAtomWithIdx(acid["c_idx"]).GetNeighbors() if n.GetIdx() in ring]
    return hits[0] if len(hits) == 1 else None

def _acid_ring(info: dict) -> set[int] | None:
    acids = info.get("carboxyls") or []
    if len(acids) not in (2, 3): return None
    for ring in _cyclo_parent_candidates(info):
        if all(_attach(info["mol"], ring, acid) is not None for acid in acids): return ring
    return None

def _alkyl_ok(mol: Mol, ring: set[int], acids: list[dict]) -> bool:
    starts = _ring_side_starts(mol, ring, {a["c_idx"] for a in acids})
    if len(starts) > 1: return False
    outside = set(_outside_carbons(mol, ring, {a["c_idx"] for a in acids}))
    if not starts: return not outside
    atoms = _walk_linear(mol, starts[0], ring)
    return atoms is not None and len(atoms) <= 4 and set(atoms) == outside

def _relative_stereo(info: dict, ring: set[int], acids: list[dict]):
    from namepredict.layer1.ring_relative_stereo import ring_relative_stereo
    mol = info["mol"]; order = next((list(r) for r in mol.GetRingInfo().AtomRings() if set(r) == ring), None)
    if order is None: return None
    acid_ids, ligands = {a["c_idx"] for a in acids}, {}
    for atom_id in ring:
        external = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_id).GetNeighbors() if n.GetIdx() not in ring and n.GetAtomicNum() != 1]
        if external: ligands[atom_id] = next((x for x in external if x in acid_ids), external[0])
    return ring_relative_stereo(mol, order, ligands)

def _is_cycloalkane_polycarboxylic(info: dict) -> bool:
    ring, acids = _acid_ring(info), info.get("carboxyls") or []
    if ring is None or any(a.get("anion") for a in acids): return False
    stereo = _relative_stereo(info, ring, acids)
    if stereo is None or stereo.invalid: return False
    allowed = set().union(*(_cooh_oxygen_idxs(info["mol"], a["c_idx"]) for a in acids))
    if not _hetero_or_ring_halo(info["mol"], ring, allowed): return False
    return _ring_halo_n(info["mol"], ring) <= 1 and _alkyl_ok(info["mol"], ring, acids)



def _polyacid_ring_like(info: dict) -> bool:
    acids, mol = info.get("carboxyls") or [], info["mol"]
    return any(sum(_attach(mol, ring, acid) is not None for acid in acids) >= 2 for ring in _cyclo_parent_candidates(info))

def cycloalkane_polycarboxylic_gate(info: dict):
    from namepredict.layer2.candidate_gate import GateScope, pass_gate, scoped_reject
    if not _polyacid_ring_like(info): return pass_gate(GateScope.CYCLOALKANE_POLYCARBOXYLIC)
    if _is_cycloalkane_polycarboxylic(info): return pass_gate(GateScope.CYCLOALKANE_POLYCARBOXYLIC)
    return scoped_reject(GateScope.CYCLOALKANE_POLYCARBOXYLIC, "unsupported_cycloalkane_polyacid")
