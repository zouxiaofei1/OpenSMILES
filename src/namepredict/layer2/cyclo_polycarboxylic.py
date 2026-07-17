"""Supported neutral cycloalkane di- and tricarboxylic acid parents."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.arene_carbonyl import _cooh_oxygen_idxs
from namepredict.layer2.cyclo_pick import _cyclo_parent_candidates
from namepredict.layer2.ring_parent import _hetero_or_ring_halo, _outside_carbons, _ring_halo_n, _ring_side_starts
from namepredict.layer2.side_alkyl import _walk_linear
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


def _has_stereo(mol: Mol) -> bool:
    return any(a.GetChiralTag().name != "CHI_UNSPECIFIED" for a in mol.GetAtoms())


def _is_cycloalkane_polycarboxylic(info: dict) -> bool:
    ring, acids = _acid_ring(info), info.get("carboxyls") or []
    if ring is None or _has_stereo(info["mol"]) or any(a.get("anion") for a in acids): return False
    allowed = set().union(*(_cooh_oxygen_idxs(info["mol"], a["c_idx"]) for a in acids))
    if not _hetero_or_ring_halo(info["mol"], ring, allowed): return False
    return _ring_halo_n(info["mol"], ring) <= 1 and _alkyl_ok(info["mol"], ring, acids)


def _cycloalkane_polycarboxylic_parent(info: dict) -> dict:
    ring, spec = _acid_ring(info) or set(), get_spec("cycloalkane_polycarboxylic")
    stems = cycloalkane_polyacid_stem(len(ring)) or (spec.stem_en, spec.stem_zh)
    attaches = [_attach(info["mol"], ring, acid) for acid in info["carboxyls"]]
    return {
        "chain": list(ring), "n_carbons": len(ring),
        "kind": "cycloalkane_polycarboxylic", "scaffold_id": spec.id,
        "stem_en": stems[0], "stem_zh": stems[1], "cooh_c_idxs": attaches,
        "acid_count": len(attaches),
        "numbering_scaffold": numbering_scaffold_facts(spec.id, len(ring)),
    }


def _try_cycloalkane_polycarboxylic(info: dict) -> dict | None:
    return _cycloalkane_polycarboxylic_parent(info) if _is_cycloalkane_polycarboxylic(info) else None


def _polyacid_ring_like(info: dict) -> bool:
    acids, mol = info.get("carboxyls") or [], info["mol"]
    return any(sum(_attach(mol, ring, acid) is not None for acid in acids) >= 2 for ring in _cyclo_parent_candidates(info))


def cycloalkane_polycarboxylic_gate(info: dict):
    from namepredict.layer2.candidate_gate import GateScope, pass_gate, scoped_reject
    if not _polyacid_ring_like(info): return pass_gate(GateScope.CYCLOALKANE_POLYCARBOXYLIC)
    if _is_cycloalkane_polycarboxylic(info): return pass_gate(GateScope.CYCLOALKANE_POLYCARBOXYLIC)
    return scoped_reject(GateScope.CYCLOALKANE_POLYCARBOXYLIC, "unsupported_cycloalkane_polyacid")
