"""Retained scaffold registry (P-22 / P-25) — topology match; stems from Spec.

Phase-0 topology predicates using ring_systems. Stem en/zh come from
scaffold.specs (single source). Existing _try_* modules remain authoritative
until fully migrated.
"""
from __future__ import annotations

from namepredict.layer2.specs import get_spec

# scaffold_id → topology entry (en/zh resolved via get_spec at read time)
# topology keys: kind, n_rings, n_atoms, hetero_Z (sorted), topology, aromatic
RetainedEntry = dict

# Topology-only table; ids MUST be ⊆ ScaffoldSpec registry.
_TOPOLOGY: dict[str, dict] = {
    "benzene": {
        "kind": "benzene",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (),
        "topology": "mono",
        "aromatic": True,
    },
    "pyridine": {
        "kind": "pyridine",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (7,),
        "topology": "mono",
        "aromatic": True,
    },
    "naphthalene": {
        "kind": "naphthalene",
        "n_rings": 2,
        "n_atoms": 10,
        "hetero_Z": (),
        "topology": "fused",
        "aromatic": True,
    },
    "indole": {
        "kind": "indole",
        "n_rings": 2,
        "n_atoms": 9,
        "hetero_Z": (7,),
        "topology": "fused",
        "aromatic": True,
    },
}


def _entry_with_stems(sid: str, topo: dict) -> RetainedEntry:
    sp = get_spec(sid)
    en = sp.stem_en if sp else None
    zh = sp.stem_zh if sp else None
    return {**topo, "en": en, "zh": zh}


def registry() -> dict[str, RetainedEntry]:
    return {sid: _entry_with_stems(sid, t) for sid, t in _TOPOLOGY.items()}


def get_entry(scaffold_id: str) -> RetainedEntry | None:
    topo = _TOPOLOGY.get(scaffold_id)
    return None if topo is None else _entry_with_stems(scaffold_id, topo)


def _hetero_Z_tuple(system: dict) -> tuple[int, ...]:
    return tuple(sorted(h["Z"] for h in system.get("hetero_atoms") or []))


def _system_matches(system: dict, entry: RetainedEntry) -> bool:
    if system.get("n_rings") != entry["n_rings"]:
        return False
    if system.get("n_atoms") != entry["n_atoms"]:
        return False
    if system.get("topology") != entry["topology"]:
        return False
    if bool(system.get("is_aromatic_mancude")) != bool(entry.get("aromatic")):
        return False
    return _hetero_Z_tuple(system) == tuple(entry["hetero_Z"])


def match_systems(info: dict) -> list[tuple[str, dict, RetainedEntry]]:
    """Return (scaffold_id, ring_system, entry) for matching systems."""
    out: list[tuple[str, dict, RetainedEntry]] = []
    for system in info.get("ring_systems") or []:
        for sid, topo in _TOPOLOGY.items():
            entry = _entry_with_stems(sid, topo)
            if _system_matches(system, entry):
                out.append((sid, system, entry))
    return out


def match_scaffold_ids(info: dict) -> list[str]:
    return [sid for sid, _, _ in match_systems(info)]
