"""Retained scaffold registry (P-22 / P-25) — match metadata for L2.

Phase-0 skeleton: entries describe topology predicates using ring_systems.
Existing _try_* modules remain authoritative until migrated.
"""
from __future__ import annotations

from typing import Callable

from rdkit.Chem import Mol

# scaffold_id → entry
# entry keys: kind, n_rings, n_atoms, hetero_Z (sorted), topology, en, zh
RetainedEntry = dict

_REGISTRY: dict[str, RetainedEntry] = {
    "benzene": {
        "kind": "benzene",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (),
        "topology": "mono",
        "aromatic": True,
        "en": "benzene",
        "zh": "苯",
    },
    "pyridine": {
        "kind": "pyridine",
        "n_rings": 1,
        "n_atoms": 6,
        "hetero_Z": (7,),
        "topology": "mono",
        "aromatic": True,
        "en": "pyridine",
        "zh": "吡啶",
    },
    "naphthalene": {
        "kind": "naphthalene",
        "n_rings": 2,
        "n_atoms": 10,
        "hetero_Z": (),
        "topology": "fused",
        "aromatic": True,
        "en": "naphthalene",
        "zh": "萘",
    },
    "indole": {
        "kind": "indole",
        "n_rings": 2,
        "n_atoms": 9,
        "hetero_Z": (7,),
        "topology": "fused",
        "aromatic": True,
        "en": "1H-indole",
        "zh": "1H-吲哚",
    },
}


def registry() -> dict[str, RetainedEntry]:
    return _REGISTRY


def get_entry(scaffold_id: str) -> RetainedEntry | None:
    return _REGISTRY.get(scaffold_id)


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
        for sid, entry in _REGISTRY.items():
            if _system_matches(system, entry):
                out.append((sid, system, entry))
    return out


def match_scaffold_ids(info: dict) -> list[str]:
    return [sid for sid, _, _ in match_systems(info)]
