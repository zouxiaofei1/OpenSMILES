"""Generalized arene FG parent: OH/NH2 on any aromatic fused/monocyclic ring.

Data-driven: _FG_SPEC maps FG type to suffix/rank; _SCAFFOLD_DETECT maps
scaffold_id to detection functions.  _try_arene_fg_parent is the ONE
generalized entry point.

IUPAC P-25 numbering via scaffold chain.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _outside_ok, _ring_halo_n, _ring_side_starts


def _detect_naphthalene(info: dict) -> tuple | None:
    """Return (ring_atoms, chain) or None."""
    from namepredict.layer2.naphthalene import _is_naphthalene_core, _naph_chains
    if not _is_naphthalene_core(info):
        return None
    rings = info.get("rings") or []
    r1, r2 = list(rings[0]["atom_ids"]), list(rings[1]["atom_ids"])
    chains = _naph_chains(info)
    return (set(r1) | set(r2), chains[0] if chains else [])


def _detect_quinoline(info: dict) -> tuple | None:
    """Return (ring_atoms, chain) or None."""
    from namepredict.layer2.quinoline import _q_core, _build_chain, _ring_set
    parts = _q_core(info)
    if parts is None:
        return None
    return _ring_set(parts), _build_chain(info, parts)


_FG_BLOCK_KEYS = (
    "has_acid", "has_aldehyde", "has_ketone", "has_ester",
    "has_amide", "has_nitrile", "has_acyl_chloride", "has_anhydride",
    "has_nitro", "has_thiol",
)

_FG_SPEC: dict[str, dict] = {
    "OH": {
        "ekey": "hydroxyls", "ckey": "c_idx",
        "fk": "oh_c_idxs", "rank": 5,
        "suffix_en": {1: "ol", 2: "diol"},
        "suffix_zh": {1: "醇", 2: "二酚"},
        "kind_suffix": {1: "ol", 2: "diol"},
        "extra_block": ("has_amine",),
    },
    "NH2": {
        "ekey": "amines", "ckey": "c_idx",
        "fk": "amine_c_idxs", "rank": 3,
        "suffix_en": {1: "amine", 2: "diamine"},
        "suffix_zh": {1: "胺", 2: "二胺"},
        "kind_suffix": {1: "amine", 2: "diamine"},
        "extra_block": ("has_alcohol",),
    },
}

_SCAFFOLD_DETECT = {
    "naphthalene": (_detect_naphthalene, "naphthalene"),
    "quinoline": (_detect_quinoline, "quinoline"),
}

_ARENE_FG_KINDS: dict[tuple, str] = {
    ("naphthalene", "OH", 1): "naphthalenol",
    ("naphthalene", "OH", 2): "naphthalenediol",
    ("naphthalene", "NH2", 1): "naphthalenamine",
    ("quinoline", "OH", 2): "quinolinediol",
}


def _arene_fg_conflict(info: dict, fg_type: str) -> bool:
    spec = _FG_SPEC.get(fg_type) or {}
    extra = spec.get("extra_block") or ()
    return any(info.get(k) for k in _FG_BLOCK_KEYS + extra)


def _fg_allowed_atoms(on_ring: list[dict], spec: dict) -> set[int]:
    ckey = spec["ckey"]
    if spec.get("ekey") == "hydroxyls":
        return {e["o_idx"] for e in on_ring}
    if spec.get("ekey") == "amines":
        return {e["n_idx"] for e in on_ring}
    return set()


def _arene_fg_subs_ok(
    info: dict, mol: Mol, ring: set[int], fg_idxs: set[int], allowed: set[int],
) -> bool:
    """Gate: only halo + simple alkyl side chains; no extra heteroatoms."""
    if not _outside_ok(mol, ring, allowed):
        return False
    starts = _ring_side_starts(mol, ring, fg_idxs)
    return _ring_halo_n(mol, ring) + len(starts) <= 2


def _resolve_arene_core(info, scaffold_id, fg_type):
    """(ring_atoms, chain, spec) or None."""
    spec = _FG_SPEC.get(fg_type)
    if spec is None: return None
    entry = _SCAFFOLD_DETECT.get(scaffold_id)
    if entry is None: return None
    core = entry[0](info)
    if core is None: return None
    ring_atoms, chain = core
    return (ring_atoms, chain, spec) if chain else None


def _check_fg_count(info, spec, ring_atoms, n_fg):
    """Return validated FG entries on ring, or None."""
    fg_entries = info.get(spec["ekey"]) or []
    on_ring = [e for e in fg_entries if e[spec["ckey"]] in ring_atoms]
    if len(on_ring) != n_fg or len(on_ring) != len(fg_entries): return None
    return on_ring


def _gate_scaffold_fg(info, fg_type, ring_atoms, fg_idxs, spec):
    """True if no conflicting FGs and ring substituents are acceptable."""
    if _arene_fg_conflict(info, fg_type): return False
    on_ring = [e for e in (info.get(spec["ekey"]) or []) if e[spec["ckey"]] in ring_atoms]
    allowed = _fg_allowed_atoms(on_ring, spec)
    return _arene_fg_subs_ok(info, info["mol"], ring_atoms, fg_idxs, allowed)


def _assemble_fg_parent(chain, kind, spec, fg_idxs, n_fg):
    """Build parent dict with FG locant keys."""
    parent = {"chain": chain, "n_carbons": len(chain), "kind": kind, "scaffold_id": kind}
    if n_fg == 1:
        parent[spec["fk"].replace("_idxs", "_idx")] = fg_idxs[0]
    else:
        parent[spec["fk"]] = fg_idxs
    return parent


def _try_arene_fg_parent(info, scaffold_id, fg_type, n_fg=1):
    """Generalized: OH/NH2 on aromatic ring -> systematic parent dict."""
    if (resolved := _resolve_arene_core(info, scaffold_id, fg_type)) is None: return None
    ring_atoms, chain, spec = resolved
    if (on_ring := _check_fg_count(info, spec, ring_atoms, n_fg)) is None: return None
    fg_idxs_set = {e[spec["ckey"]] for e in on_ring}
    if not _gate_scaffold_fg(info, fg_type, ring_atoms, fg_idxs_set, spec): return None
    fg_idxs_list = [e[spec["ckey"]] for e in on_ring]
    kind = _ARENE_FG_KINDS.get((scaffold_id, fg_type, n_fg))
    return _assemble_fg_parent(chain, kind, spec, fg_idxs_list, n_fg) if kind else None


# ---- Thin wrappers for per-scaffold FG parents ----

def _try_naphthalenol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "OH")


def _try_naphthalenediol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "OH", n_fg=2)


def _try_naphthalenamine_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "NH2")


def _try_quinolinediol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "quinoline", "OH", n_fg=2)
