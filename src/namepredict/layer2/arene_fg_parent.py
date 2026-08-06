"""Generalized arene FG parent: OH/NH2/CHO/CN on any aromatic fused/monocyclic ring.

Data-driven: _FG_SPEC maps FG type to suffix/rank; _SCAFFOLD_DETECT maps
scaffold_id to detection functions.  _try_arene_fg_parent is the ONE
generalized entry point.

IUPAC P-25 numbering via scaffold chain.
"""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.scaffold.ring_parent import _dbl_o_idx, _outside_ok, _ring_halo_n, _ring_side_starts
from namepredict.layer2.arene_carbonyl import _nitrile_n_idx


# --- Extra atoms for exocyclic FG (CHO carbonyl O, CN nitrile N) ---

def _carbonyl_o(mol: Mol, c_idx: int) -> set[int]:
    o = _dbl_o_idx(mol, c_idx)
    return {o} if o is not None else set()


def _nitrile_n(mol: Mol, c_idx: int) -> set[int]:
    n = _nitrile_n_idx(mol, c_idx)
    return {n} if n is not None else set()


_EXTRA_ATOMS_FNS = {"_carbonyl_o": _carbonyl_o, "_nitrile_n": _nitrile_n}


def _detect_naphthalene(info: dict) -> tuple | None:
    """Return (ring_atoms, chain) or None."""
    from namepredict.layer2.scaffold.naphthalene import _is_naphthalene_core, _naph_chains
    if not _is_naphthalene_core(info):
        return None
    rings = info.get("rings") or []
    r1, r2 = list(rings[0]["atom_ids"]), list(rings[1]["atom_ids"])
    chains = _naph_chains(info)
    return (set(r1) | set(r2), chains[0] if chains else [])


def _detect_quinoline(info: dict) -> tuple | None:
    """Return (ring_atoms, chain) or None."""
    from namepredict.layer2.scaffold.quinoline import _q_core, _build_chain, _ring_set
    parts = _q_core(info)
    if parts is None:
        return None
    return _ring_set(parts), _build_chain(info, parts)


def _detect_pyrazole(info: dict):
    """Return (ring_atoms, chain, meta) for pyrazole-amine parent."""
    from namepredict.layer2.scaffold.heteroarene5 import _pyrazole_n_pair
    pair = _pyrazole_n_pair(info)
    if pair is None: return None
    nh, n = pair
    ring = list(info["rings"][0]["atom_ids"])
    chain = [nh] + [i for i in ring if i != nh]
    return set(ring), chain, {"nh_idx": nh, "n_idx": n}


def _detect_thiazole(info: dict):
    """Return (ring_atoms, chain, meta) for thiazole-amine parent."""
    from namepredict.layer2.scaffold.azole13 import _azole13_hetero_pair
    pair = _azole13_hetero_pair(info)
    if pair is None: return None
    h, n = pair
    ring = list(info["rings"][0]["atom_ids"])
    chain = [h] + [i for i in ring if i != h]
    return set(ring), chain, {"hetero_idx": h, "n_idx": n}


def _detect_quinazoline(info: dict):
    """Return (ring_atoms, chain, meta) for quinazoline-amine parent."""
    from namepredict.layer2.scaffold.benzodiazine import _qz_core, _n1
    core = _qz_core(info)
    if core is None: return None
    atoms, ns, bridge = core
    n1_idx = _n1(info["mol"], ns, bridge)
    chain = [n1_idx] + [i for i in sorted(atoms) if i != n1_idx]
    return atoms, chain, {"n_idxs": ns}


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
    "CHO": {
        "ekey": "aldehydes", "ckey": "c_idx", "fk": "aldehyde_c_idxs", "rank": 7,
        "suffix_en": {1: "carbaldehyde"}, "suffix_zh": {1: "甲醛"},
        "kind_suffix": {1: "carbaldehyde"},
        "extra_atoms_fn": "_carbonyl_o",
        "self_allow": ("has_aldehyde",),
    },
    "CN": {
        "ekey": "nitriles", "ckey": "c_idx", "fk": "nitrile_c_idxs", "rank": 8,
        "suffix_en": {1: "carbonitrile"}, "suffix_zh": {1: "甲腈"},
        "kind_suffix": {1: "carbonitrile"},
        "extra_atoms_fn": "_nitrile_n",
        "self_allow": ("has_nitrile",),
    },
}

_SCAFFOLD_DETECT = {
    "naphthalene": (_detect_naphthalene, "naphthalene"),
    "quinoline": (_detect_quinoline, "quinoline"),
    "pyrazole": (_detect_pyrazole, "pyrazole"),
    "thiazole": (_detect_thiazole, "thiazole"),
    "quinazoline": (_detect_quinazoline, "quinazoline"),
}

_ARENE_FG_KINDS: dict[tuple, str] = {
    ("naphthalene", "OH", 1): "naphthalenol",
    ("naphthalene", "OH", 2): "naphthalenediol",
    ("naphthalene", "NH2", 1): "naphthalenamine",
    ("quinoline", "OH", 2): "quinolinediol",
    ("naphthalene", "CHO", 1): "naphthalenecarbaldehyde",
    ("quinoline", "CHO", 1): "quinolinecarbaldehyde",
    ("naphthalene", "CN", 1): "naphthalenecarbonitrile",
    ("quinoline", "CN", 1): "quinolinecarbonitrile",
    ("pyrazole", "NH2", 1): "pyrazolamine",
    ("thiazole", "NH2", 1): "thiazolamine",
    ("quinazoline", "NH2", 1): "quinazolinamine",
}


def _arene_fg_conflict(info: dict, fg_type: str) -> bool:
    spec = _FG_SPEC.get(fg_type) or {}
    extra = spec.get("extra_block") or ()
    self_allow = set(spec.get("self_allow") or ())
    return any(info.get(k) for k in _FG_BLOCK_KEYS + extra if k not in self_allow)


def _fg_allowed_atoms(on_ring: list[dict], spec: dict, mol: Mol | None = None) -> set[int]:
    ckey = spec["ckey"]
    if spec.get("ekey") == "hydroxyls":
        return {e["o_idx"] for e in on_ring}
    if spec.get("ekey") == "amines":
        return {e["n_idx"] for e in on_ring}
    extra_fn = spec.get("extra_atoms_fn")
    fn = _EXTRA_ATOMS_FNS.get(extra_fn) if (extra_fn and mol) else None
    base = set().union(*(fn(mol, e[ckey]) for e in on_ring)) if fn else set()
    return base | {e[ckey] for e in on_ring} if fn else set()


def _arene_fg_subs_ok(
    info: dict, mol: Mol, ring: set[int], fg_idxs: set[int], allowed: set[int],
) -> bool:
    """Gate: only halo + simple alkyl side chains; no extra heteroatoms."""
    if not _outside_ok(mol, ring, allowed):
        return False
    starts = _ring_side_starts(mol, ring, fg_idxs)
    return _ring_halo_n(mol, ring) + len(starts) <= 2


def _resolve_arene_core(info, scaffold_id, fg_type):
    """(ring_atoms, chain, spec, [meta]) or None."""
    spec = _FG_SPEC.get(fg_type)
    if spec is None: return None
    entry = _SCAFFOLD_DETECT.get(scaffold_id)
    if entry is None: return None
    core = entry[0](info)
    if core is None: return None
    if len(core) == 3:
        ring_atoms, chain, meta = core
        return (ring_atoms, chain, spec, meta) if chain else None
    ring_atoms, chain = core
    return (ring_atoms, chain, spec) if chain else None


def _rg_attach(mol: Mol, c_idx: int, ring_atoms: set[int]) -> int | None:
    """Return ring atom attached to exocyclic carbon c_idx, or None."""
    nbs = [n.GetIdx() for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
           if n.GetAtomicNum() == 6 and n.GetIdx() in ring_atoms]
    return nbs[0] if len(nbs) == 1 else None


def _check_fg_count(info, spec, ring_atoms, n_fg):
    """Return validated FG entries on ring, or None."""
    fg_entries = info.get(spec["ekey"]) or []
    if spec.get("extra_atoms_fn"):
        mol = info["mol"]
        on_ring = [e for e in fg_entries
                   if _rg_attach(mol, e[spec["ckey"]], ring_atoms) is not None]
    else:
        ckey = spec["ckey"]
        on_ring = [e for e in fg_entries if ckey in e and e[ckey] in ring_atoms]
    if len(on_ring) != n_fg or len(on_ring) != len(fg_entries): return None
    return on_ring


def _gate_scaffold_fg(info, fg_type, ring_atoms, fg_idxs, spec):
    """True if no conflicting FGs and ring substituents are acceptable."""
    if _arene_fg_conflict(info, fg_type): return False
    mol = info["mol"]
    fg_entries = info.get(spec["ekey"]) or []
    if spec.get("extra_atoms_fn"):
        allowed = _fg_allowed_atoms(fg_entries, spec, mol)
    else:
        ckey = spec["ckey"]
        on_ring = [e for e in fg_entries if ckey in e and e[ckey] in ring_atoms]
        allowed = _fg_allowed_atoms(on_ring, spec)
    return _arene_fg_subs_ok(info, mol, ring_atoms, fg_idxs, allowed)


def _arene_fg_attach(info, spec, on_ring, ring_atoms):
    """Compute ring_attach_idx for exocyclic FG; None for direct FG."""
    if not spec.get("extra_atoms_fn") or not on_ring:
        return None
    return _rg_attach(info["mol"], on_ring[0][spec["ckey"]], ring_atoms)


def _assemble_fg_parent(chain, kind, spec, fg_idxs, n_fg, ring_attach=None):
    """Build parent dict with FG locant keys."""
    parent = {"chain": chain, "n_carbons": len(chain), "kind": kind, "scaffold_id": kind}
    if n_fg == 1:
        parent[spec["fk"].replace("_idxs", "_idx")] = fg_idxs[0]
    else:
        parent[spec["fk"]] = fg_idxs
    if ring_attach is not None:
        parent["ring_attach_idx"] = ring_attach
    return parent


def _try_arene_fg_parent(info, scaffold_id, fg_type, n_fg=1):
    """Generalized: FG on aromatic ring -> systematic parent dict."""
    resolved = _resolve_arene_core(info, scaffold_id, fg_type)
    if resolved is None: return None
    if len(resolved) == 4:
        ring_atoms, chain, spec, meta = resolved
    else:
        ring_atoms, chain, spec = resolved
        meta = {}
    if (on_ring := _check_fg_count(info, spec, ring_atoms, n_fg)) is None: return None
    fg_idxs_set, fg_idxs_list = {e[spec["ckey"]] for e in on_ring}, [e[spec["ckey"]] for e in on_ring]
    if not _gate_scaffold_fg(info, fg_type, ring_atoms, fg_idxs_set, spec): return None
    kind = _ARENE_FG_KINDS.get((scaffold_id, fg_type, n_fg))
    if kind is None: return None
    parent = _assemble_fg_parent(chain, kind, spec, fg_idxs_list, n_fg,
                                  _arene_fg_attach(info, spec, on_ring, ring_atoms))
    if meta:
        parent.update(meta)
    return parent


# ---- Thin wrappers for per-scaffold FG parents ----

def _try_naphthalenol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "OH")


def _try_naphthalenediol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "OH", n_fg=2)


def _try_naphthalenamine_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "naphthalene", "NH2")


def _try_quinolinediol_parent(info: dict) -> dict | None:
    return _try_arene_fg_parent(info, "quinoline", "OH", n_fg=2)


def _try_arene_fg_aldehyde(info: dict) -> dict | None:
    for sid in ("naphthalene", "quinoline"):
        if (p := _try_arene_fg_parent(info, sid, "CHO")) is not None:
            return p
    return None


def _try_arene_fg_nitrile(info: dict) -> dict | None:
    for sid in ("naphthalene", "quinoline"):
        if (p := _try_arene_fg_parent(info, sid, "CN")) is not None:
            return p
    return None


def _try_arene_fg_amine(info: dict) -> dict | None:
    """Generalized: NH2 on aromatic ring -> amine parent (naphthalene / heterocycles)."""
    for sid in ("naphthalene", "pyrazole", "thiazole", "quinazoline"):
        if (p := _try_arene_fg_parent(info, sid, "NH2")) is not None:
            return p
    return None
