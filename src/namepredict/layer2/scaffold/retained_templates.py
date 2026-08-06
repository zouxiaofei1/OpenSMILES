"""Retained mother-ring templates matched by subgraph isomorphism (P-22 / P-25).

Each retained mother (benzene, quinoline, benzimidazole, ...) is registered as
one aromatic SMILES template. `match_retained(info, atom_ids)` finds the mother
whose template subgraph-isomorphically covers the given ring-atom set exactly
(`set(match) == atom_ids`), which is robust to exocyclic substituents and to
double-bond placement on rings (both are graph-isomorphic).

Positional isomers that are graph-isomorphic (quinoline/isoquinoline,
pyridazine/pyrimidine/pyrazine, imidazole/pyrazole, benzimidazole/indazole,
quinazoline/quinoxaline) hit several templates at once and are disambiguated by
a small positional check on the graph, reusing the existing core detectors.
"""
from __future__ import annotations

from collections import Counter

from rdkit.Chem import Mol, MolFromSmiles

# ── template registry: scaffold_id -> aromatic SMILES ────────────────────────
_TEMPLATES: dict[str, str] = {
    # carbocycles
    "benzene":      "c1ccccc1",
    "naphthalene":  "c1ccc2ccccc2c1",
    "anthracene":   "c1ccc2cc3ccccc3cc2c1",
    # monocyclic heteroarenes
    "furan":        "c1ccoc1",
    "thiophene":    "c1ccsc1",
    "pyrrole":      "c1cc[nH]c1",
    "pyridine":     "n1ccccc1",
    "pyridazine":   "c1ccnnc1",
    "pyrimidine":   "c1cncnc1",
    "pyrazine":     "c1cnccn1",
    "imidazole":    "c1cnc[nH]1",
    "pyrazole":     "c1ccn[nH]1",
    "oxazole":      "c1cocn1",
    "thiazole":     "c1cscn1",
    # fused 5+6
    "indole":        "c1ccc2[nH]ccc2c1",
    "indazole":      "c1ccc2cn[nH]c2c1",
    "benzimidazole": "c1ccc2[nH]cnc2c1",
    "benzofuran":    "c1ccc2occc2c1",
    "benzothiophene":"c1ccc2sccc2c1",
    "benzothiazole": "c1ccc2scnc2c1",
    "benzoxazole":   "c1ccc2ocnc2c1",
    "quinoline":     "c1ccc2ncccc2c1",
    "isoquinoline":  "c1nccc2ccccc21",
    "quinazoline":   "c1ccc2ncncc2c1",
    "quinoxaline":   "c1ccc2nccnc2c1",
    # NOTE: carbonyl mothers (benzoquinone / anthraquinone / chromenone /
    # ortho_benzoquinone) are NOT here: their templates include exocyclic =O
    # atoms, so the match set would exceed the ring-system atom set. They keep
    # dedicated core adapters in ring_core.
}

# Query mols and element signatures, built once at import.
_Q: dict[str, Mol] = {sid: MolFromSmiles(smi) for sid, smi in _TEMPLATES.items()}


def _elem_sig(mol: Mol) -> frozenset:
    return frozenset(Counter(a.GetAtomicNum() for a in mol.GetAtoms()).items())


_TEMPLATE_ELEM: dict[str, frozenset] = {sid: _elem_sig(q) for sid, q in _Q.items()}


# ── positional-isomer disambiguation (graph-isomorphic groups) ───────────────
def _n_ring_dist(info: dict, atoms, dist_to_sid: dict[int, str]) -> str | None:
    """Single aromatic ring: return sid by the ring distance of its two N."""
    mol = info["mol"]
    for r in mol.GetRingInfo().AtomRings():
        if set(r) != atoms:
            continue
        from namepredict.layer2.scaffold.heteroarene5 import (
            _ring_n_idxs,
            _ring_nn_dist,
        )
        ns = _ring_n_idxs(mol, list(r))
        if len(ns) == 2:
            return dist_to_sid.get(_ring_nn_dist(list(r), ns))
    return None


def _dis_diazine(info: dict, atoms) -> str | None:
    return _n_ring_dist(info, atoms, {1: "pyridazine", 2: "pyrimidine", 3: "pyrazine"})


def _dis_diazole(info: dict, atoms) -> str | None:
    return _n_ring_dist(info, atoms, {1: "pyrazole", 2: "imidazole"})


def _dis_quinoline(info: dict, atoms) -> str | None:
    from namepredict.layer2.scaffold.quinoline import _q_core

    parts = _q_core(info)
    if parts is None or set(parts[0]) | set(parts[1]) != atoms:
        return None
    return parts[5]  # "quinoline" | "isoquinoline"


def _dis_bim_iz(info: dict, atoms) -> str | None:
    from namepredict.layer2.scaffold.benzimidazole import _bim_parts
    from namepredict.layer2.scaffold.indazole import _iz_parts

    parts = _bim_parts(info)
    if parts is not None and set(parts[0]) | set(parts[1]) == atoms:
        return "benzimidazole"
    parts = _iz_parts(info)
    if parts is not None and set(parts[0]) | set(parts[1]) == atoms:
        return "indazole"
    return None


def _dis_qz_qx(info: dict, atoms) -> str | None:
    from namepredict.layer2.scaffold.benzodiazine import _qx_core, _qz_core

    parts = _qz_core(info)
    if parts is not None and parts[0] == atoms:
        return "quinazoline"
    parts = _qx_core(info)
    if parts is not None and parts[0] == atoms:
        return "quinoxaline"
    return None


_DISAMBIGUATE: dict[tuple[str, ...], callable] = {
    ("quinoline", "isoquinoline"): _dis_quinoline,
    ("pyridazine", "pyrimidine", "pyrazine"): _dis_diazine,
    ("imidazole", "pyrazole"): _dis_diazole,
    ("benzimidazole", "indazole"): _dis_bim_iz,
    ("quinazoline", "quinoxaline"): _dis_qz_qx,
}


# ── entry point ───────────────────────────────────────────────────────────────
def match_retained(info: dict, atom_ids) -> str | None:
    """Return the retained mother sid whose template covers atom_ids exactly.

    Robust to exocyclic substituents (templates ignore atoms outside the ring)
    and ring double-bond placement (graph-isomorphic). Positional isomers are
    resolved by `_DISAMBIGUATE`. Element signature prefilter skips templates
    whose composition cannot match before running subgraph isomorphism.
    """
    mol = info["mol"]
    atoms = frozenset(atom_ids)
    elem = frozenset(Counter(mol.GetAtomWithIdx(i).GetAtomicNum() for i in atom_ids).items())
    hits: list[str] = []
    for sid, q in _Q.items():
        if _TEMPLATE_ELEM[sid] != elem:
            continue
        if any(set(m) == atoms for m in mol.GetSubstructMatches(q, uniquify=True)):
            hits.append(sid)
    if not hits:
        return None
    if len(hits) == 1:
        return hits[0]
    for group, fn in _DISAMBIGUATE.items():
        if set(hits) == set(group):
            chosen = fn(info, atoms)
            if chosen in hits:
                return chosen
    return hits[0]
