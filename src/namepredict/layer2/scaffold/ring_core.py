"""Core-only mother-ring identification (no simple/substituent gate).

Every registered fn(info) -> list[tuple[set[int], str]] maps ring-atom sets to
scaffold ids WITHOUT the substituent / outside / FG gates that the full
`_try_*_parent` producers apply. `ring_scaffold._producer_scaffold_ids` uses
these so a substituted mother (e.g. 2,3,6-trimethylquinoline) still resolves
to its mother scaffold instead of being dropped by `_*_subs_ok`.

The full producers remain authoritative for FG-variant retained names
(benzofuranamine, benzothiophenol, ...): `_producer_id` falls back to them
lazily for skeletons the core table did not tag.
"""
from __future__ import annotations

_CORE_FNS: list = []


def _register(fn):
    _CORE_FNS.append(fn)
    return fn


def ring_core_fns() -> list:
    return list(_CORE_FNS)


def _fused56_mono(info: dict, spec, sid: str) -> list:
    from namepredict.layer2.scaffold.fused56 import _mono_parts

    parts = _mono_parts(info, spec)
    return [(set(parts[0]) | set(parts[1]), sid)] if parts else []


def _fused56_di13(info: dict, spec, sid: str) -> list:
    from namepredict.layer2.scaffold.fused56 import _di13_parts

    parts = _di13_parts(info, spec)
    return [(set(parts[0]) | set(parts[1]), sid)] if parts else []


@_register
def _benzene_core(info: dict) -> list:
    from namepredict.layer2.aryl_sub import _arom_c6_ring_lists, _is_unfused_benzene_ring

    mol = info["mol"]
    return [(set(r), "benzene") for r in _arom_c6_ring_lists(mol)
            if _is_unfused_benzene_ring(mol, set(r))]


@_register
def _cycloalkane_core(info: dict) -> list:
    from namepredict.layer2.scaffold.cyclo_pick import _cyclo_parent_candidates

    return [(rs, "cycloalkane") for rs in _cyclo_parent_candidates(info)]


@_register
def _cyclopolyene_core(info: dict) -> list:
    from namepredict.layer2.scaffold.ring_parent import _is_carbocycle_ring, _is_cyclopolyene_core

    if not _is_cyclopolyene_core(info):
        return []
    atoms = _is_carbocycle_ring(info)
    return [(set(atoms), "cyclopolyene")] if atoms else []


@_register
def _pyridine_core(info: dict) -> list:
    from namepredict.layer2.scaffold.pyridine import _pyridine_ring_lists

    return [(set(r), "pyridine") for r in _pyridine_ring_lists(info["mol"])]


@_register
def _diazine_core(info: dict) -> list:
    from namepredict.layer2.scaffold.heteroarene5 import (
        _DIAZINE_KIND,
        _diazine_rings,
        _ring_n_idxs,
        _ring_nn_dist,
    )

    mol = info["mol"]
    out = []
    for r in _diazine_rings(mol):
        ns = _ring_n_idxs(mol, r)
        if len(ns) == 2:
            kind = _DIAZINE_KIND.get(_ring_nn_dist(r, ns))
            if kind:
                out.append((set(r), kind))
    return out


@_register
def _hetero5_core(info: dict) -> list:
    from namepredict.layer2.scaffold.heteroarene5 import _hetero5_kind_of, _hetero5_rings

    mol = info["mol"]
    return [(set(r), k) for r in _hetero5_rings(mol) if (k := _hetero5_kind_of(mol, r))]


@_register
def _azole13_core(info: dict) -> list:
    from namepredict.layer2.scaffold.azole13 import _azole13_kind
    from namepredict.layer2.scaffold.heteroarene5 import _ring_atoms_if_mono

    kind = _azole13_kind(info)
    if kind is None:
        return []
    atom_ids = _ring_atoms_if_mono(info)
    return [(set(atom_ids), kind)] if atom_ids else []


@_register
def _diazole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.heteroarene5 import (
        _imidazole_n_pair,
        _pyrazole_n_pair,
        _ring_atoms_if_mono,
    )

    atom_ids = _ring_atoms_if_mono(info)
    out = []
    if atom_ids:
        if _imidazole_n_pair(info):
            out.append((set(atom_ids), "imidazole"))
        if _pyrazole_n_pair(info):
            out.append((set(atom_ids), "pyrazole"))
    return out


@_register
def _naphthalene_core(info: dict) -> list:
    from namepredict.layer2.scaffold.naphthalene import _is_naphthalene_core, _two_six_rings

    if not _is_naphthalene_core(info):
        return []
    r1, r2 = _two_six_rings(info)
    return [(set(r1) | set(r2), "naphthalene")]


@_register
def _quinoline_core(info: dict) -> list:
    from namepredict.layer2.scaffold.quinoline import _q_core

    parts = _q_core(info)
    if parts is None:
        return []
    r1, r2, _n, _ba, _bb, kind = parts
    return [(set(r1) | set(r2), kind)]  # quinoline | isoquinoline


@_register
def _indole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.indole import _indole_parts

    parts = _indole_parts(info)
    return [(set(parts[0]) | set(parts[1]), "indole")] if parts else []


@_register
def _indazole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.indazole import _iz_parts

    parts = _iz_parts(info)
    return [(set(parts[0]) | set(parts[1]), "indazole")] if parts else []


@_register
def _benzimidazole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.benzimidazole import _bim_parts

    parts = _bim_parts(info)
    return [(set(parts[0]) | set(parts[1]), "benzimidazole")] if parts else []


@_register
def _benzofuran_core(info: dict) -> list:
    from namepredict.layer2.scaffold.builders.fused56 import BF_MONO

    return _fused56_mono(info, BF_MONO, "benzofuran")


@_register
def _benzothiophene_core(info: dict) -> list:
    from namepredict.layer2.scaffold.builders.fused56 import BT_MONO

    return _fused56_mono(info, BT_MONO, "benzothiophene")


@_register
def _benzothiazole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.builders.fused56 import BTZ_DI13

    return _fused56_di13(info, BTZ_DI13, "benzothiazole")


@_register
def _benzoxazole_core(info: dict) -> list:
    from namepredict.layer2.scaffold.builders.fused56 import BOX_DI13

    return _fused56_di13(info, BOX_DI13, "benzoxazole")


@_register
def _quinazoline_core(info: dict) -> list:
    from namepredict.layer2.scaffold.benzodiazine import _qz_core

    parts = _qz_core(info)
    return [(parts[0], "quinazoline")] if parts else []


@_register
def _quinoxaline_core(info: dict) -> list:
    from namepredict.layer2.scaffold.benzodiazine import _qx_core

    parts = _qx_core(info)
    return [(parts[0], "quinoxaline")] if parts else []


@_register
def _anthracene_core(info: dict) -> list:
    from namepredict.layer2.scaffold.anthracene import _anthracene_system

    s = _anthracene_system(info)
    return [(set(s["atom_ids"]), "anthracene")] if s else []


@_register
def _anthraquinone_core(info: dict) -> list:
    from namepredict.layer2.scaffold.anthraquinone import _aq_system

    s = _aq_system(info)
    return [(set(s["atom_ids"]), "anthraquinone")] if s else []


@_register
def _benzoquinone_core(info: dict) -> list:
    from namepredict.layer2.scaffold.benzoquinone import _bq_ring, _ketone_idxs, _ketones_para

    ring = _bq_ring(info)
    if ring is None:
        return []
    ket = _ketone_idxs(info)
    return [(set(ring), "benzoquinone")] if len(ket) == 2 and _ketones_para(ring, ket) else []


@_register
def _ortho_benzoquinone_core(info: dict) -> list:
    from namepredict.layer2.scaffold.benzoquinone import _bq_ring, _ketone_idxs
    from namepredict.layer2.scaffold.ortho_benzoquinone import _ketones_ortho

    ring = _bq_ring(info)
    if ring is None:
        return []
    ket = _ketone_idxs(info)
    return [(set(ring), "ortho_benzoquinone")] if len(ket) == 2 and _ketones_ortho(ring, ket) else []


@_register
def _chromenone_core(info: dict) -> list:
    from namepredict.layer2.scaffold.chromenone import _chrom_core

    parts = _chrom_core(info)
    return [(set(parts[0]) | set(parts[1]), "chromenone")] if parts else []


@_register
def _bridged_core(info: dict) -> list:
    from namepredict.layer2.scaffold.bridged_parent import _is_bridged_system

    return [(set(r["atom_ids"]), "bridged") for r in (info.get("ring_systems") or [])
            if _is_bridged_system(r)]


@_register
def _spiro_core(info: dict) -> list:
    from namepredict.layer2.scaffold.spiro_parent import _is_simple_spiro

    s = _is_simple_spiro(info)
    return [(set(s["atom_ids"]), "spiro")] if s else []


@_register
def _sat_hetero_core(info: dict) -> list:
    from namepredict.layer2.scaffold.sat_hetero import _is_sat_hetero_core, _kind_of, _ring_atoms_if_mono

    if not _is_sat_hetero_core(info):
        return []
    kind = _kind_of(info)
    atom_ids = _ring_atoms_if_mono(info)
    return [(set(atom_ids), kind)] if kind and atom_ids else []


@_register
def _sat_hetero_repl_core(info: dict) -> list:
    from namepredict.layer2.scaffold.sat_hetero import _ring_atoms_if_mono
    from namepredict.layer2.scaffold.builders.sat_hetero_repl import _is_repl_core

    if not _is_repl_core(info):
        return []
    atom_ids = _ring_atoms_if_mono(info)
    return [(set(atom_ids), "sat_hetero_repl")] if atom_ids else []
