"""Core-only mother-ring identification (no simple/substituent gate).

Splits 'is this ring system mother X' from 'is molecule simple enough to use
retained name X'. Every registered fn(info) -> list[tuple[set[int], str]]
maps ring-atom sets to scaffold ids WITHOUT the substituent / outside / FG
gates that the full `_try_*_parent` producers apply.

Two sources:
  * `_template_mother` — retained mothers registered as SMILES templates and
    matched by subgraph isomorphism (see retained_templates), covering the
    aromatic / heteroarene / fused mothers.
  * dedicated adapters — carbonyl mothers (benzoquinone / anthraquinone /
    chromenone / ortho_benzoquinone) whose templates would include exocyclic
    =O atoms, and rule mothers (cycloalkane / cyclopolyene / bridged / spiro /
    saturated hetero) whose composition is open-ended.

The full producers are NOT consulted for scaffold identity. A former fallback
over `ring_try_fns` (in `_producer_id`) never hit on real data — the template
table already covers every mother ring the producers can name, including the
FG-variant skeletons (benzofuranamine's ring set IS benzofuran) — and it
clobbered the carbonyl-mother adapters above, so it was removed.
"""
from __future__ import annotations

_CORE_FNS: list = []


def _register(fn):
    _CORE_FNS.append(fn)
    return fn


def ring_core_fns() -> list:
    return list(_CORE_FNS)


@_register
def _template_mother(info: dict) -> list:
    """Retained mothers via SMILES-template subgraph isomorphism."""
    from namepredict.layer2.scaffold.retained_templates import match_retained

    out: list = []
    for s in info.get("ring_systems") or []:
        atoms = tuple(sorted(s["atom_ids"]))
        sid = match_retained(info, atoms)
        if sid:
            out.append((set(atoms), sid))
    return out


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
def _anthraquinone_core(info: dict) -> list:
    from namepredict.layer2.scaffold.anthraquinone import _aq_system

    s = _aq_system(info)
    return [(set(s["atom_ids"]), "anthraquinone")] if s else []


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
    from namepredict.layer2.scaffold.builders.sat_hetero_repl import _is_repl_core
    from namepredict.layer2.scaffold.sat_hetero import _ring_atoms_if_mono

    if not _is_repl_core(info):
        return []
    atom_ids = _ring_atoms_if_mono(info)
    return [(set(atom_ids), "sat_hetero_repl")] if atom_ids else []
