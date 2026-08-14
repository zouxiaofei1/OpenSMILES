from __future__ import annotations
from namepredict.layer4._chain_orient import _chain_pos, _edge_min_locant, _pair_locants
from namepredict.layer4.locants.adapt import effective_sub_locant, plan_from_chain
from namepredict.layer4.omit_locants import (
    omit_amine as _omit_amine, omit_ketone as _omit_ketone, omit_sh as _omit_sh,
)
from namepredict.layer4.orienters import _typed_group_atoms
from namepredict.layer4.polyene import ene_locants
def _atom_locant(chain: list[int], atom: int | None, kind: str | None, facts=None, required=False) -> int | None:
    """Scaffold facts use a plan; otherwise retain ordinary chain numbering."""
    if atom is None or atom not in chain:
        return None
    plan = plan_from_chain(chain, kind, facts, required=required)
    loc = effective_sub_locant(plan, atom) if plan else None
    return loc if loc is not None else chain.index(atom) + 1
def _fg_locant(oriented: dict, kinds: tuple, key: str) -> int | None:
    if oriented.get("kind") not in kinds:
        return None
    return _atom_locant(oriented.get("chain") or [], oriented.get(key), oriented.get("kind"), oriented.get("numbering_scaffold"), oriented.get("numbering_scaffold_required", False))
_OH_KINDS = ("alcohol", "cycloalcohol", "benzothiophenol", "naphthalenol")
_AMINE_KINDS = (
    "amine", "cycloamine", "sec_amine", "tert_amine",
    "benzofuranamine", "benzothiazolamine",
    "benzoxazolamine", "benzimidazolamine",
    "pyrazolamine", "thiazolamine", "quinazolinamine",
)
def _typed_atom_locants(oriented: dict, group: str) -> list[int]:
    chain = oriented.get("chain") or []
    atoms = _typed_group_atoms(oriented, group)
    kind, facts = oriented.get("kind"), oriented.get("numbering_scaffold")
    required = oriented.get("numbering_scaffold_required", False)
    return sorted(loc for atom in atoms
                  if (loc := _atom_locant(chain, atom, kind, facts, required)) is not None)


def _typed_alcohol_atoms(oriented: dict) -> list[int]:
    return _typed_group_atoms(oriented, "alcohol")


def _oh_locant(oriented: dict) -> int | None:
    locs = _typed_atom_locants(oriented, "alcohol")
    return locs[0] if len(locs) == 1 else _fg_locant(oriented, _OH_KINDS, "oh_c_idx")
def _sh_locant(oriented: dict) -> int | None:
    return _fg_locant(oriented, ("thiol",), "sh_c_idx")
def _oriented_pair_locants(oriented: dict, kinds, key: str) -> list[int] | None:
    if oriented.get("kind") not in kinds:
        return None
    locs = _pair_locants(oriented.get("chain") or [], oriented.get(key))
    return list(locs) if locs else None


def _scaffold_pair_locants(oriented: dict, key: str) -> list[int] | None:
    """Multi-FG locants with scaffold-aware numbering (fused rings)."""
    chain = oriented.get("chain") or []
    cs = oriented.get(key)
    if not cs or not chain:
        return None
    kind = oriented.get("kind")
    facts = oriented.get("numbering_scaffold")
    required = oriented.get("numbering_scaffold_required", False)
    locs = [_atom_locant(chain, c, kind, facts, required) for c in cs if c in chain]
    return sorted(locs) if len(locs) == len(cs) and all(l is not None for l in locs) else None
def _oh_locants(oriented: dict) -> list[int] | None:
    locs = _typed_atom_locants(oriented, "alcohol")
    if locs:
        return locs
    if oriented.get("kind") in ("naphthalenediol",):
        return _scaffold_pair_locants(oriented, "oh_c_idxs")
    return _oriented_pair_locants(oriented, ("alcohol", "diol", "triol", "benzenediol", "cycloalkanediol"), "oh_c_idxs")
def _typed_amine_atoms(oriented: dict) -> list[int]:
    return _typed_group_atoms(oriented, "amine")


def _amine_pair_locants(oriented: dict) -> list[int] | None:
    locs = _typed_atom_locants(oriented, "amine")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("diamine", "triamine", "tetraamine"), "amine_c_idxs")
def _amine_locant(oriented: dict) -> int | None:
    locs = _typed_atom_locants(oriented, "amine")
    return locs[0] if len(locs) == 1 else _fg_locant(oriented, _AMINE_KINDS, "amine_c_idx")
def _ketone_locant(oriented: dict) -> int | None:
    atoms = _typed_group_atoms(oriented, "ketone")
    return _chain_pos(oriented.get("chain") or [], atoms[0]) if len(atoms) == 1 else _fg_locant(oriented, ("ketone", "cycloketone"), "ketone_c_idx")
def _ketone_pair_locants(oriented: dict) -> list[int] | None:
    atoms = _typed_group_atoms(oriented, "ketone")
    if len(atoms) > 1:
        return _pair_locants(oriented.get("chain") or [], atoms)
    return _oriented_pair_locants(
        oriented, ("dione", "benzoquinone", "ortho_benzoquinone", "cycloalkanedione"), "ketone_c_idxs",
    )
def _has_parent_ene(oriented: dict) -> bool:
    return bool(oriented.get("double_bond") or oriented.get("double_bonds"))
def _ene_locant(oriented: dict) -> int | None:
    kind = oriented.get("kind")
    if kind in ("alkene", "cycloalkene") or _has_parent_ene(oriented):
        return _edge_min_locant(
            oriented.get("chain") or [], oriented.get("double_bond")
        )
    return None
def _yne_locant(oriented: dict) -> int | None:
    locs = oriented.get("yne_locants") or []
    if locs:
        return locs[0]
    if oriented.get("kind") == "alkyne" or oriented.get("triple_bond"):
        return _edge_min_locant(oriented.get("chain") or [], oriented.get("triple_bond"))
    return None
def _has_parent_yne(oriented: dict) -> bool:
    return bool(oriented.get("triple_bond") or oriented.get("triple_bonds"))
def _omit_oh(oh_pos, n_carbons, kind=None, parent=None, n_subs=0):
    from namepredict.layer4.omit_locants import omit_oh as _core
    return _core(oh_pos, n_carbons, kind, parent, n_subs,
                 has_ene=_has_parent_ene, has_yne=_has_parent_yne)
def _omit_unsat(n_carbons, kind=None, parent=None):
    from namepredict.layer4.omit_locants import omit_unsat as _core
    return _core(n_carbons, kind, parent, has_ene=_has_parent_ene, has_yne=_has_parent_yne)
def _sub_locant(chain: list[int], attach: int, kind: str | None, facts=None) -> int:
    loc = _atom_locant(chain, attach, kind, facts)
    return 0 if loc is None else loc
def _with_locants(chain: list[int], substituents: list, kind: str | None = None, facts=None) -> list:
    return [{**s, "locant": _sub_locant(chain, s["attach_idx"], kind, facts)} for s in substituents]
def _unsat_locants(oriented: dict, n: int) -> dict:
    kind = oriented.get("kind")
    return {
        "ene_locant": _ene_locant(oriented),
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": _omit_unsat(n, kind, oriented),
        "yne_locant": _yne_locant(oriented),
        "omit_yne_locant": _omit_unsat(n, kind, oriented),
    }
def _typed_ring_numbering_kind(oriented: dict, atoms: list[int], cyclo_kind: str) -> str:
    return cyclo_kind if len(atoms) == 1 and oriented.get("kind") == "cycloalkane" else oriented.get("kind")


def _oh_am_locants(oriented: dict, n: int, n_subs: int = 0) -> dict:
    oh, am = _oh_locant(oriented), _amine_locant(oriented)
    kind = _typed_ring_numbering_kind(oriented, _typed_alcohol_atoms(oriented), "cycloalcohol")
    amine_kind = _typed_ring_numbering_kind(oriented, _typed_amine_atoms(oriented), "cycloamine")
    return {
        "oh_locant": oh, "oh_locants": _oh_locants(oriented),
        "omit_oh_locant": _omit_oh(oh, n, kind, oriented, n_subs),
        "amine_locant": am, "amine_locants": _amine_pair_locants(oriented),
        "omit_amine_locant": _omit_amine(am, n, amine_kind, n_subs),
    }
def _sh_locants(oriented: dict, n: int) -> dict:
    sh = _sh_locant(oriented)
    return {"sh_locant": sh, "omit_sh_locant": _omit_sh(sh, n)}
def _cooh_locants(oriented: dict) -> list[int] | None:
    if oriented.get("kind") not in {"polycarboxylic", "benzene_polycarboxylic"}: return None
    return _pair_locants(oriented.get("chain") or [], oriented.get("cooh_c_idxs"))
def _omit_ket_loc(oriented: dict, n_subs: int) -> bool:
    return _omit_ketone(oriented.get("kind"), n_subs, oriented, has_ene=_has_parent_ene)
def _fg_locants(oriented: dict, n_subs: int = 0) -> dict:
    n = oriented.get("n_carbons", 0)
    return {
        **_oh_am_locants(oriented, n, n_subs), **_sh_locants(oriented, n),
        "ketone_locant": _ketone_locant(oriented),
        "ketone_locants": _ketone_pair_locants(oriented),
        "omit_ketone_locant": _omit_ket_loc(oriented, n_subs),
        "cooh_locants": _cooh_locants(oriented),
        **_unsat_locants(oriented, n),
    }
def _pack(oriented: dict, substituents: list) -> dict:
    from namepredict.layer4.cyclo_relative_stereo import relative_stereo_facts
    facts = relative_stereo_facts(oriented)
    n_subs = len(substituents or [])
    return {
        "parent": {**oriented, **facts}, "substituents": substituents,
        **_fg_locants({**oriented, **facts}, n_subs), **facts,
    }
