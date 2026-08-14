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
# 组合 kind 已根除（scaffold×FG 正交化）：集合只含活的 FG 类别 kind。
_OH_KINDS = ("alcohol",)
_AMINE_KINDS = ("amine", "sec_amine", "tert_amine")
def _typed_atom_locants(oriented: dict, group: str) -> list[int]:
    chain = oriented.get("chain") or []
    atoms = _typed_group_atoms(oriented, group)
    kind, facts = oriented.get("kind"), oriented.get("numbering_scaffold")
    required = oriented.get("numbering_scaffold_required", False)
    return sorted(loc for atom in atoms
                  if (loc := _atom_locant(chain, atom, kind, facts, required)) is not None)


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


def _oh_locants(oriented: dict) -> list[int] | None:
    locs = _typed_atom_locants(oriented, "alcohol")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("alcohol", "benzenediol"), "oh_c_idxs")
def _amine_pair_locants(oriented: dict) -> list[int] | None:
    locs = _typed_atom_locants(oriented, "amine")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("amine",), "amine_c_idxs")
def _amine_locant(oriented: dict) -> int | None:
    locs = _typed_atom_locants(oriented, "amine")
    return locs[0] if len(locs) == 1 else _fg_locant(oriented, _AMINE_KINDS, "amine_c_idx")
def _ketone_locant(oriented: dict) -> int | None:
    atoms = _typed_group_atoms(oriented, "ketone")
    return _chain_pos(oriented.get("chain") or [], atoms[0]) if len(atoms) == 1 else _fg_locant(oriented, ("ketone",), "ketone_c_idx")
def _ketone_pair_locants(oriented: dict) -> list[int] | None:
    atoms = _typed_group_atoms(oriented, "ketone")
    if len(atoms) > 1:
        return _pair_locants(oriented.get("chain") or [], atoms)
    return _oriented_pair_locants(oriented, ("dione",), "ketone_c_idxs")
def _has_parent_ene(oriented: dict) -> bool:
    return bool(oriented.get("double_bond") or oriented.get("double_bonds"))
def _ene_locant(oriented: dict) -> int | None:
    if not _has_parent_ene(oriented):
        return None
    db = oriented.get("double_bond")
    if db is None:
        return None  # 多烯：位次由 ene_locants 列表承载
    return _edge_min_locant(oriented.get("chain") or [], db)
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
def _sh_locants_list(oriented: dict) -> list[int] | None:
    loc = _sh_locant(oriented)
    return [loc] if loc is not None else None
def _omit_ket_loc(oriented: dict, n_subs: int) -> bool:
    single = len(_typed_group_atoms(oriented, "ketone")) == 1
    return _omit_ketone(oriented.get("kind"), n_subs, oriented,
                        has_ene=_has_parent_ene, single=single)
def _omit_for(kind: str, oriented: dict, n: int, n_subs: int) -> bool:
    """FG 记录 omit 标志:环状判断由 omit_locants 基于 scaffold_id 完成（不虚构 cyclo* kind）。"""
    if kind == "oh":
        return _omit_oh(_oh_locant(oriented), n, oriented.get("kind"), oriented, n_subs)
    if kind == "amine":
        return _omit_amine(_amine_locant(oriented), n, oriented.get("kind"), n_subs, oriented)
    if kind == "ketone":
        return _omit_ket_loc(oriented, n_subs)
    if kind == "sh":
        return _omit_sh(_sh_locant(oriented), n)
    return False
def _amine_fg_locants(oriented: dict) -> list[int] | None:
    locs = _amine_pair_locants(oriented) or []
    if locs:
        return locs
    loc = _amine_locant(oriented)
    return [loc] if loc is not None else None
def _ketone_fg_locants(oriented: dict) -> list[int] | None:
    locs = _ketone_pair_locants(oriented) or []
    if locs:
        return locs
    loc = _ketone_locant(oriented)
    return [loc] if loc is not None else None
def _acid_fg_locants(oriented: dict) -> list[int] | None:
    # exocyclic 酸的羧基碳在环外；locant 取环上附着原子 ring_attach_idx。
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


_FG_LOCANTS = (
    ("oh", _oh_locants),
    ("amine", _amine_fg_locants),
    ("ketone", _ketone_fg_locants),
    ("sh", _sh_locants_list),
    ("acid", _acid_fg_locants),
)
def _fg_locants(oriented: dict, n_subs: int = 0) -> list[dict]:
    """FG 位次记录: [{kind, locants, omit}] — 稀疏,只产实际存在的 principal FG."""
    n = oriented.get("n_carbons", 0)
    records = []
    for kind, locs_fn in _FG_LOCANTS:
        locs = locs_fn(oriented)
        if not locs:
            continue
        records.append({
            "kind": kind, "locants": sorted(locs), "omit": _omit_for(kind, oriented, n, n_subs),
        })
    return records
def _pack(oriented: dict, substituents: list) -> dict:
    from namepredict.layer4.cyclo_relative_stereo import relative_stereo_facts
    facts = relative_stereo_facts(oriented)
    merged = {**oriented, **facts}
    n_subs = len(substituents or [])
    return {
        "parent": merged, "substituents": substituents,
        "fg_locants": _fg_locants(merged, n_subs),
        **_unsat_locants(merged, merged.get("n_carbons", 0)), **facts,
    }
