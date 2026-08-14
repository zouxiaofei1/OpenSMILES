from __future__ import annotations
from namepredict.layer4._chain_orient import (
    _chain_pos,
    _edge_locants,
    _orient_table_cands,
    _pair_locants,
    _pick_ring_by_pair_locants,
    _prefer as _prefer_chain,
    _prefer_lowest_bond_locs,
    _prefer_lowest_locs,
    _rotate_to,
    _table_loc_key,
)

from namepredict.layer4.locants.generate import ring_candidates
from namepredict.layer4.polyene import (
    orient_alkenol, orient_alkenedioic, orient_cycloalkene,
    orient_cyclopolyene, orient_polyene, orient_ring_fg_ene, prefer_unsat_if_fg_tie,
)
from namepredict.layer4.sat_hetero_orient import sat_hetero_orienters as _sat_hetero_orienters
def _maybe_reverse(chain: list[int], pos: int) -> list[int]:
    if pos > (len(chain) + 1) // 2:
        return list(reversed(chain))
    return chain
def _orient_alkane(chain: list[int], substituents: list) -> list[int]:
    if not substituents or not chain:
        return chain
    return _prefer_chain(chain, list(reversed(chain)), substituents)
def _orient_by_single_fg(chain: list[int], parent: dict, substituents: list, key: str) -> list[int]:
    pos = _chain_pos(chain, parent.get(key))
    if pos is None: return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    b, r = _chain_pos(base, parent.get(key)), _chain_pos(rev, parent.get(key))
    if b is not None and r is not None and b == r:
        return _prefer_chain(base, rev, substituents)
    return base
def _orient_alcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    if parent.get("double_bonds"): return orient_alkenol(chain, parent, substituents, _prefer_chain)
    cs = parent.get("oh_c_idxs")
    if cs is not None:
        return _orient_pair(chain, parent, "oh_c_idxs", substituents)
    return _orient_by_single_fg(chain, parent, substituents, "oh_c_idx")
def _orient_thiol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "sh_c_idx")
def _orient_amine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "amine_c_idx")
def _orient_ketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    base = _orient_by_single_fg(chain, parent, substituents, "ketone_c_idx")
    return prefer_unsat_if_fg_tie(base, parent, substituents, "ketone_c_idx", _prefer_ene_orient)
def _better_pair_orient(a: list[int], b: list[int], cs, subs: list) -> list[int]:
    return _prefer_lowest_locs(
        a, b, lambda x: _pair_locants(x, cs), lambda x, y: _prefer_chain(x, y, subs),
    )
def _orient_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    return chain if not cs else _better_pair_orient(chain, list(reversed(chain)), cs, subs)
def _orient_polyol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_diamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "amine_c_idxs", substituents)
def _orient_dione(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "ketone_c_idxs", substituents)
def _orient_to_terminal(chain: list[int], c_idx: int | None) -> list[int]:
    pos = _chain_pos(chain, c_idx)
    return chain if pos is None or pos == 1 else list(reversed(chain))
def _acid_anchor(parent: dict) -> int | None:
    facts = parent.get("principal_expression_facts")
    atoms = sorted(facts.attachment_atoms) if facts and facts.group_class.value == "acid" else []
    return atoms[0] if len(atoms) == 1 else parent.get("cooh_c_idx")


def _orient_acid(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, _acid_anchor(parent))
def _orient_aldehyde(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("aldehyde_c_idx"))
def _orient_ester(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("ester_c_idx"))
def _orient_amide(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("amide_c_idx"))
def _orient_nitrile(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("nitrile_c_idx"))
def _orient_acyl_chloride(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("acyl_c_idx") or parent.get("c_attach"))
def _prefer_ene_orient(a, b, ends0, subs):
    """Pick chain direction with lower unsaturation locant (P-31.1)."""
    return _prefer_lowest_bond_locs(a, b, [ends0], subs, _prefer_chain)
def _orient_by_bond(chain, parent, substituents, key):
    ends0 = parent.get(key)
    if _edge_locants(chain, ends0) is None: return chain
    return _prefer_ene_orient(chain, list(reversed(chain)), ends0, substituents)
def _orient_alkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "double_bond")
def _orient_cycloalkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_cycloalkene(chain, parent, substituents, _prefer_chain)
def _orient_alkyne(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "triple_bond")
def _term_fn(fn):
    return lambda c, p, s: fn(c, p)
def _terminal_orienters() -> dict:
    ac = _term_fn(_orient_acyl_chloride)
    return {
        "acid": _term_fn(_orient_acid), "aldehyde": _term_fn(_orient_aldehyde),
        "ester": _term_fn(_orient_ester), "amide": _term_fn(_orient_amide),
        "nitrile": _term_fn(_orient_nitrile), "sulfonyl_chloride": ac,
        "acyl_chloride": ac, "acyl_bromide": ac,
    }
def _carbonyl_orienters() -> dict:
    return {**_terminal_orienters(), "ketone": _orient_ketone, "dione": _orient_dione}
def _best_ring(chain: list[int], substituents: list) -> list[int]:
    best = chain
    for cand in ring_candidates(chain):
        best = _prefer_chain(best, cand, substituents)
    return best
def _orient_cycloalkane(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return chain if not chain or not substituents else _best_ring(chain, substituents)
def _orient_ring_fixed(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    c = parent.get(key)
    if c is None or c not in chain:
        return chain
    base = _rotate_to(chain, c)
    rev = _rotate_to(list(reversed(chain)), c)
    return _prefer_chain(base, rev, substituents)
def _orient_ring_fg_ene(chain, parent, substituents, key):
    return orient_ring_fg_ene(
        chain, parent, substituents, key, _prefer_chain, _prefer_ene_orient,
    )
def _orient_cycloalcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fg_ene(chain, parent, substituents, "oh_c_idx")
def _orient_ring_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    if not cs or not chain:
        return chain
    return _pick_ring_by_pair_locants(chain, cs, subs, ring_candidates, _prefer_chain)
def _orient_benzene_polycarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "cooh_c_idxs", substituents)
def _orient_benzenediol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_diazine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "n_idxs", substituents)
def _orient_cycloketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fg_ene(chain, parent, substituents, "ketone_c_idx")
def _orient_cycloamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "amine_c_idx")
def _orient_benzoic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "ring_attach_idx")
def _orient_pyridine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "n_idx")
def _orient_hetero5(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "hetero_idx")
def _n_loc_on(chain: list[int], n: int | None) -> int:
    return chain.index(n) + 1 if n is not None and n in chain else 99
def _orient_imidazole(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """NH = 1; choose direction so the other N has the lowest locant (3 vs 4; 2 vs 5)."""
    nh = parent.get("nh_idx")
    if nh is None or nh not in chain:
        return chain
    base = _rotate_to(chain, nh)
    rev = _rotate_to(list(reversed(chain)), nh)
    n = parent.get("n_idx")
    return base if _n_loc_on(base, n) <= _n_loc_on(rev, n) else rev
_NAPH_LOCANTS = (1, 2, 3, 4, None, 5, 6, 7, 8, None)


def _orient_naphthalene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_table_cands(
        chain, parent, substituents, "naph_chains",
        lambda c, s: _table_loc_key(c, s, _NAPH_LOCANTS),
    )
def _orient_indole(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return chain
def _aza_orienters() -> dict:
    return {
        "pyridine": _orient_pyridine,
        "imidazole": _orient_imidazole, "pyrazole": _orient_imidazole,
        "pyrazolamine": _orient_imidazole,
        "oxazole": _orient_imidazole, "thiazole": _orient_imidazole,
        "pyrimidine": _orient_diazine, "pyrazine": _orient_diazine,
        "pyridazine": _orient_diazine}
_FIXED_FUSED = (
    "indole", "indazole",
    "benzofuran", "benzofuranamine", "benzothiophene", "benzothiophenol", "benzothiazole",
    "benzothiazolamine", "benzoxazole", "benzoxazolamine", "benzimidazole", "benzimidazolamine",
    "quinoline", "isoquinoline", "chromenone",
    "quinazoline", "quinazolinamine",
)
def _fused_orienters() -> dict:
    d = {kind: _orient_indole for kind in _FIXED_FUSED}
    d["naphthalene"] = _orient_naphthalene
    d["anthraquinone"] = None
    return d
def _arene_orienters() -> dict:
    return {
        "furan": _orient_hetero5, "thiophene": _orient_hetero5, "pyrrole": _orient_hetero5,
        "thiazolamine": _orient_hetero5,
        **_fused_orienters(), **_aza_orienters(), **_sat_hetero_orienters(),
    }
def _orient_ring_ketone_pair(c, p, s):
    return _orient_ring_pair(c, p, "ketone_c_idxs", s)
def _hetero_orienters() -> dict:
    bq, di = _orient_ring_ketone_pair, _orient_benzenediol
    return {
        "alcohol": _orient_alcohol, "thiol": _orient_thiol,
        "diol": _orient_polyol, "triol": _orient_polyol, "diamine": _orient_diamine, "triamine": _orient_diamine, "tetraamine": _orient_diamine,
        "cycloalcohol": _orient_cycloalcohol, "phenol": _orient_cycloalcohol,
        "benzenediol": di, "cycloalkanediol": di,
        "benzoquinone": bq, "ortho_benzoquinone": bq, "amine": _orient_amine,
        "cycloamine": _orient_cycloamine,
        "aniline": _orient_cycloamine, "sec_amine": _orient_amine,
        "tert_amine": _orient_amine, **_arene_orienters()}
def _orient_polyene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_polyene(chain, parent, substituents, _prefer_chain)
def _orient_alkenedioic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_alkenedioic(chain, parent, substituents, _prefer_chain, _orient_alkene)
def _orient_diacid(chain: list[int], parent: dict, substituents: list) -> list[int]:
    if parent.get("double_bond") or parent.get("double_bonds"):
        return _orient_alkenedioic(chain, parent, substituents)
    return _orient_alkane(chain, substituents)
def _orient_polycarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return None
def _benzoic_orienters() -> dict:
    b = _orient_benzoic
    return {"benzoic": b, "benzaldehyde": b, "acetophenone": b, "benzoate": b,
            "benzonitrile": b, "benzoyl_chloride": b, "benzoyl_bromide": b,
            "benzamide": b}
def _cyclo_exo_orienters() -> dict:
    """Monocyclic cycloalkane + one exocyclic carbonyl FG (ring_attach_idx = 1)."""
    b = _orient_benzoic
    return {"cycloalkanecarboxylic": b}
def _append_bridge_paths(order: list[int], paths: list[list[int]], b: int) -> None:
    """Append bridge paths to order, alternating direction; append b after first path."""
    for i, path in enumerate(paths):
        order.extend(path if i % 2 == 0 else list(reversed(path)))
        if i == 0:
            order.append(b)


def _orient_bridged(chain: list[int], parent: dict, substituents: list) -> list[int]:
    """von Baeyer numbering: bh1, longest bridge → bh2, next longest ← bh1, ..."""
    bh, paths = parent.get("bridgeheads"), parent.get("bridge_paths")
    if not bh or len(bh) < 2 or not paths:
        return chain
    order = [bh[0]]
    _append_bridge_paths(order, paths, bh[1])
    return order


def _bridged_orienters() -> dict:
    return {"bridged": _orient_bridged}


def _unsat_orienters() -> dict:
    bq = _orient_ring_ketone_pair
    return {"cycloketone": _orient_cycloketone, "cycloalkanedione": bq,
            "alkene": _orient_alkene, "diacid": _orient_diacid, "diester": _orient_diacid,
            "polycarboxylic": _orient_polycarboxylic,
            "benzene_polycarboxylic": _orient_benzene_polycarboxylic,
            "polyene": _orient_polyene,
            "cyclopolyene": orient_cyclopolyene, "cycloalkene": _orient_cycloalkene,
            "alkyne": _orient_alkyne, "cycloalkane": _orient_cycloalkane, "benzene": _orient_cycloalkane,
            **_benzoic_orienters(), **_cyclo_exo_orienters(), **_bridged_orienters()}
def _orient_radical(chain: list[int], parent: dict, subs: list) -> list[int]:
    """Free-radical parent: anchor carbon (radical_c_idx) is locant 1; choose the
    lower-locant direction round the ring."""
    return _orient_ring_fixed(chain, parent, subs, "radical_c_idx")
def _kind_orienters() -> dict:
    return {**_hetero_orienters(), **_unsat_orienters(), **_carbonyl_orienters(),
            "phenyl": _orient_radical}
def _orient_by_kind(kind: str, chain: list[int], parent: dict, subs: list) -> list[int]:
    
    fn = _kind_orienters().get(kind)
    print("_orient_by_kind: ",fn)
    #print("orient",kind,chain,parent)
    return _orient_alkane(chain, subs) if fn is None else fn(chain, parent, subs)
def _typed_ring_acid(parent: dict) -> bool:
    facts = parent.get("principal_expression_facts")
    return bool(facts and facts.group_class.value == "acid" and
                facts.relation.value == "exocyclic")


def _typed_group_atoms(parent: dict, group: str) -> list[int]:
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts and facts.group_class.value == group else []


def _typed_alcohol_parent(parent: dict, atoms: list[int]) -> dict:
    fields = {"oh_c_idxs": atoms}
    if len(atoms) == 1:
        fields["oh_c_idx"] = atoms[0]
    return {**parent, **fields}


_TYPED_ALCOHOL_ORIENT_KINDS = frozenset({
    "alcohol", "diol", "triol", "cycloalkane", "cycloalcohol",
    "cycloalkanediol", "benzene", "phenol", "benzenediol",
})


def _typed_alcohol_orient(parent: dict, chain: list[int], subs: list) -> list[int] | None:
    print("_typed_alcohol_orient")
    atoms = _typed_group_atoms(parent, "alcohol")
    if not atoms or parent.get("kind") not in _TYPED_ALCOHOL_ORIENT_KINDS:
        return None
    typed = _typed_alcohol_parent(parent, atoms)
    if parent.get("kind") in {"cycloalkane", "benzene", "phenol", "benzenediol"}:
        return (_orient_ring_fixed(chain, typed, subs, "oh_c_idx") if len(atoms) == 1
                else _orient_ring_pair(chain, typed, "oh_c_idxs", subs))
    return _orient_alcohol(chain, typed, subs) if len(atoms) == 1 else _orient_polyol(chain, typed, subs)


_TYPED_AMINE_ORIENT_KINDS = frozenset({
    "amine", "diamine", "triamine", "tetraamine", "cycloalkane", "cycloamine",
})


def _typed_amine_orient(parent: dict, chain: list[int], subs: list) -> list[int] | None:
    atoms = _typed_group_atoms(parent, "amine")
    if not atoms or parent.get("kind") not in _TYPED_AMINE_ORIENT_KINDS:
        return None
    typed = {**parent, "amine_c_idx": atoms[0], "amine_c_idxs": atoms}
    if parent.get("kind") == "cycloalkane" and len(atoms) == 1:
        return _orient_ring_fixed(chain, typed, subs, "amine_c_idx")
    return _orient_amine(chain, typed, subs) if len(atoms) == 1 else _orient_diamine(chain, typed, subs)


def _typed_principal_orient(parent: dict, chain: list[int], subs: list) -> list[int] | None:
    print("using _typed_principal_orient ")
    for orient in (_typed_alcohol_orient, _typed_amine_orient):
        result = orient(parent, chain, subs)
        if result is not None:
            return result
    return None


def _orient_chain(parent: dict, substituents: list) -> list[int]:
    chain = list(parent.get("chain") or [])
    if not chain:
        return chain
    typed = _typed_principal_orient(parent, chain, substituents)
    print("_orient_chain: ",parent.get("kind"), chain, parent, substituents)
    if typed is not None:
        return typed
    if _typed_ring_acid(parent):
        return _orient_ring_fixed(chain, parent, substituents, "ring_attach_idx")

    return _orient_by_kind(parent.get("kind"), chain, parent, substituents)
