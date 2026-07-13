from __future__ import annotations
from namepredict.layer3.substituent_extractor import alkyl_alpha_key
from namepredict.layer4.polyene import ene_locants, orient_polyene
def _pos_on(chain: list[int], c: int | None) -> int | None:
    if c is None or c not in chain:
        return None
    return chain.index(c) + 1
def _maybe_reverse(chain: list[int], pos: int) -> list[int]:
    if pos > (len(chain) + 1) // 2:
        return list(reversed(chain))
    return chain
def _locants_on(chain: list[int], substituents: list) -> list[int]:
    return sorted(chain.index(s["attach_idx"]) + 1 for s in substituents)
def _locant_key(locs: list[int]) -> tuple:
    return (locs, len(locs))
def _stem_loc_pairs(chain: list[int], substituents: list) -> list[tuple]:
    return sorted(
        (alkyl_alpha_key(s.get("en") or ""), chain.index(s["attach_idx"]) + 1)
        for s in substituents
    )
def _orient_key(chain: list[int], substituents: list) -> tuple:
    return (_locant_key(_locants_on(chain, substituents)), _stem_loc_pairs(chain, substituents))
def _prefer_chain(a: list[int], b: list[int], substituents: list) -> list[int]:
    return a if _orient_key(a, substituents) <= _orient_key(b, substituents) else b
def _orient_alkane(chain: list[int], substituents: list) -> list[int]:
    if not substituents or not chain:
        return chain
    return _prefer_chain(chain, list(reversed(chain)), substituents)
def _orient_by_single_fg(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    pos = _pos_on(chain, parent.get(key))
    if pos is None:
        return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    b, r = _pos_on(base, parent.get(key)), _pos_on(rev, parent.get(key))
    if b is not None and r is not None and b == r:
        return _prefer_chain(base, rev, substituents)
    return base
def _orient_alcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "oh_c_idx")
def _orient_thiol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "sh_c_idx")
def _orient_amine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "amine_c_idx")
def _orient_ketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_single_fg(chain, parent, substituents, "ketone_c_idx")
def _pair_locs_on(chain: list[int], cs) -> tuple[int, ...] | None:
    if not cs:
        return None
    locs = sorted(chain.index(c) + 1 for c in cs if c in chain)
    return tuple(locs) if len(locs) == len(cs) else None
def _better_pair_orient(a: list[int], b: list[int], cs, subs: list) -> list[int]:
    la, lb = _pair_locs_on(a, cs), _pair_locs_on(b, cs)
    if la is None:
        return b
    if lb is None or la < lb:
        return a
    if lb < la:
        return b
    return _prefer_chain(a, b, subs)
def _orient_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    if not cs:
        return chain
    return _better_pair_orient(chain, list(reversed(chain)), cs, subs)
def _orient_polyol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_diamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "amine_c_idxs", substituents)
def _orient_dione(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pair(chain, parent, "ketone_c_idxs", substituents)
def _orient_to_terminal(chain: list[int], c_idx: int | None) -> list[int]:
    pos = _pos_on(chain, c_idx)
    if pos is None or pos == 1:
        return chain
    return list(reversed(chain))
def _orient_acid(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("cooh_c_idx"))
def _orient_aldehyde(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("aldehyde_c_idx"))
def _orient_ester(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("ester_c_idx"))
def _orient_amide(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("amide_c_idx"))
def _orient_nitrile(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("nitrile_c_idx"))
def _orient_acyl_chloride(chain: list[int], parent: dict) -> list[int]:
    return _orient_to_terminal(chain, parent.get("acyl_c_idx"))
def _ene_ends_on(chain: list[int], ends: tuple[int, int] | None) -> tuple[int, int] | None:
    if not ends or ends[0] not in chain or ends[1] not in chain:
        return None
    return chain.index(ends[0]) + 1, chain.index(ends[1]) + 1
def _ene_locant_of(ends: tuple[int, int] | None) -> int | None:
    if ends is None:
        return None
    return min(ends)
def _tie_break_orient(
    base: list[int], ends0: tuple[int, int] | None, substituents: list
) -> list[int]:
    rev = list(reversed(base))
    if _ene_locant_of(_ene_ends_on(base, ends0)) == _ene_locant_of(
        _ene_ends_on(rev, ends0)
    ):
        return _prefer_chain(base, rev, substituents)
    return base
def _orient_by_bond(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    ends0 = parent.get(key)
    ends = _ene_ends_on(chain, ends0)
    if ends is None:
        return chain
    base = _maybe_reverse(chain, _ene_locant_of(ends) or 1)
    return _tie_break_orient(base, ends0, substituents)
def _orient_alkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "double_bond")
def _orient_cycloalkene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    ends = parent.get("double_bond")
    if not ends or ends[0] not in chain:
        return chain
    base = _rotate_to_front(chain, ends[0])
    # prefer orientation with second double-bond atom at locant 2
    if len(base) > 1 and base[1] != ends[1]:
        base = _rotate_to_front(list(reversed(chain)), ends[0])
    return base
def _orient_alkyne(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "triple_bond")
def _term_fn(fn):
    return lambda c, p, s: fn(c, p)
def _terminal_orienters() -> dict:
    return {
        "acid": _term_fn(_orient_acid), "alkenoic_acid": _term_fn(_orient_acid),
        "alkenal": _term_fn(_orient_aldehyde), "aldehyde": _term_fn(_orient_aldehyde),
        "ester": _term_fn(_orient_ester), "alkenoate": _term_fn(_orient_ester),
        "amide": _term_fn(_orient_amide), "nitrile": _term_fn(_orient_nitrile),
        "alkenenitrile": _term_fn(_orient_nitrile),
        "acyl_chloride": _term_fn(_orient_acyl_chloride),
    }
def _carbonyl_orienters() -> dict:
    return {
        **_terminal_orienters(),
        "ketone": _orient_ketone,
        "dione": _orient_dione,
    }
def _rotate_to_front(chain: list[int], atom: int) -> list[int]:
    if atom not in chain:
        return chain
    i = chain.index(atom)
    return chain[i:] + chain[:i]
def _rotations(chain: list[int]) -> list[list[int]]:
    return [chain[i:] + chain[:i] for i in range(len(chain))]
def _ring_candidates(chain: list[int]) -> list[list[int]]:
    out: list[list[int]] = []
    for base in (chain, list(reversed(chain))):
        out.extend(_rotations(base))
    return out
def _best_ring(chain: list[int], substituents: list) -> list[int]:
    best = chain
    for cand in _ring_candidates(chain):
        best = _prefer_chain(best, cand, substituents)
    return best
def _orient_cycloalkane(chain: list[int], parent: dict, substituents: list) -> list[int]:
    if not chain or not substituents:
        return chain
    return _best_ring(chain, substituents)
def _orient_ring_fixed(
    chain: list[int], parent: dict, substituents: list, key: str
) -> list[int]:
    c = parent.get(key)
    if c is None or c not in chain:
        return chain
    base = _rotate_to_front(chain, c)
    rev = _rotate_to_front(list(reversed(chain)), c)
    return _prefer_chain(base, rev, substituents)
def _orient_cycloalcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "oh_c_idx")
def _better_ring_pair(best, best_locs, cand, cs, subs):
    locs = _pair_locs_on(cand, cs)
    if locs is None:
        return best, best_locs
    if best_locs is None or locs < best_locs:
        return cand, locs
    if locs == best_locs:
        return _prefer_chain(best, cand, subs), best_locs
    return best, best_locs
def _orient_ring_pair(chain: list[int], parent: dict, key: str, subs: list) -> list[int]:
    cs = parent.get(key)
    if not cs or not chain:
        return chain
    best, best_locs = chain, _pair_locs_on(chain, cs)
    for cand in _ring_candidates(chain):
        best, best_locs = _better_ring_pair(best, best_locs, cand, cs, subs)
    return best
def _orient_benzenediol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "oh_c_idxs", substituents)
def _orient_diazine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_pair(chain, parent, "n_idxs", substituents)
def _orient_cycloketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_ring_fixed(chain, parent, substituents, "ketone_c_idx")
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
    """NH = 1; choose direction so the other N is at 3 (not 4)."""
    nh = parent.get("nh_idx")
    if nh is None or nh not in chain:
        return chain
    base = _rotate_to_front(chain, nh)
    rev = _rotate_to_front(list(reversed(chain)), nh)
    n = parent.get("n_idx")
    return base if _n_loc_on(base, n) <= _n_loc_on(rev, n) else rev
def _orient_pyridinecarboxylic(chain: list[int], parent: dict, substituents: list) -> list[int]:
    attach, virtual = parent.get("ring_attach_idx"), list(substituents)
    if attach is not None:
        virtual = virtual + [{"attach_idx": attach, "en": ""}]
    return _orient_pyridine(chain, parent, virtual)
def _orient_pyridin_fg(chain: list[int], parent: dict, substituents: list, key: str) -> list[int]:
    attach, virtual = parent.get(key), list(substituents)
    if attach is not None:
        virtual = virtual + [{"attach_idx": attach, "en": ""}]
    return _orient_pyridine(chain, parent, virtual)
def _orient_pyridinamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pyridin_fg(chain, parent, substituents, "amine_c_idx")
def _orient_pyridinol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_pyridin_fg(chain, parent, substituents, "oh_c_idx")
def _arene_orienters() -> dict:
    return {
        "pyridine": _orient_pyridine, "pyridinecarboxylic": _orient_pyridinecarboxylic,
        "pyridinamine": _orient_pyridinamine, "pyridinol": _orient_pyridinol,
        "furan": _orient_hetero5, "thiophene": _orient_hetero5, "pyrrole": _orient_hetero5,
        "imidazole": _orient_imidazole,
        "pyrimidine": _orient_diazine, "pyrazine": _orient_diazine,
        "pyridazine": _orient_diazine,
    }
def _hetero_orienters() -> dict:
    return {
        "alcohol": _orient_alcohol, "alkenol": _orient_alcohol, "thiol": _orient_thiol,
        "diol": _orient_polyol, "triol": _orient_polyol, "diamine": _orient_diamine,
        "cycloalcohol": _orient_cycloalcohol, "phenol": _orient_cycloalcohol,
        "benzenediol": _orient_benzenediol,
        "cycloamine": _orient_cycloamine, "aniline": _orient_cycloamine,
        "amine": _orient_amine, "sec_amine": _orient_amine, "tert_amine": _orient_amine,
        **_arene_orienters(),
    }
def _orient_polyene(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return orient_polyene(chain, parent, substituents, _prefer_chain)
def _unsat_orienters() -> dict:
    return {
        "cycloketone": _orient_cycloketone, "alkene": _orient_alkene,
        "alkenedioic": _orient_alkene, "polyene": _orient_polyene,
        "cycloalkene": _orient_cycloalkene, "alkyne": _orient_alkyne,
        "cycloalkane": _orient_cycloalkane, "benzene": _orient_cycloalkane,
        "benzoic": _orient_benzoic, "benzaldehyde": _orient_benzoic,
        "acetophenone": _orient_benzoic, "benzoate": _orient_benzoic,
        "benzonitrile": _orient_benzoic, "benzoyl_chloride": _orient_benzoic,
        "cycloalkanecarboxylic": _orient_benzoic,
    }
def _kind_orienters() -> dict:
    return {**_hetero_orienters(), **_unsat_orienters(), **_carbonyl_orienters()}
def _orient_by_kind(kind: str, chain: list[int], parent: dict, subs: list) -> list[int]:
    fn = _kind_orienters().get(kind)
    return _orient_alkane(chain, subs) if fn is None else fn(chain, parent, subs)
def _orient_chain(parent: dict, substituents: list) -> list[int]:
    chain = list(parent.get("chain") or [])
    if not chain:
        return chain
    return _orient_by_kind(parent.get("kind"), chain, parent, substituents)
def _fg_locant(oriented: dict, kinds: tuple, key: str) -> int | None:
    if oriented.get("kind") not in kinds:
        return None
    return _pos_on(oriented.get("chain") or [], oriented.get(key))
def _oh_locant(oriented: dict) -> int | None:
    return _fg_locant(
        oriented, ("alcohol", "alkenol", "cycloalcohol", "pyridinol"), "oh_c_idx",
    )
def _sh_locant(oriented: dict) -> int | None:
    return _fg_locant(oriented, ("thiol",), "sh_c_idx")
def _pair_locants(oriented: dict, kind: str, key: str) -> list[int] | None:
    if oriented.get("kind") != kind:
        return None
    locs = _pair_locs_on(oriented.get("chain") or [], oriented.get(key))
    return list(locs) if locs else None
def _oh_locants(oriented: dict) -> list[int] | None:
    kind = oriented.get("kind")
    if kind not in ("diol", "triol", "benzenediol"):
        return None
    locs = _pair_locs_on(oriented.get("chain") or [], oriented.get("oh_c_idxs"))
    return list(locs) if locs else None
def _amine_pair_locants(oriented: dict) -> list[int] | None:
    return _pair_locants(oriented, "diamine", "amine_c_idxs")
def _amine_locant(oriented: dict) -> int | None:
    return _fg_locant(
        oriented,
        ("amine", "cycloamine", "sec_amine", "tert_amine", "pyridinamine"),
        "amine_c_idx",
    )
def _ketone_locant(oriented: dict) -> int | None:
    return _fg_locant(oriented, ("ketone", "cycloketone"), "ketone_c_idx")
def _ketone_pair_locants(oriented: dict) -> list[int] | None:
    return _pair_locants(oriented, "dione", "ketone_c_idxs")
def _bond_locant(oriented: dict, kind: str, key: str) -> int | None:
    if oriented.get("kind") != kind:
        return None
    return _ene_locant_of(
        _ene_ends_on(oriented.get("chain") or [], oriented.get(key))
    )
def _ene_locant(oriented: dict) -> int | None:
    kind = oriented.get("kind")
    if kind in (
        "alkene", "alkenoic_acid", "alkenal", "alkenenitrile", "alkenol",
        "alkenoate", "cycloalkene", "alkenedioic",
    ):
        return _bond_locant(oriented, kind, "double_bond")
    return None
def _yne_locant(oriented: dict) -> int | None:
    return _bond_locant(oriented, "alkyne", "triple_bond")
def _omit_oh(oh_pos: int | None, n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloalcohol":
        return True
    if kind == "alkenol":
        return False
    return oh_pos == 1 and n_carbons <= 2
def _omit_sh(sh_pos: int | None, n_carbons: int) -> bool:
    return sh_pos == 1 and n_carbons <= 2
def _omit_amine(am_pos: int | None, n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloamine":
        return True
    return am_pos == 1 and n_carbons <= 2
def _omit_unsat(n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloalkene":
        return True
    if kind in (
        "alkenoic_acid", "alkenal", "alkenenitrile", "alkenol",
        "alkenoate", "alkenedioic",
    ):
        return False
    return n_carbons <= 3
def _with_locants(chain: list[int], substituents: list) -> list:
    out: list = []
    for s in substituents:
        loc = chain.index(s["attach_idx"]) + 1
        out.append({**s, "locant": loc})
    return out
def _unsat_locants(oriented: dict, n: int) -> dict:
    kind = oriented.get("kind")
    return {
        "ene_locant": _ene_locant(oriented),
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": _omit_unsat(n, kind),
        "yne_locant": _yne_locant(oriented),
        "omit_yne_locant": _omit_unsat(n),
    }
def _oh_am_locants(oriented: dict, n: int) -> dict:
    oh, am = _oh_locant(oriented), _amine_locant(oriented)
    kind = oriented.get("kind")
    return {
        "oh_locant": oh, "oh_locants": _oh_locants(oriented),
        "omit_oh_locant": _omit_oh(oh, n, kind),
        "amine_locant": am, "amine_locants": _amine_pair_locants(oriented),
        "omit_amine_locant": _omit_amine(am, n, kind),
    }
def _sh_locants(oriented: dict, n: int) -> dict:
    sh = _sh_locant(oriented)
    return {"sh_locant": sh, "omit_sh_locant": _omit_sh(sh, n)}
def _fg_locants(oriented: dict) -> dict:
    n = oriented.get("n_carbons", 0)
    base = {
        **_oh_am_locants(oriented, n), **_sh_locants(oriented, n),
        "ketone_locant": _ketone_locant(oriented),
        "ketone_locants": _ketone_pair_locants(oriented),
    }
    return {**base, **_unsat_locants(oriented, n)}
def _pack(oriented: dict, substituents: list) -> dict:
    base = {"parent": oriented, "substituents": substituents}
    return {**base, **_fg_locants(oriented)}
def number(parent: dict, substituents: list) -> dict:
    chain = _orient_chain(parent, substituents)
    oriented = {**parent, "chain": chain}
    return _pack(oriented, _with_locants(chain, substituents))
