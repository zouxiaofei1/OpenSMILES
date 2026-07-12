from __future__ import annotations


def _oh_position(parent: dict) -> int | None:
    chain = parent.get("chain") or []
    oh_c = parent.get("oh_c_idx")
    if oh_c is None or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


def _amine_position(parent: dict) -> int | None:
    chain = parent.get("chain") or []
    am_c = parent.get("amine_c_idx")
    if am_c is None or am_c not in chain:
        return None
    return chain.index(am_c) + 1


def _maybe_reverse(chain: list[int], pos: int) -> list[int]:
    if pos > (len(chain) + 1) // 2:
        return list(reversed(chain))
    return chain


def _locants_on(chain: list[int], substituents: list) -> list[int]:
    return sorted(chain.index(s["attach_idx"]) + 1 for s in substituents)


def _locant_key(locs: list[int]) -> tuple:
    return (locs, len(locs))


def _prefer_chain(a: list[int], b: list[int], substituents: list) -> list[int]:
    ka = _locant_key(_locants_on(a, substituents))
    kb = _locant_key(_locants_on(b, substituents))
    return a if ka <= kb else b


def _orient_alkane(chain: list[int], substituents: list) -> list[int]:
    if not substituents or not chain:
        return chain
    return _prefer_chain(chain, list(reversed(chain)), substituents)


def _oh_pos_on(chain: list[int], oh_c: int | None) -> int | None:
    if oh_c is None or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


def _orient_alcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    pos = _oh_position(parent)
    if pos is None:
        return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    oh_b = _oh_pos_on(base, parent.get("oh_c_idx"))
    oh_r = _oh_pos_on(rev, parent.get("oh_c_idx"))
    if oh_b is not None and oh_r is not None and oh_b == oh_r:
        return _prefer_chain(base, rev, substituents)
    return base


def _amine_pos_on(chain: list[int], am_c: int | None) -> int | None:
    if am_c is None or am_c not in chain:
        return None
    return chain.index(am_c) + 1


def _orient_amine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    pos = _amine_position(parent)
    if pos is None:
        return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    am_b = _amine_pos_on(base, parent.get("amine_c_idx"))
    am_r = _amine_pos_on(rev, parent.get("amine_c_idx"))
    if am_b is not None and am_r is not None and am_b == am_r:
        return _prefer_chain(base, rev, substituents)
    return base


def _cooh_pos_on(chain: list[int], cooh_c: int | None) -> int | None:
    if cooh_c is None or cooh_c not in chain:
        return None
    return chain.index(cooh_c) + 1


def _orient_acid(chain: list[int], parent: dict) -> list[int]:
    pos = _cooh_pos_on(chain, parent.get("cooh_c_idx"))
    if pos is None:
        return chain
    if pos == 1:
        return chain
    return list(reversed(chain))


def _ald_pos_on(chain: list[int], ald_c: int | None) -> int | None:
    if ald_c is None or ald_c not in chain:
        return None
    return chain.index(ald_c) + 1


def _orient_aldehyde(chain: list[int], parent: dict) -> list[int]:
    pos = _ald_pos_on(chain, parent.get("aldehyde_c_idx"))
    if pos is None:
        return chain
    if pos == 1:
        return chain
    return list(reversed(chain))


def _ester_pos_on(chain: list[int], ester_c: int | None) -> int | None:
    if ester_c is None or ester_c not in chain:
        return None
    return chain.index(ester_c) + 1


def _orient_ester(chain: list[int], parent: dict) -> list[int]:
    pos = _ester_pos_on(chain, parent.get("ester_c_idx"))
    if pos is None:
        return chain
    if pos == 1:
        return chain
    return list(reversed(chain))


def _amide_pos_on(chain: list[int], amide_c: int | None) -> int | None:
    if amide_c is None or amide_c not in chain:
        return None
    return chain.index(amide_c) + 1


def _orient_amide(chain: list[int], parent: dict) -> list[int]:
    pos = _amide_pos_on(chain, parent.get("amide_c_idx"))
    if pos is None:
        return chain
    if pos == 1:
        return chain
    return list(reversed(chain))


def _nitrile_pos_on(chain: list[int], nit_c: int | None) -> int | None:
    if nit_c is None or nit_c not in chain:
        return None
    return chain.index(nit_c) + 1


def _orient_nitrile(chain: list[int], parent: dict) -> list[int]:
    pos = _nitrile_pos_on(chain, parent.get("nitrile_c_idx"))
    if pos is None:
        return chain
    if pos == 1:
        return chain
    return list(reversed(chain))


def _ketone_pos_on(chain: list[int], ket_c: int | None) -> int | None:
    if ket_c is None or ket_c not in chain:
        return None
    return chain.index(ket_c) + 1


def _ketone_position(parent: dict) -> int | None:
    return _ketone_pos_on(parent.get("chain") or [], parent.get("ketone_c_idx"))


def _orient_ketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    pos = _ketone_position(parent)
    if pos is None:
        return chain
    base = _maybe_reverse(chain, pos)
    rev = list(reversed(base))
    k_b = _ketone_pos_on(base, parent.get("ketone_c_idx"))
    k_r = _ketone_pos_on(rev, parent.get("ketone_c_idx"))
    if k_b is not None and k_r is not None and k_b == k_r:
        return _prefer_chain(base, rev, substituents)
    return base


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


def _orient_alkyne(chain: list[int], parent: dict, substituents: list) -> list[int]:
    return _orient_by_bond(chain, parent, substituents, "triple_bond")


def _terminal_orienters() -> dict:
    return {
        "acid": lambda c, p, s: _orient_acid(c, p),
        "aldehyde": lambda c, p, s: _orient_aldehyde(c, p),
        "ester": lambda c, p, s: _orient_ester(c, p),
        "amide": lambda c, p, s: _orient_amide(c, p),
        "nitrile": lambda c, p, s: _orient_nitrile(c, p),
    }


def _carbonyl_orienters() -> dict:
    return {**_terminal_orienters(), "ketone": _orient_ketone}


def _rotate_to_front(chain: list[int], atom: int) -> list[int]:
    if atom not in chain:
        return chain
    i = chain.index(atom)
    return chain[i:] + chain[:i]


def _orient_cycloalcohol(chain: list[int], parent: dict, substituents: list) -> list[int]:
    oh_c = parent.get("oh_c_idx")
    if oh_c is None or oh_c not in chain:
        return chain
    base = _rotate_to_front(chain, oh_c)
    rev = _rotate_to_front(list(reversed(chain)), oh_c)
    return _prefer_chain(base, rev, substituents)


def _orient_cycloketone(chain: list[int], parent: dict, substituents: list) -> list[int]:
    ket_c = parent.get("ketone_c_idx")
    if ket_c is None or ket_c not in chain:
        return chain
    base = _rotate_to_front(chain, ket_c)
    rev = _rotate_to_front(list(reversed(chain)), ket_c)
    return _prefer_chain(base, rev, substituents)


def _orient_cycloamine(chain: list[int], parent: dict, substituents: list) -> list[int]:
    am_c = parent.get("amine_c_idx")
    if am_c is None or am_c not in chain:
        return chain
    base = _rotate_to_front(chain, am_c)
    rev = _rotate_to_front(list(reversed(chain)), am_c)
    return _prefer_chain(base, rev, substituents)


def _kind_orienters() -> dict:
    base = {
        "alcohol": _orient_alcohol,
        "cycloalcohol": _orient_cycloalcohol,
        "cycloketone": _orient_cycloketone,
        "cycloamine": _orient_cycloamine,
        "amine": _orient_amine,
        "alkene": _orient_alkene,
        "alkyne": _orient_alkyne,
    }
    return {**base, **_carbonyl_orienters()}


def _orient_by_kind(kind: str, chain: list[int], parent: dict, subs: list) -> list[int]:
    fn = _kind_orienters().get(kind)
    if fn is None:
        return _orient_alkane(chain, subs)
    return fn(chain, parent, subs)


def _orient_chain(parent: dict, substituents: list) -> list[int]:
    chain = list(parent.get("chain") or [])
    if not chain:
        return chain
    return _orient_by_kind(parent.get("kind"), chain, parent, substituents)


def _oh_locant(oriented: dict) -> int | None:
    chain = oriented.get("chain") or []
    oh_c = oriented.get("oh_c_idx")
    if oriented.get("kind") not in ("alcohol", "cycloalcohol") or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


def _amine_locant(oriented: dict) -> int | None:
    chain = oriented.get("chain") or []
    am_c = oriented.get("amine_c_idx")
    if oriented.get("kind") not in ("amine", "cycloamine") or am_c not in chain:
        return None
    return chain.index(am_c) + 1


def _ketone_locant(oriented: dict) -> int | None:
    chain = oriented.get("chain") or []
    ket_c = oriented.get("ketone_c_idx")
    if oriented.get("kind") not in ("ketone", "cycloketone") or ket_c not in chain:
        return None
    return chain.index(ket_c) + 1


def _bond_locant(oriented: dict, kind: str, key: str) -> int | None:
    if oriented.get("kind") != kind:
        return None
    return _ene_locant_of(
        _ene_ends_on(oriented.get("chain") or [], oriented.get(key))
    )


def _ene_locant(oriented: dict) -> int | None:
    return _bond_locant(oriented, "alkene", "double_bond")


def _yne_locant(oriented: dict) -> int | None:
    return _bond_locant(oriented, "alkyne", "triple_bond")


def _omit_oh(oh_pos: int | None, n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloalcohol":
        return True
    return oh_pos == 1 and n_carbons <= 2


def _omit_amine(am_pos: int | None, n_carbons: int, kind: str | None = None) -> bool:
    if kind == "cycloamine":
        return True
    return am_pos == 1 and n_carbons <= 2


def _omit_unsat(n_carbons: int) -> bool:
    return n_carbons <= 3


def _with_locants(chain: list[int], substituents: list) -> list:
    out: list = []
    for s in substituents:
        loc = chain.index(s["attach_idx"]) + 1
        out.append({**s, "locant": loc})
    return out


def _unsat_locants(oriented: dict, n: int) -> dict:
    return {
        "ene_locant": _ene_locant(oriented),
        "omit_ene_locant": _omit_unsat(n),
        "yne_locant": _yne_locant(oriented),
        "omit_yne_locant": _omit_unsat(n),
    }


def _oh_am_locants(oriented: dict, n: int) -> dict:
    oh, am = _oh_locant(oriented), _amine_locant(oriented)
    kind = oriented.get("kind")
    return {
        "oh_locant": oh,
        "omit_oh_locant": _omit_oh(oh, n, kind),
        "amine_locant": am,
        "omit_amine_locant": _omit_amine(am, n, kind),
    }


def _fg_locants(oriented: dict) -> dict:
    n = oriented.get("n_carbons", 0)
    base = {**_oh_am_locants(oriented, n), "ketone_locant": _ketone_locant(oriented)}
    return {**base, **_unsat_locants(oriented, n)}


def _pack(oriented: dict, substituents: list) -> dict:
    base = {"parent": oriented, "substituents": substituents}
    return {**base, **_fg_locants(oriented)}


def number(parent: dict, substituents: list) -> dict:
    chain = _orient_chain(parent, substituents)
    oriented = {**parent, "chain": chain}
    return _pack(oriented, _with_locants(chain, substituents))
