from __future__ import annotations


def _oh_position(parent: dict) -> int | None:
    chain = parent.get("chain") or []
    oh_c = parent.get("oh_c_idx")
    if oh_c is None or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


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


def _orient_chain(parent: dict, substituents: list) -> list[int]:
    chain = list(parent.get("chain") or [])
    if not chain:
        return chain
    if parent.get("kind") == "alcohol":
        return _orient_alcohol(chain, parent, substituents)
    if parent.get("kind") == "acid":
        return _orient_acid(chain, parent)
    return _orient_alkane(chain, substituents)


def _oh_locant(oriented: dict) -> int | None:
    chain = oriented.get("chain") or []
    oh_c = oriented.get("oh_c_idx")
    if oriented.get("kind") != "alcohol" or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


def _omit_oh(oh_pos: int | None, n_carbons: int) -> bool:
    return oh_pos == 1 and n_carbons <= 2


def _with_locants(chain: list[int], substituents: list) -> list:
    out: list = []
    for s in substituents:
        loc = chain.index(s["attach_idx"]) + 1
        out.append({**s, "locant": loc})
    return out


def _pack(oriented: dict, substituents: list, oh_pos: int | None) -> dict:
    return {
        "parent": oriented,
        "substituents": substituents,
        "oh_locant": oh_pos,
        "omit_oh_locant": _omit_oh(oh_pos, oriented.get("n_carbons", 0)),
    }


def number(parent: dict, substituents: list) -> dict:
    chain = _orient_chain(parent, substituents)
    oriented = {**parent, "chain": chain}
    numbered = _with_locants(chain, substituents)
    return _pack(oriented, numbered, _oh_locant(oriented))
