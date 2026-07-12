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


def _orient_chain(parent: dict) -> list[int]:
    chain = list(parent.get("chain") or [])
    if parent.get("kind") != "alcohol" or not chain:
        return chain
    pos = _oh_position(parent)
    if pos is None:
        return chain
    return _maybe_reverse(chain, pos)


def _oh_locant(oriented: dict) -> int | None:
    chain = oriented.get("chain") or []
    oh_c = oriented.get("oh_c_idx")
    if oriented.get("kind") != "alcohol" or oh_c not in chain:
        return None
    return chain.index(oh_c) + 1


def number(parent: dict, substituents: list) -> dict:
    chain = _orient_chain(parent)
    oriented = {**parent, "chain": chain}
    oh_pos = _oh_locant(oriented)
    omit = oh_pos == 1 and oriented.get("n_carbons", 0) <= 2
    return {
        "parent": oriented,
        "substituents": substituents,
        "oh_locant": oh_pos,
        "omit_oh_locant": omit,
    }
