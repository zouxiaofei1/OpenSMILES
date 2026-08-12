"""Orient sat_hetero* chains: plain, carboxylic, lactone/lactam, repl (L4)."""
from __future__ import annotations

from namepredict.layer4._chain_orient import (
    _pair_locants,
    _rotate_to,
    _stem_loc_pairs,
    _sub_locants,
)
from namepredict.layer4.locants.engine import choose_numbering


def _attach_idxs(substituents: list) -> list[int]:
    return [s["attach_idx"] for s in substituents if s.get("attach_idx") is not None]


def orient_sat_hetero_repl(
    chain: list[int], parent: dict, substituents: list,
) -> list[int]:
    """multi_hetero: hetero set + a-order + sub set."""
    hs = list(parent.get("hetero_idxs") or [])
    plan = choose_numbering(
        chain, "multi_hetero", hetero_atoms=hs, hetero_z=parent.get("hetero_z"),
        sub_attach=_attach_idxs(substituents), scaffold_id="sat_hetero_repl",
    )
    return list(plan.atom_order) if plan.atom_order else chain


SAT_HETERO_PLAIN = (
    "aziridine", "oxirane", "oxolane", "oxane", "pyrrolidine", "piperidine",
    "morpholine", "piperazine", "dioxolane", "dioxane", "thiolane",
)


def _orient_key(chain: list[int], subs: list) -> tuple:
    """Locant set first, then stem-alpha pairs when sets tie (P-14.5)."""
    return (tuple(_sub_locants(chain, subs)), _stem_loc_pairs(chain, subs))


def _attach_loc(chain: list[int], attach: int | None) -> int:
    if attach is None or attach not in chain:
        return 99
    return chain.index(attach) + 1


def _attach_key(chain: list[int], subs: list, attach: int | None) -> tuple:
    """Principal attach, then loc set + stem pairs (COOH/one lowest)."""
    return (_attach_loc(chain, attach), _orient_key(chain, subs))


def _prefer(a: list[int], b: list[int], subs: list, attach: int | None = None) -> list[int]:
    if attach is not None:
        return a if _attach_key(a, subs, attach) <= _attach_key(b, subs, attach) else b
    return a if _orient_key(a, subs) <= _orient_key(b, subs) else b


def _fixed(
    chain: list[int], atom: int | None, subs: list, attach: int | None = None,
) -> list[int]:
    if atom is None or atom not in chain:
        return chain
    base, rev = _rotate_to(chain, atom), _rotate_to(list(reversed(chain)), atom)
    return _prefer(base, rev, subs, attach)


def _ring_cands(chain: list[int]) -> list[list[int]]:
    n, out = len(chain), []
    for i in range(n):
        rot = chain[i:] + chain[:i]
        out.append(rot)
        out.append(list(reversed(rot)))
    return out


def _orient_pair(chain: list[int], cs, subs: list) -> list[int]:
    best, best_locs = chain, _pair_locants(chain, cs)
    for cand in _ring_cands(chain):
        locs = _pair_locants(cand, cs)
        if locs is None:
            continue
        if best_locs is None or locs < best_locs:
            best, best_locs = cand, locs
        elif locs == best_locs:
            best = _prefer(best, cand, subs)
    return best


def orient_sat_hetero(
    chain: list[int], parent: dict, substituents: list,
) -> list[int]:
    """Mono hetero = 1; di-hetero use pair locs (O before N for morpholine)."""
    hs = parent.get("hetero_idxs") or []
    if len(hs) >= 2:
        return _orient_pair(chain, hs, substituents)
    return _fixed(chain, parent.get("hetero_idx"), substituents)


def sat_hetero_orienters() -> dict:
    """kind → orient-fn for plain / repl sat_hetero families."""
    d = {k: orient_sat_hetero for k in SAT_HETERO_PLAIN}
    return {**d, "sat_hetero_repl": orient_sat_hetero_repl}
