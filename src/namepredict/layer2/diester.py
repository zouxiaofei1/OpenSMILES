"""Open-chain symmetric dialkyl alkanedioate parents (P-65.1.1 / P-65.6)."""
from __future__ import annotations


_DIESTER_BAD = (
    "has_acid", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_amine", "has_alcohol",
    "has_anhydride", "has_thiol", "has_phosphate", "has_phosphonic",
    "has_carbamate",
)


def _two_esters(info: dict) -> list | None:
    xs = info.get("esters") or []
    return xs if len(xs) == 2 else None


def _open_sat_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs

    ring_unsat = info.get("has_ring") or info.get("has_alkene") or info.get("has_alkyne")
    return (not ring_unsat) and _no_fgs(info, _DIESTER_BAD)


def _side_of(info: dict, e: dict) -> dict:
    from namepredict.layer2.alkoxy_side import classify_alkoxy

    return classify_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"])


def _heavy_beyond(mol, i: int, seen: set[int]) -> list[int]:
    atom = mol.GetAtomWithIdx(i)
    return [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetIdx() not in seen and n.GetAtomicNum() != 1
    ]


def _push_carbon(mol, i: int, seen: set[int], cs: list[int], stack: list[int]) -> bool:
    """Add carbon i to arm; False if non-carbon heavy."""
    if mol.GetAtomWithIdx(i).GetAtomicNum() != 6:
        return False
    seen.add(i)
    cs.append(i)
    stack.extend(_heavy_beyond(mol, i, seen))
    return True


def _arm_carbons(mol, o_idx: int, start: int) -> list[int] | None:
    """Carbons of O–R arm; None if hetero / non-carbon heavy atoms."""
    seen, stack, cs = {o_idx}, [start], []
    while stack:
        i = stack.pop()
        if i in seen:
            continue
        if not _push_carbon(mol, i, seen, cs, stack):
            return None
    return cs


def _c_degree(mol, i: int) -> int:
    return sum(1 for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 6)


def _is_unbranched(mol, cs: list[int]) -> bool:
    """True iff every arm carbon has ≤2 carbon neighbors (n-alkyl)."""
    return all(_c_degree(mol, i) <= 2 for i in cs)


def _true_linear_n(mol, o_idx: int, alkoxy_c: int) -> int | None:
    cs = _arm_carbons(mol, o_idx, alkoxy_c)
    return len(cs) if cs and _is_unbranched(mol, cs) else None


def _linear_n(info: dict, e: dict) -> int | None:
    side = _side_of(info, e)
    if side.get("alkoxy_en") or side.get("alkoxy_n") is None:
        return None
    return _true_linear_n(info["mol"], e["o_idx"], e["alkoxy_c_idx"])


def _sym_alkoxy(info: dict, esters: list) -> int | None:
    n0, n1 = _linear_n(info, esters[0]), _linear_n(info, esters[1])
    return n0 if n0 is not None and n0 == n1 else None


def _carbonyls_open(mol, c_idxs: list[int]) -> bool:
    return all(not mol.GetAtomWithIdx(int(c)).IsInRing() for c in c_idxs)


def _cover_chain(info: dict, c_idxs: list[int]) -> list[int] | None:
    from namepredict.layer2.parent_selector import _best_cover_pair

    if not _carbonyls_open(info["mol"], c_idxs):
        return None
    chain = _best_cover_pair(info["mol"], c_idxs)
    return chain if chain and set(c_idxs) <= set(chain) else None


def _diester_meta(c_idxs: list[int], alkoxy_n: int) -> dict:
    return {"ester_c_idxs": list(c_idxs), "alkoxy_n": alkoxy_n, "symmetric": True}


def _build_parent(info: dict, esters: list, alkoxy_n: int) -> dict | None:
    from namepredict.layer2.parent_selector import _parent_dict

    c_idxs = [esters[0]["c_idx"], esters[1]["c_idx"]]
    chain = _cover_chain(info, c_idxs)
    if not chain:
        return None
    return _parent_dict(chain, "diester", **_diester_meta(c_idxs, alkoxy_n))


def _diester_parent(info: dict) -> dict | None:
    """Symmetric open-chain diester: chain through both ester carbonyls."""
    if not _open_sat_ok(info):
        return None
    esters = _two_esters(info)
    if esters is None:
        return None
    n = _sym_alkoxy(info, esters)
    return _build_parent(info, esters, n) if n is not None else None
