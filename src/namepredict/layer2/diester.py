"""Open-chain symmetric dialkyl alkane/alkenedioate parents (P-65.1.1 / P-31.1)."""
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


def _open_diester_ok(info: dict) -> bool:
    """Open-chain diester: no ring/alkyne/polyene; mono C=C allowed."""
    from namepredict.layer2.parent_core import _no_fgs

    if info.get("has_ring") or info.get("has_alkyne"):
        return False
    if len(info.get("double_bonds") or []) > 1:
        return False
    return _no_fgs(info, _DIESTER_BAD)


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


def _cover_atoms(info: dict, c_idxs: list[int]) -> list[int]:
    atoms = list(c_idxs)
    for db in info.get("double_bonds") or []:
        atoms.extend([db["c1"], db["c2"]])
    return atoms


def _cover_chain(info: dict, c_idxs: list[int]) -> list[int] | None:
    from namepredict.layer2.parent_core import _best_cover_pair

    if not _carbonyls_open(info["mol"], c_idxs):
        return None
    chain = _best_cover_pair(info["mol"], _cover_atoms(info, c_idxs))
    return chain if chain and set(c_idxs) <= set(chain) else None


def _ene_on_chain(chain: list[int], info: dict) -> tuple[int, int] | None:
    dbs = info.get("double_bonds") or []
    if len(dbs) != 1:
        return None
    c1, c2 = dbs[0]["c1"], dbs[0]["c2"]
    return (c1, c2) if c1 in chain and c2 in chain else None


def _diester_meta(info: dict, c_idxs: list[int], alkoxy_n: int, chain: list[int]) -> dict:
    meta = {
        "ester_c_idxs": list(c_idxs), "alkoxy_n": alkoxy_n,
        "symmetric": True, "mol": info["mol"],
    }
    ene = _ene_on_chain(chain, info)
    if ene is not None:
        meta["double_bond"] = ene
    return meta


def _build_parent(info: dict, esters: list, alkoxy_n: int) -> dict | None:
    from namepredict.layer2.parent_core import _parent_dict

    c_idxs = [esters[0]["c_idx"], esters[1]["c_idx"]]
    chain = _cover_chain(info, c_idxs)
    if not chain:
        return None
    if info.get("has_alkene") and _ene_on_chain(chain, info) is None:
        return None
    return _parent_dict(chain, "diester", **_diester_meta(info, c_idxs, alkoxy_n, chain))


def _diester_parent(info: dict) -> dict | None:
    """Symmetric open-chain diester: chain through both ester carbonyls."""
    if not _open_diester_ok(info):
        return None
    esters = _two_esters(info)
    if esters is None:
        return None
    n = _sym_alkoxy(info, esters)
    return _build_parent(info, esters, n) if n is not None else None
