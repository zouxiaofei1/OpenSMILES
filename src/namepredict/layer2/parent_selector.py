from __future__ import annotations

from rdkit.Chem import Mol


def _carbon_neighbors(mol: Mol, idx: int) -> list[int]:
    atom = mol.GetAtomWithIdx(idx)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]


def _extend_best(mol: Mol, node: int, path: list[int], forbid: set[int], best: list[int]) -> list[int]:
    for nb in _carbon_neighbors(mol, node):
        if nb in path or nb in forbid:
            continue
        cand = _dfs_path(mol, nb, path + [nb], forbid)
        if len(cand) > len(best):
            best = cand
    return best


def _dfs_path(mol: Mol, node: int, path: list[int], forbid: set[int]) -> list[int]:
    return _extend_best(mol, node, path, forbid, path)


def _longest_from(mol: Mol, start: int, forbidden: set[int] | None = None) -> list[int]:
    return _dfs_path(mol, start, [start], forbidden or set())


def _all_carbons(mol: Mol) -> list[int]:
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]


def _side_count(mol: Mol, chain: list[int]) -> int:
    chain_set = set(chain)
    n = 0
    for c in chain:
        atom = mol.GetAtomWithIdx(c)
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() != 1 and nb.GetIdx() not in chain_set:
                n += 1
    return n


def _chain_key(mol: Mol, path: list[int]) -> tuple:
    return (len(path), _side_count(mol, path))


def _better(mol: Mol, cand: list[int], best: list[int]) -> bool:
    if not best:
        return True
    return _chain_key(mol, cand) > _chain_key(mol, best)


def _best_among(mol: Mol, seeds: list[int]) -> list[int]:
    best: list[int] = []
    for c in seeds:
        path = _longest_from(mol, c)
        if _better(mol, path, best):
            best = path
    return best


def _longest_chain(mol: Mol, seeds: list[int] | None = None) -> list[int]:
    return _best_among(mol, seeds or _all_carbons(mol))


def _arms_from(mol: Mol, center: int) -> list[list[int]]:
    forbid = {center}
    return [_longest_from(mol, nb, forbid) for nb in _carbon_neighbors(mol, center)]


def _join_through(center: int, arms: list[list[int]]) -> list[int]:
    arms = sorted(arms, key=len, reverse=True)
    if not arms:
        return [center]
    if len(arms) == 1:
        return list(reversed(arms[0])) + [center]
    return list(reversed(arms[0])) + [center] + arms[1]


def _chain_through(info: dict, c_idx: int) -> list[int]:
    mol: Mol = info["mol"]
    return _join_through(c_idx, _arms_from(mol, c_idx))


def _parent_core(chain: list[int], kind: str) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind}


def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    return {**_parent_core(chain, kind), **kw}


def _alcohol_parent(info: dict) -> dict:
    oh_c = info["hydroxyls"][0]["c_idx"]
    return _parent_dict(_chain_through(info, oh_c), "alcohol", oh_c_idx=oh_c)


def _acid_parent(info: dict) -> dict:
    cooh_c = info["carboxyls"][0]["c_idx"]
    return _parent_dict(_chain_through(info, cooh_c), "acid", cooh_c_idx=cooh_c)


def _ketone_parent(info: dict) -> dict:
    ket_c = info["ketones"][0]["c_idx"]
    return _parent_dict(_chain_through(info, ket_c), "ketone", ketone_c_idx=ket_c)


def _aldehyde_parent(info: dict) -> dict:
    ald_c = info["aldehydes"][0]["c_idx"]
    return _parent_dict(_chain_through(info, ald_c), "aldehyde", aldehyde_c_idx=ald_c)


def _ester_parent(info: dict) -> dict:
    e = info["esters"][0]
    chain = _chain_through(info, e["c_idx"])
    alkoxy_n = len(_longest_from(info["mol"], e["alkoxy_c_idx"]))
    return _parent_dict(
        chain, "ester",
        ester_c_idx=e["c_idx"], o_idx=e["o_idx"], alkoxy_c_idx=e["alkoxy_c_idx"],
        alkoxy_n=alkoxy_n,
    )


def _is_mono_ester(info: dict) -> bool:
    esters = info.get("esters") or []
    return bool(info.get("has_ester")) and len(esters) == 1


def _best_arm_away(mol: Mol, from_c: int, forbid: set[int]) -> list[int]:
    best: list[int] = []
    for nb in _carbon_neighbors(mol, from_c):
        if nb in forbid:
            continue
        path = _longest_from(mol, nb, forbid | {from_c})
        if len(path) > len(best):
            best = path
    return best


def _chain_through_bond(mol: Mol, c1: int, c2: int) -> list[int]:
    left = _best_arm_away(mol, c1, {c2})
    right = _best_arm_away(mol, c2, {c1})
    return list(reversed(left)) + [c1, c2] + right


def _alkene_parent(info: dict) -> dict:
    db = info["double_bonds"][0]
    c1, c2 = db["c1"], db["c2"]
    chain = _chain_through_bond(info["mol"], c1, c2)
    return _parent_dict(chain, "alkene", double_bond=(c1, c2))


def _alkyne_parent(info: dict) -> dict:
    tb = info["triple_bonds"][0]
    c1, c2 = tb["c1"], tb["c2"]
    chain = _chain_through_bond(info["mol"], c1, c2)
    return _parent_dict(chain, "alkyne", triple_bond=(c1, c2))


def _is_mono_alkene(info: dict) -> bool:
    bonds = info.get("double_bonds") or []
    return bool(info.get("has_alkene")) and len(bonds) == 1


def _is_mono_alkyne(info: dict) -> bool:
    triples = info.get("triple_bonds") or []
    doubles = info.get("double_bonds") or []
    return len(triples) == 1 and len(doubles) == 0


def _carbonyl_parent(info: dict) -> dict | None:
    if info.get("has_acid") and info.get("carboxyls"):
        return _acid_parent(info)
    if _is_mono_ester(info):
        return _ester_parent(info)
    if info.get("has_aldehyde") and info.get("aldehydes"):
        return _aldehyde_parent(info)
    if info.get("has_ketone") and info.get("ketones"):
        return _ketone_parent(info)
    return None


def _fg_parent(info: dict) -> dict | None:
    carb = _carbonyl_parent(info)
    if carb is not None:
        return carb
    if info.get("has_alcohol") and info.get("hydroxyls"):
        return _alcohol_parent(info)
    return None


def select_parent(info: dict) -> dict:
    fg = _fg_parent(info)
    if fg is not None:
        return fg
    if _is_mono_alkyne(info):
        return _alkyne_parent(info)
    if _is_mono_alkene(info):
        return _alkene_parent(info)
    return _parent_dict(_longest_chain(info["mol"]), "alkane")

