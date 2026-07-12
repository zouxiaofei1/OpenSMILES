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


def _chain_with_oh(info: dict) -> list[int]:
    mol: Mol = info["mol"]
    oh_c = info["hydroxyls"][0]["c_idx"]
    return _join_through(oh_c, _arms_from(mol, oh_c))


def _parent_dict(chain: list[int], kind: str, oh_c_idx: int | None) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, "oh_c_idx": oh_c_idx}


def _alcohol_parent(info: dict) -> dict:
    oh_c = info["hydroxyls"][0]["c_idx"]
    return _parent_dict(_chain_with_oh(info), "alcohol", oh_c)


def select_parent(info: dict) -> dict:
    if info.get("has_alcohol") and info.get("hydroxyls"):
        return _alcohol_parent(info)
    return _parent_dict(_longest_chain(info["mol"]), "alkane", None)
