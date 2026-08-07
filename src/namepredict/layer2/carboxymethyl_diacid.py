"""P-65.1.2.2.3 terminal carboxyalkyl-substituted diacids."""
from itertools import combinations
from namepredict.tools.side_facts import carboxyalkyl_arms

def _acids(info): return {int(x["c_idx"]) for x in info.get("carboxyls") or []}
def _graph(info):
    return {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 6}
            for a in info["mol"].GetAtoms() if a.GetAtomicNum() == 6}
def _path(graph, start, end):
    todo = [(start, [start])]
    while todo:
        atom, path = todo.pop()
        if atom == end: return path
        todo += [(n, path + [n]) for n in graph[atom] if n not in path]
def _paths(info):
    graph, acids = _graph(info), sorted(_acids(info))
    return [p for a, b in combinations(acids, 2) if (p := _path(graph, a, b))]
def _arms(info, chain):
    return carboxyalkyl_arms(info["mol"], chain, sorted(_acids(info) - {chain[0], chain[-1]}))
def _claims_all(info, chain):
    arms = _arms(info, chain)
    claimed = set(chain) | {x for arm in arms for x in (*arm.atoms, arm.carboxyl)}
    return len(arms) == len(_acids(info)) - 2 and claimed == set(_graph(info))
def _score(info, chain):
    arms = _arms(info, chain)
    return min(tuple(sorted(chain.index(a.attachment) + 1 for a in arms)), tuple(sorted(len(chain)-chain.index(a.attachment) for a in arms)))
def _main(info):
    paths = [p for p in _paths(info) if len(p) >= 7 and _claims_all(info, p)]
    return min(paths, key=lambda p: (len(p), _score(info, p)), default=None)
def _ok(info):
    bad = ("has_ring", "has_alkyne", "has_ester", "has_amide", "has_ketone", "has_aldehyde", "has_nitrile", "has_alcohol", "has_amine", "has_halo", "has_thiol")
    return not any(info.get(k) for k in bad) and not any(x.get("anion") for x in info.get("carboxyls") or [])
def _eligible_unsat(info, chain):
    bonds = info.get("double_bonds") or []
    return len(bonds) <= 1 and (not bonds or {bonds[0]["c1"], bonds[0]["c2"]} <= set(chain))
def _double(info, chain):
    return next(((x["c1"], x["c2"]) for x in info.get("double_bonds") or [] if {x["c1"], x["c2"]} <= set(chain)), None)
def carboxymethyl_diacid_parent(info):
    chain = _main(info)
    if not chain or len(chain) < 5 or len(_acids(info)) < 3 or not _ok(info): return None
    arms = _arms(info, chain)
    if not _claims_all(info, chain) or not _eligible_unsat(info, chain): return None
    return {"chain": chain, "n_carbons": len(chain), "kind": "diacid", "cooh_c_idxs": [chain[0],chain[-1]], "anion": False, "double_bond": _double(info, chain), "carboxymethyl_arms": arms}
def is_carboxymethyl_diacid(info): return carboxymethyl_diacid_parent(info) is not None
