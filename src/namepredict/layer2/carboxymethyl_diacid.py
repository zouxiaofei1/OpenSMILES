"""L2 P-65.1.2.2.3 parent selection for carboxymethyl-substituted diacids."""
from __future__ import annotations

from itertools import combinations

from namepredict.layer2.side_facts import carboxymethyl_arms


def _acids(info: dict) -> set[int]:
    return {int(x["c_idx"]) for x in info.get("carboxyls") or []}


def _carbon_graph(info: dict) -> dict[int, set[int]]:
    mol = info["mol"]
    return {a.GetIdx(): {n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 6}
            for a in mol.GetAtoms() if a.GetAtomicNum() == 6}


def _path(graph: dict[int, set[int]], start: int, end: int) -> list[int] | None:
    todo = [(start, [start])]
    while todo:
        atom, path = todo.pop()
        if atom == end:
            return path
        todo.extend((n, path + [n]) for n in graph[atom] if n not in path)
    return None


def _acid_paths(info: dict) -> list[list[int]]:
    graph, acids = _carbon_graph(info), sorted(_acids(info))
    return [p for a, b in combinations(acids, 2) if (p := _path(graph, a, b))]


def _score_path(info: dict, chain: list[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Orient and score a path by its lowest carboxymethyl locant set."""
    graph, acids = _carbon_graph(info), _remaining_acids(info, chain)
    choices = []
    for oriented in (chain, list(reversed(chain))):
        parent = set(oriented)
        attachments = [next(iter(graph[_methylene_on_chain(graph, acid, parent)] & parent)) for acid in acids]
        locants = tuple(sorted(oriented.index(attachment) + 1 for attachment in attachments))
        choices.append((locants, tuple(oriented)))
    return min(choices)


def _main_path(info: dict) -> list[int] | None:
    paths = [path for path in _acid_paths(info) if _all_carboxymethyl_arms(info, path)]
    if not paths:
        return None
    longest = max(map(len, paths))
    return list(min(_score_path(info, path) for path in paths if len(path) == longest)[1])


def _is_open_saturated(info: dict) -> bool:
    return not (info.get("has_ring") or info.get("has_alkene") or info.get("has_alkyne"))


def _only_acid_fgs(info: dict) -> bool:
    bad = ("has_ester", "has_amide", "has_ketone", "has_aldehyde", "has_nitrile",
           "has_alcohol", "has_amine", "has_halo", "has_thiol")
    return not any(info.get(key) for key in bad)


def _remaining_acids(info: dict, chain: list[int]) -> set[int]:
    return _acids(info) - {chain[0], chain[-1]}


def _methylene_on_chain(graph: dict[int, set[int]], acid: int, chain: set[int]) -> int | None:
    arms = graph[acid] - chain
    return next(iter(arms)) if len(arms) == 1 else None


def _is_methylene_arm(graph: dict[int, set[int]], acid: int, chain: set[int]) -> bool:
    arm = _methylene_on_chain(graph, acid, chain)
    return arm is not None and len(graph[arm]) == 2 and bool((graph[arm] - {acid}) & chain)


def _all_carboxymethyl_arms(info: dict, chain: list[int]) -> bool:
    graph, main = _carbon_graph(info), set(chain)
    return all(_is_methylene_arm(graph, acid, main) for acid in _remaining_acids(info, chain))


def _all_neutral(info: dict) -> bool:
    """Round B accepts only wholly neutral acids; ion-state prefixes are Round D."""
    acids = info.get("carboxyls") or []
    return bool(acids) and not any(item.get("anion") for item in acids)


def carboxymethyl_diacid_parent(info: dict) -> dict | None:
    """Return the longest two-terminal-acid chain with only -CH2-COOH arms."""
    chain = _main_path(info)
    # Keep the Round B narrow subset disjoint from established short-chain
    # tricarboxylic parents, whose retained parent is selected by polycarboxylic.
    if not chain or len(chain) < 5 or len(_acids(info)) < 3 or not _all_neutral(info):
        return None
    if not (_is_open_saturated(info) and _only_acid_fgs(info) and _all_carboxymethyl_arms(info, chain)):
        return None
    return _parent(info, chain)


def _parent(info: dict, chain: list[int]) -> dict:
    acids = sorted(_remaining_acids(info, chain))
    facts = carboxymethyl_arms(info["mol"], chain, acids)
    return {"chain": chain, "n_carbons": len(chain), "kind": "diacid",
            "cooh_c_idxs": [chain[0], chain[-1]], "anion": False,
            "carboxymethyl_arms": facts}


def is_carboxymethyl_diacid(info: dict) -> bool:
    return carboxymethyl_diacid_parent(info) is not None
