"""Open-chain polycarboxylic parent selection (IUPAC P-65.1.1)."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit.Chem.rdchem import Mol
_ALKANE_NAMES = {
    1: ("methane", "甲烷"), 2: ("ethane", "乙烷"), 3: ("propane", "丙烷"),
    4: ("butane", "丁烷"), 5: ("pentane", "戊烷"), 6: ("hexane", "己烷"),
    7: ("heptane", "庚烷"), 8: ("octane", "辛烷"), 9: ("nonane", "壬烷"),
    10: ("decane", "癸烷"),
}

def _alkane_names(n: int) -> tuple[str, str] | None:
    return _ALKANE_NAMES.get(n)



@dataclass(frozen=True)
class PolycarboxylicEligibility:
    """L2 eligibility fact for structures with three or more carboxyl groups."""

    supported: bool
    reason: str | None = None


def _carboxyl_carbons(info: dict) -> set[int]:
    return {int(e["c_idx"]) for e in info.get("carboxyls") or []}


def _core_carbons(mol: Mol, acids: set[int]) -> set[int]:
    return {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIdx() not in acids}


def _carbon_neighbors(mol: Mol, idx: int, allowed: set[int]) -> list[int]:
    return [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors() if n.GetIdx() in allowed]


def _path_from(mol: Mol, start: int, allowed: set[int]) -> list[int]:
    stack, paths = [(start, [start])], []
    while stack:
        cur, path = stack.pop(); paths.append(path)
        stack.extend((n, path + [n]) for n in _carbon_neighbors(mol, cur, allowed) if n not in path)
    return max(paths, key=len)


def _longest_core_path(mol: Mol, core: set[int]) -> list[int]:
    paths = [_path_from(mol, atom, core) for atom in core]
    return max(paths, key=len, default=[])


def _acid_attachment(mol: Mol, acid: int, core: set[int]) -> int | None:
    hits = _carbon_neighbors(mol, acid, core)
    return hits[0] if len(hits) == 1 else None


def _acid_attachments(mol: Mol, acids: set[int], core: set[int]) -> list[int] | None:
    attached = [_acid_attachment(mol, acid, core) for acid in acids]
    return attached if all(atom is not None for atom in attached) else None


def _is_acyclic_carbon_set(mol: Mol, core: set[int]) -> bool:
    return not any(mol.GetAtomWithIdx(atom).IsInRing() for atom in core)


def _has_only_polyacid_fgs(info: dict) -> bool:
    forbidden = ("has_ester", "has_amide", "has_ketone", "has_aldehyde", "has_nitrile")
    return not any(info.get(key) for key in forbidden)


def _is_supported_count(acids: set[int]) -> bool:
    return 3 <= len(acids) <= 10


def _all_attachments_on_chain(chain: list[int], attachments: list[int]) -> bool:
    return all(atom in chain for atom in attachments)


def _core_path(info: dict) -> tuple[list[int], list[int]] | None:
    mol, acids = info["mol"], _carboxyl_carbons(info)
    core = _core_carbons(mol, acids); chain = _longest_core_path(mol, core)
    attached = _acid_attachments(mol, acids, core)
    return (chain, attached) if attached and set(chain) == core and _all_attachments_on_chain(chain, attached) else None


def _all_deprotonated(info: dict) -> bool:
    acids = info.get("carboxyls") or []
    return bool(acids) and all(acid.get("anion") for acid in acids)


def _core_bonds(info: dict, key: str, core: set[int]) -> list[tuple[int, int]]:
    return [(b["c1"], b["c2"]) for b in info.get(key) or [] if b["c1"] in core and b["c2"] in core]


def _bond_meta(info: dict, core: set[int]) -> dict:
    doubles = _core_bonds(info, "double_bonds", core)
    triples = _core_bonds(info, "triple_bonds", core)
    return {"double_bonds": doubles, "triple_bonds": triples,
            "double_bond": doubles[0] if len(doubles) == 1 else None,
            "triple_bond": triples[0] if len(triples) == 1 else None}


def _has_external_unsaturation(info: dict, core: set[int]) -> bool:
    bonds = _core_bonds(info, "double_bonds", core) + _core_bonds(info, "triple_bonds", core)
    return len(bonds) != len(info.get("double_bonds") or []) + len(info.get("triple_bonds") or [])


def _polyacid_ok(info: dict) -> bool:
    acids = _carboxyl_carbons(info); mol = info["mol"]
    core = _core_carbons(mol, acids)
    return _is_supported_count(acids) and _has_only_polyacid_fgs(info) and _is_acyclic_carbon_set(mol, core) and not _has_external_unsaturation(info, core)


def _polycarboxylic_block_reason(info: dict) -> str | None:
    if has_partial_deprotonation(info):
        return "partial_deprotonation"
    if not _polyacid_ok(info) or _core_path(info) is None:
        return "unsupported_skeleton"
    return None


def polycarboxylic_eligibility(info: dict) -> PolycarboxylicEligibility | None:
    """Classify >=3-carboxyl structures before candidate production."""
    if len(_carboxyl_carbons(info)) < 3:
        return None
    acids = _carboxyl_carbons(info)
    if len(acids) == 3 and not _is_acyclic_carbon_set(info["mol"], _core_carbons(info["mol"], acids)):
        return None
    reason = _polycarboxylic_block_reason(info)
    return PolycarboxylicEligibility(reason is None, reason)


def try_polycarboxylic_parent(info: dict) -> dict | None:
    eligibility = polycarboxylic_eligibility(info)
    if eligibility is None or not eligibility.supported:
        return None
    got = _core_path(info)
    return _parent(info, got, len(_carboxyl_carbons(info)), _all_deprotonated(info)) if got else None


def _parent(info: dict, got: tuple[list[int], list[int]], count: int, anion: bool) -> dict:
    chain, attach = got; core = set(chain)
    names = _alkane_names(len(chain))
    stems = {} if names is None else {"stem_en": names[0], "stem_zh": names[1]}
    return {"chain": chain, "n_carbons": len(chain), "kind": "polycarboxylic", "cooh_c_idxs": attach,
            "acid_count": count, "anion": anion, "mol": info["mol"], **stems, **_bond_meta(info, core)}


def has_partial_deprotonation(info: dict) -> bool:
    acids = info.get("carboxyls") or []
    states = {bool(acid.get("anion")) for acid in acids}
    return len(acids) >= 3 and len(states) == 2
