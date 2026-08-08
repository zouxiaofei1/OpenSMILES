"""Secondary amino substituents on non-amine parents (e.g. N-benzyl on ethanol)."""
from __future__ import annotations

from rdkit.Chem import Mol
import namepredict.layer3.side_facts as side_facts
from namepredict.layer3.aryl_names import aryl_arm_name


def _make_amino(attach: int, n_idx: int, en: str = "amino", zh: str = "氨基",
                atoms: list[int] | None = None, paren: bool = False) -> dict:
    return {
        "kind": "amino", "attach_idx": attach,
        "atoms": [n_idx] if atoms is None else list(atoms),
        "en": en, "zh": zh, "paren": paren,
    }


def _n_aryl_named(mol: Mol, start: int, n_idx: int, is_ch2: bool, owned):
    return None



def _sec_n_side_name(mol: Mol, start: int, n_idx: int, owned):
    return _n_aryl_named(mol, start, n_idx, True, owned) or _n_aryl_named(mol, start, n_idx, False, owned)


def _sec_on_off(a: dict, chain_set: set[int]) -> tuple[int, int] | None:
    cs = a.get("c_idxs") or []
    if a.get("degree") != 2 or len(cs) != 2:
        return None
    on = [c for c in cs if c in chain_set]
    off = [c for c in cs if c not in chain_set]
    return (on[0], off[0]) if len(on) == 1 and len(off) == 1 else None


def _sec_amino_off_chain(mol: Mol, a: dict, chain_set: set[int], owned) -> dict | None:
    pair = _sec_on_off(a, chain_set)
    if pair is None:
        return None
    named = _sec_n_side_name(mol, pair[1], a["n_idx"], owned)
    if named is None:
        return None
    en, zh, atoms = named
    return _make_amino(pair[0], a["n_idx"], en, zh, atoms, True)


def _one_amino(info: dict, a: dict, chain_set: set[int], owned) -> dict | None:
    if "c_idx" in a and a["c_idx"] in chain_set:
        return _make_amino(a["c_idx"], a["n_idx"])
    return _sec_amino_off_chain(info["mol"], a, chain_set, owned)


def _principal_amine_attachments(parent: dict) -> frozenset[int]:
    facts = parent.get("principal_expression_facts")
    return facts.attachment_atoms if facts and facts.group_class.value == "amine" else frozenset()


def _extract_aminos(info: dict, parent: dict, parent_nh2_kinds: set) -> list[dict]:
    principal = _principal_amine_attachments(parent)
    if parent.get("kind") in parent_nh2_kinds and not principal:
        return []
    chain_set = set(parent.get("chain") or [])
    owned = parent.get("owned_atoms")
    out: list[dict] = []
    for a in info.get("amines") or []:
        one = _one_amino(info, a, chain_set, owned)
        if one is not None and one["attach_idx"] not in principal:
            out.append(one)
    return out
