"""N-block extract: complex N-substituents via SubstituentNamer (amide/benzamide)."""
from __future__ import annotations

from namepredict.layer2.claimable_block import ClaimedBlock, SideSlot
from namepredict.layer3.substituent_namer import SubstituentName, SubstituentNamer


def _n_attach(parent: dict) -> int | None:
    if parent.get("kind") == "benzamide":
        return parent.get("ring_attach_idx")
    return parent.get("amide_c_idx")


def _amide_n(parent: dict, mol) -> int:
    c = parent.get("amide_c_idx")
    if c is None:
        return -1
    for n in mol.GetAtomWithIdx(c).GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return -1


def _n_claim(mol, parent: dict) -> ClaimedBlock | None:
    from namepredict.layer2.block_cut import cut_block, parent_atom_set

    root = parent.get("n_block_root")
    if root is None:
        return None
    atoms = cut_block(mol, root, parent_atom_set(parent, mol))
    if atoms is None:
        return None
    return ClaimedBlock(
        slot=SideSlot.AMIDE_N,
        attach_parent=_amide_n(parent, mol),
        root=root,
        atoms=frozenset(atoms),
    )


def _wrap_n(en: str, zh: str, paren: bool) -> tuple[str, str]:
    if paren and not (en.startswith("(") and en.endswith(")")):
        return f"N-({en})", f"N-({zh})"
    return f"N-{en}", f"N-{zh}"


def _n_block_dict(parent: dict, named: SubstituentName) -> dict:
    en, zh = _wrap_n(named.en, named.zh, named.requires_parentheses)
    return {
        "kind": "n_block",
        "n_carbons": 0,
        "attach_idx": _n_attach(parent),
        "atoms": sorted(named.claim.atoms),
        "en": en,
        "zh": zh,
        "paren": False,
        "backend": named.backend,
    }


def extract_n_blocks(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") not in {"amide", "benzamide"} or not parent.get("n_block"):
        return []
    mol = info["mol"]
    claim = _n_claim(mol, parent)
    if claim is None or claim.attach_parent < 0:
        return []
    named = SubstituentNamer().name(mol, claim, depth=0)
    return [] if named is None else [_n_block_dict(parent, named)]
