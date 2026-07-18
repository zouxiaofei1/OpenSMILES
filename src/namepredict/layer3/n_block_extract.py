"""N-block extract: complex N-substituents via name_as_substituent (amide/benzamide)."""
from __future__ import annotations


def extract_n_blocks(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") not in {"amide", "benzamide"} or not parent.get("n_block"):
        return []
    mol, root = info["mol"], parent["n_block_root"]
    from namepredict.layer2.block_cut import cut_block, parent_atom_set
    from namepredict.layer3.as_substituent import name_as_substituent

    patoms = parent_atom_set(parent, mol)
    atoms = cut_block(mol, root, patoms)
    if atoms is None:
        return []
    named = name_as_substituent(mol, root, atoms, depth=0)
    if named is None:
        return []
    en, zh, paren = named
    attach = (
        parent.get("ring_attach_idx")
        if parent["kind"] == "benzamide"
        else parent.get("amide_c_idx")
    )
    # Parens live in en/zh (like n_benzyl); do not set paren — L5 would double-wrap.
    return [{
        "kind": "n_block",
        "n_carbons": 0,
        "attach_idx": attach,
        "atoms": list(atoms),
        "en": f"N-({en})",
        "zh": f"N-({zh})",
        "paren": False,
    }]
