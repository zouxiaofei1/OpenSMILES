"""Terminal parent ownership: immutable owned_atoms + ordered candidates."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.tools.block_cut import parent_atom_set
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import iter_parent_candidates, select_parent


def _info(smiles: str):
    mol = preprocess(smiles)
    return mol, analyze(mol)


def test_benzamide_owns_core_not_n_phenyl():
    """Benzamide owns aryl + amide C/N/O; N-phenyl carbons stay outside."""
    mol, info = _info("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(info)
    assert parent["kind"] == "benzamide"
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)

    am = info["amides"][0]
    assert am["c_idx"] in owned
    assert am["n_idx"] in owned
    # carbonyl O is owned
    for n in mol.GetAtomWithIdx(am["c_idx"]).GetNeighbors():
        if n.GetAtomicNum() == 8:
            assert n.GetIdx() in owned
    # aryl core (chain) owned
    for idx in parent["chain"]:
        assert idx in owned
    # N-phenyl carbons not owned
    for c in am.get("n_c_idxs") or []:
        assert c not in owned
        for n in mol.GetAtomWithIdx(c).GetNeighbors():
            if n.GetAtomicNum() == 6 and n.GetIsAromatic() and n.GetIdx() != c:
                assert n.GetIdx() not in owned


def test_acid_owns_carboxyl_c_and_both_oxygens():
    """Acid owns carboxyl carbon + both oxygens (carbonyl O and OH O)."""
    mol, info = _info("CC(=O)O")
    parent = select_parent(info)
    assert parent["kind"] == "acid"
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)

    c_idx = parent["cooh_c_idx"]
    assert c_idx in owned
    oxygens = [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 8
    ]
    assert len(oxygens) == 2
    for o in oxygens:
        assert o in owned


def test_owned_atoms_is_frozenset_and_immutable_after_ops():
    """owned_atoms is frozenset and unchanged after extraction / side ops."""
    mol, info = _info("c1ccccc1C(=O)O")
    parent = select_parent(info)
    owned = parent["owned_atoms"]
    assert isinstance(owned, frozenset)
    snapshot = frozenset(owned)

    # side operations must not mutate ownership truth
    _ = parent_atom_set(parent, mol)
    _ = finalize_parent_ownership(parent, mol)
    _ = list(owned)
    copy_out = set(owned)
    copy_out.add(999)

    assert parent["owned_atoms"] == snapshot
    assert parent["owned_atoms"] is owned
    assert isinstance(parent["owned_atoms"], frozenset)


def test_parent_atom_set_adapter_prefers_owned_atoms():
    """Compatibility adapter returns owned_atoms when present."""
    mol, info = _info("CC(=O)O")
    parent = select_parent(info)
    assert parent_atom_set(parent, mol) == parent["owned_atoms"]


def test_iter_parent_candidates_finalizes_ordered():
    """iter_parent_candidates returns ranked finalized parents; first matches select."""
    mol, info = _info("c1ccccc1C(=O)Nc2ccccc2")
    cands = iter_parent_candidates(info)
    assert cands
    assert all(isinstance(c.get("owned_atoms"), frozenset) for c in cands)
    best = select_parent(info)
    assert best["kind"] == cands[0]["kind"]
    assert best["owned_atoms"] == cands[0]["owned_atoms"]


def test_finalize_copies_once_with_owned_atoms():
    """finalize_parent_ownership copies candidate once with owned_atoms frozenset."""
    mol, info = _info("CC(=O)O")
    raw = {
        "chain": [0, 1],
        "n_carbons": 2,
        "kind": "acid",
        "cooh_c_idx": 1,
    }
    fin = finalize_parent_ownership(raw, mol)
    assert fin is not raw
    assert "owned_atoms" not in raw
    assert isinstance(fin["owned_atoms"], frozenset)
    assert 1 in fin["owned_atoms"]
    assert len([a for a in fin["owned_atoms"] if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]) == 2


def _ester_o_count(mol, owned):
    """Count O atoms in owned set that are neighbours of owned C atoms (ester O)."""
    return sum(
        1 for i in owned
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 8
        and any(n.GetIdx() in owned and n.GetAtomicNum() == 6 for n in mol.GetAtomWithIdx(i).GetNeighbors())
    )
