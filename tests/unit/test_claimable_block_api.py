"""Ownership-only claimable blocks: topology claims, no naming mode."""
from __future__ import annotations

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer3.claimable_block import SideSlot, claim_block, iter_claims
from namepredict.layer2.parent_selector import select_parent


def _info(smiles: str):
    mol = preprocess(smiles)
    return mol, analyze(mol)


def _benzene_owned(smiles: str):
    mol, info = _info(smiles)
    for cand in select_parent(info, all_candidates=True):
        if cand.get("scaffold_id") == "benzene":
            return mol, cand["owned_atoms"]
    raise AssertionError("no benzene parent candidate")


def test_n_phenyl_benzamide_amide_n_claim():
    """N-phenyl benzamide: one AMIDE_N claim with six phenyl atoms."""
    mol, info = _info("c1ccccc1C(=O)Nc2ccccc2")
    parent = select_parent(info)
    assert parent["kind"] == "amide"
    assert parent.get("scaffold_id") == "benzene"
    claims = iter_claims(mol, parent["owned_atoms"])
    assert len(claims) == 1
    c = claims[0]
    assert c.slot == SideSlot.AMIDE_N
    assert len(c.atoms) == 6
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in c.atoms)
    assert c.attach_parent in parent["owned_atoms"]
    assert c.root not in parent["owned_atoms"]
    assert c.root in c.atoms


def test_methylbutylbenzene_ring_c_claim():
    """2-methylbutylbenzene with ring ownership: one RING_C claim, five side C."""
    mol, owned = _benzene_owned("c1ccc(cc1)CC(C)CC")
    claims = iter_claims(mol, owned)
    assert len(claims) == 1
    c = claims[0]
    assert c.slot == SideSlot.RING_C
    assert len(c.atoms) == 5
    assert all(mol.GetAtomWithIdx(i).GetAtomicNum() == 6 for i in c.atoms)
    assert all(not mol.GetAtomWithIdx(i).IsInRing() for i in c.atoms)
    assert mol.GetAtomWithIdx(c.attach_parent).IsInRing()


def test_multi_attachment_component_rejected():
    """Component touching two owned atoms (tetralin aliphatic bridge) is rejected."""
    mol = preprocess("c1cccc2c1CCCC2")
    owned = frozenset(
        i
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetIsAromatic()
    )
    claims = iter_claims(mol, owned)
    assert claims == []
    # explicit claim_block also returns None for either attachment
    roots = [
        n.GetIdx()
        for p in owned
        for n in mol.GetAtomWithIdx(p).GetNeighbors()
        if n.GetAtomicNum() != 1 and n.GetIdx() not in owned
    ]
    assert roots
    for root in roots:
        attach = next(
            p
            for p in owned
            for n in mol.GetAtomWithIdx(p).GetNeighbors()
            if n.GetIdx() == root
        )
        assert (
            claim_block(
                mol,
                owned_atoms=owned,
                attach_parent=attach,
                root=root,
                slot=SideSlot.RING_C,
            )
            is None
        )


def test_claim_block_valid_returns_atoms():
    """claim_block returns ClaimedBlock for a valid single-attachment side."""
    mol, owned = _benzene_owned("c1ccc(cc1)CC")
    claims = iter_claims(mol, owned)
    assert len(claims) == 1
    c = claims[0]
    again = claim_block(
        mol,
        owned_atoms=owned,
        attach_parent=c.attach_parent,
        root=c.root,
        slot=c.slot,
    )
    assert again == c


def test_ketone_carbonyl_o_not_claimed():
    """外部羰基氧（双键连所属碳）被跳过，不作侧链 claim（主 FG 已处理）。"""
    mol, info = _info("OC(=O)C(=O)C")
    parent = select_parent(info)
    assert parent["kind"] == "acid"
    claims = iter_claims(mol, parent["owned_atoms"])
    assert all(
        mol.GetAtomWithIdx(c.root).GetAtomicNum() != 8 for c in claims
    )


def test_sulfonyl_o_not_claimed():
    """砜双键氧同样被跳过，避免 cut 出 *O 污染成羟基。"""
    mol, info = _info("CS(=O)(=O)C")
    parent = select_parent(info)
    claims = iter_claims(mol, parent["owned_atoms"])
    assert all(
        mol.GetAtomWithIdx(c.root).GetAtomicNum() != 8 for c in claims
    )
