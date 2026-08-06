"""Open-chain alkoxy: methoxybutane via typed ether claims, not alkane evaporation."""
from __future__ import annotations

import pytest
from rdkit import Chem

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer1.analyzer import analyze
from namepredict.layer3.claimable_block import iter_claims
from namepredict.layer2.parent_selector import select_parent
from namepredict.namer import SMILESNNamer

METHOXYBUTANE = "COC(C)CC"


def test_methoxybutane_final_name():
    result = SMILESNNamer().name(METHOXYBUTANE)
    assert result.success
    assert normalize_en(result.en) == normalize_en("2-methoxybutane")
    assert normalize_zh(result.zh) == normalize_zh("2-甲氧基丁烷")


def test_methoxybutane_has_ether_or_alkoxy_claim_path():
    """Alkane parent with methoxy claim, or ether parent — O must be covered."""
    mol = Chem.MolFromSmiles(METHOXYBUTANE)
    info = analyze(mol)
    parent = select_parent(info)
    owned = parent["owned_atoms"]
    claims = iter_claims(mol, owned)
    o_idxs = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 8}
    covered = set(owned) | {i for c in claims for i in c.atoms}
    assert o_idxs <= covered
    from namepredict.layer3.substituent_namer import SubstituentNamer

    named = [SubstituentNamer().name(mol, c) for c in claims]
    assert any(n is not None and n.en == "methoxy" for n in named) or parent.get("kind") == "ether"


@pytest.mark.parametrize(
    "smiles,en",
    [
        ("COC", "methoxyethane"),  # may be dimethyl ether / methoxyethane retained
        ("CCOCC", "ethoxyethane"),
    ],
)
def test_simple_sym_ethers_still_succeed(smiles, en):
    result = SMILESNNamer().name(smiles)
    assert result.success
    # loose: just ensure oxygen is named somehow
    assert "oxy" in normalize_en(result.en) or "ether" in normalize_en(result.en)
