# IUPAC: P-65.1.1 / P-72.2.2.1
# Layer: L2,L4,L5
"""Open-chain tricarboxylic acids and their fully deprotonated anions."""
from __future__ import annotations

from namepredict.tools.re import normalize_en
from namepredict.namer import SMILESNNamer


def test_tricarboxylate_does_not_fall_back_to_monoacid() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)[O-])CC(=O)[O-]")
    assert result.success
    assert "acetate" not in normalize_en(result.en)
    assert "hexanoate" not in normalize_en(result.en)


def test_partial_deprotonation_does_not_claim_polycarboxylate() -> None:
    result = SMILESNNamer().name("O=C([O-])CC(C(=O)O)CC(=O)O")
    assert "tricarboxylate" not in normalize_en(result.en)
