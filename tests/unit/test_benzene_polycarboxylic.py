# IUPAC: P-65.1.1/P-65.1.2
# Layer: L2,L4,L5
"""Neutral, unfused benzene di- and tricarboxylic acid parents."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.arene_carbonyl import benzene_polycarboxylic_eligibility
from namepredict.layer2.candidates import _collect_candidates
from namepredict.layer2.kind_registry import parent_names
from namepredict.namer import SMILESNNamer

# (SMILES, English PIN, strict Chinese name); last case preserves retained monoacid.
CASES = [
    ("O=C(O)c1ccccc1C(=O)O", "benzene-1,2-dicarboxylic acid", "苯-1,2-二羧酸"),
    ("O=C(O)c1cccc(C(=O)O)c1", "benzene-1,3-dicarboxylic acid", "苯-1,3-二羧酸"),
    ("O=C(O)c1ccc(C(=O)O)cc1", "benzene-1,4-dicarboxylic acid", "苯-1,4-二羧酸"),
    ("O=C(O)c1c(C(=O)O)c(C(=O)O)ccc1", "benzene-1,2,3-tricarboxylic acid", "苯-1,2,3-三羧酸"),
    ("O=C(O)c1ccccc1", "benzoic acid", "苯甲酸"),
]


def test_benzene_polyacid_has_authoritative_parent_stem() -> None:
    assert parent_names("benzene_polycarboxylic") == ("benzene", "苯")
def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def _candidate_kinds(smiles: str) -> set[str]:
    return {c["kind"] for c in _collect_candidates(_info(smiles))}


@pytest.mark.parametrize("smiles", [
    "O=C(O)C1CCCCC1C(=O)O",
    "O=C(O)c1ccncc1C(=O)O",
    "O=C(O)c1occc1C(=O)O",
])
def test_nonbenzene_polyacids_bypass_benzene_gate(smiles: str) -> None:
    info = _info(smiles)
    assert benzene_polycarboxylic_eligibility(info) is None
    assert "unsupported_polycarboxylic" not in _candidate_kinds(smiles)


def test_unsupported_benzene_polyacid_is_still_blocked() -> None:
    assert _candidate_kinds("O=C(O)c1ccccc1C(=O)OC") == {"unsupported_polycarboxylic"}


@pytest.mark.parametrize(("smiles", "en"), [
    ("O=C(O)CC(Cc1ccccc1)C(=O)O", "2-benzylbutanedioic acid"),
    ("O=C([O-])CC(Cc1ccccc1)C(=O)[O-]", "2-benzylbutanedioate"),
])
def test_benzylbutanedioic_acids_bypass_benzene_gate(smiles: str, en: str) -> None:
    info = _info(smiles)
    assert benzene_polycarboxylic_eligibility(info) is None
    assert "unsupported_polycarboxylic" not in _candidate_kinds(smiles)
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)


@pytest.mark.parametrize("smiles", [
    "O=C([O-])c1ccccc1C(=O)O",  # partial salt
    "O=C(O)c1c(C(=O)O)c(C(=O)O)c(C(=O)O)cc1",  # tetraacid
    "COC(=O)c1ccccc1C(=O)O",  # ester
])
def test_benzene_polyacid_exclusions_are_rejected(smiles: str) -> None:
    assert not SMILESNNamer().name(smiles).success
