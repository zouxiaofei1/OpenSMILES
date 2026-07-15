# IUPAC: P-22.2.1 / P-25
# Layer: L2
"""Data-driven fused 5+6 1,3-dihetero engine (benzothiazole / benzoxazole)."""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.fused56 import (
    Fused56Di13Spec,
    _di13_parts,
    _try_di13_amine,
    _try_di13_fused56,
)
from namepredict.namer import SMILESNNamer

BTZ = Fused56Di13Spec(
    kind="benzothiazole",
    hetero_z=16,
    hetero_key="s_idx",
    amine_kind="benzothiazolamine",
)
BOX = Fused56Di13Spec(
    kind="benzoxazole",
    hetero_z=8,
    hetero_key="o_idx",
    amine_kind="benzoxazolamine",
)


def _info(smiles: str) -> dict:
    mol = preprocess(smiles)
    assert mol is not None
    return analyze(mol)


def test_btz_parts_unsubstituted() -> None:
    p = _di13_parts(_info("c1nc2ccccc2s1"), BTZ)
    assert p is not None
    five, six, h, n, ba, bb = p
    assert len(five) == 5 and len(six) == 6
    assert h is not None and n is not None


def test_box_parts_unsubstituted() -> None:
    p = _di13_parts(_info("c1ccc2ocnc2c1"), BOX)
    assert p is not None
    assert p[2] is not None and p[3] is not None


def test_btz_try_parent_kind() -> None:
    parent = _try_di13_fused56(_info("c1nc2ccccc2s1"), BTZ)
    assert parent is not None
    assert parent["kind"] == "benzothiazole"
    assert "s_idx" in parent and "n_idx" in parent
    assert parent.get("nh_idx") == parent["s_idx"]
    assert len(parent["chain"]) == 9


def test_box_try_parent_kind() -> None:
    parent = _try_di13_fused56(_info("c1ccc2ocnc2c1"), BOX)
    assert parent is not None
    assert parent["kind"] == "benzoxazole"
    assert "o_idx" in parent and "n_idx" in parent
    assert len(parent["chain"]) == 9


def test_btz_amine_via_engine() -> None:
    parent = _try_di13_amine(_info("Nc1nc2ccccc2s1"), BTZ)
    assert parent is not None
    assert parent["kind"] == "benzothiazolamine"
    assert parent.get("amine_c_idx") is not None
    assert len(parent["chain"]) == 9


def test_box_halo_via_engine() -> None:
    parent = _try_di13_fused56(_info("Brc1ccc2ncoc2c1"), BOX)
    assert parent is not None
    assert parent["kind"] == "benzoxazole"


def test_benzofuran_not_box() -> None:
    assert _try_di13_fused56(_info("c1coc2ccccc12"), BOX) is None
    assert _try_di13_amine(_info("c1coc2ccccc12"), BOX) is None


def test_pyridine_not_btz() -> None:
    assert _try_di13_fused56(_info("c1ccncc1"), BTZ) is None


def test_cross_spec_rejects() -> None:
    """S-core must not match O-spec and vice versa."""
    assert _try_di13_fused56(_info("c1nc2ccccc2s1"), BOX) is None
    assert _try_di13_fused56(_info("c1ccc2ocnc2c1"), BTZ) is None


# e2e bilingual via public thin wrappers
_E2E = [
    ("c1nc2ccccc2s1", "1,3-benzothiazole", "1,3-苯并噻唑"),
    ("c1ccc2ocnc2c1", "1,3-benzoxazole", "1,3-苯并噁唑"),
    ("Nc1nc2ccccc2s1", "1,3-benzothiazol-2-amine", "2-氨基苯并噻唑"),
    ("Brc1ccc2ncoc2c1", "6-bromo-1,3-benzoxazole", "6-溴-1,3-苯并噁唑"),
    ("c1coc2ccccc12", "benzofuran", "苯并呋喃"),
    ("c1ccncc1", "pyridine", "吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", _E2E)
def test_e2e_fused56_di13(smiles: str, en: str, zh: str) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
