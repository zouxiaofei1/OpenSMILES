# IUPAC: P-63.2.2
# Layer: L2,L3,L5
"""Aromatic alkoxy prefixes (methoxy/ethoxy) and retained anisole (P-63.2.2).

Parent = benzene or phenol; ring alkoxy O–R (R = Me/Et linear) is a prefix.
Unsubstituted methoxybenzene → retained anisole / 甲氧基苯（对齐 benchmark 金标）.
Do not absorb open-chain ethers or bare phenol/benzene.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: retained anisole
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    # positive: methoxyphenol
    ("COc1ccc(O)cc1", "4-methoxyphenol", "4-甲氧基苯酚"),
    # positive: methoxy + halo / methyl (lowest set of locants)
    ("COc1ccccc1I", "1-iodo-2-methoxybenzene", None),
    ("COc1ccccc1C", "1-methoxy-2-methylbenzene", None),
    # positive: ethoxybenzene (systemic; no retained name)
    ("CCOc1ccccc1", "ethoxybenzene", None),
    # negative: bare / phenol / open ether / alcohol must not break
    ("c1ccccc1", "benzene", "苯"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CCOCC", "diethyl ether", "二乙基醚"),
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCO", "ethanol", "乙醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_anisole(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_anisole_not_hexane_or_ether() -> None:
    """Aryl methoxy must not fall to hexane parent or dialkyl ether."""
    r = SMILESNNamer().name("COc1ccccc1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "anisole"
    assert "hexane" not in en
    assert "ether" not in en
