# IUPAC: P-63.2.2
# Layer: L1,L2,L5
"""Open-chain simple dialkyl ethers (symmetric retained + alkoxyalkane).

Exactly one ether O (two C neighbors, no H, not ester/anhydride). Both arms
unsubstituted linear alkyl C1–C4. Symmetric: diethyl ether / 二乙基醚.
Asymmetric: methoxyethane / 甲氧基乙烷; 1-methoxypropane / 1-甲氧基丙烷.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: symmetric retained
    ("CCOCC", "diethyl ether", "二乙基醚"),
    ("COC", "dimethyl ether", "二甲基醚"),
    ("CCCOCCC", "dipropyl ether", "二丙基醚"),
    ("CCCCOCCCC", "dibutyl ether", "二丁基醚"),
    # positive: asymmetric alkoxyalkane
    ("CCOC", "methoxyethane", "甲氧基乙烷"),
    ("CCCOC", "1-methoxypropane", "1-甲氧基丙烷"),
    ("CCOCCC", "1-ethoxypropane", "1-乙氧基丙烷"),
    # negative: alcohol / ester / alkenoate must not break
    ("CCO", "ethanol", "乙醇"),
    ("CC(=O)OC", "methyl acetate", "乙酸甲酯"),
    ("CCOC(=O)C", "ethyl acetate", "乙酸乙酯"),
    ("C=CC(=O)OC", "methyl prop-2-enoate", "丙-2-烯酸甲酯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_dialkyl_ether(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_diethyl_ether_not_ethane() -> None:
    """CCOCC must not collapse to ethane parent."""
    r = SMILESNNamer().name("CCOCC")
    assert r.success
    assert normalize_en(r.en) == "diethyl ether"
    assert normalize_en(r.en) != "ethane"
