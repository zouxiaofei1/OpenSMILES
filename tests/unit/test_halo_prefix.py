# IUPAC: P-61.3.1
# Layer: L3,L4,L5
"""Generic halogen substituent prefixes (fluoro/chloro/bromo/iodo) on alkane parents.

Affiliated: P-14.3.4 locant omission, P-14.5 alphabetical order, multiplicative di/tri/tetra.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: mono-halo, multi-halo, with alkyl, long chain
    ("CCCl", "chloroethane", "氯乙烷"),
    ("CCBr", "bromoethane", "溴乙烷"),
    ("CF", "fluoromethane", None),
    ("ClCCCl", "1,2-dichloroethane", "1,2-二氯乙烷"),
    ("CC(Cl)C", "2-chloropropane", "2-氯丙烷"),
    ("CC(Br)CC", "2-bromobutane", "2-溴丁烷"),
    ("CCCCF", "1-fluorobutane", "1-氟丁烷"),
    ("CCCCCI", "1-iodopentane", "1-碘戊烷"),
    ("CC(C)CCl", "1-chloro-2-methylpropane", "1-氯-2-甲基丙烷"),
    ("BrC(C)CC(C)C", "2-bromo-4-methylpentane", "2-溴-4-甲基戊烷"),
    ("BrCBr", "dibromomethane", None),
    ("ClCCCCCCCCCC", "1-chlorodecane", "1-氯癸烷"),
    # negative: alcohols / unsubstituted alkane must stay correct
    ("CCO", "ethanol", "乙醇"),
    ("CC(C)CO", "2-methylpropan-1-ol", "2-甲基丙-1-醇"),
    ("CC(C)C", "2-methylpropane", "2-甲基丙烷"),
    ("CCCO", "propan-1-ol", "丙-1-醇"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_halo_prefix(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
