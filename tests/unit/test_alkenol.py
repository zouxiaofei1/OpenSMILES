# IUPAC: P-63.1.1 / P-31.1 / P-93.4
# Layer: L2,L4,L5
"""Open-chain unsaturated monoalcohols (alkenols / polyalkenols).

Alcohol is the principal characteristic group; non-aromatic C=C is inserted as
-ene_locant-en-OH_locant-ol (mono) or -a,b-dien-OH-ol (poly). Parent chain covers
the OH carbon and all double-bond carbons. Numbering: lowest OH locant first,
then lowest ene set. BondStereo → (E)-/(Z)- or (2E,6Z)- prefixes (P-93.4).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # C2 单烯醇：烯 1-2 + OH 1 无歧义，位次省略融合 ethenol (P-14.3.4)
    ("C=CO", "ethenol", "乙烯醇"),
    # mono-alkenols
    ("C=CCCO", "but-3-en-1-ol", "丁-3-烯-1-醇"),
    ("C=CCO", "prop-2-en-1-ol", "丙-2-烯-1-醇"),
    ("CC(O)C=C", "but-3-en-2-ol", "丁-3-烯-2-醇"),
    ("CC=CCO", "but-2-en-1-ol", "丁-2-烯-1-醇"),
    ("C=CCCCCO", "hex-5-en-1-ol", "己-5-烯-1-醇"),
    ("CC(C)=CCO", "3-methylbut-2-en-1-ol", "3-甲基丁-2-烯-1-醇"),
    # E/Z mono
    (r"C(\C=C\CCCCCCCCC)O", "(2E)-dodec-2-en-1-ol", "(2E)-十二-2-烯-1-醇"),
    (r"CCCC/C=C\CCCCCCCCCCCCO", "(13Z)-octadec-13-en-1-ol", None),
    (r"C/C=C/CCO", "(3E)-pent-3-en-1-ol", "(3E)-戊-3-烯-1-醇"),
    # polyalkenols
    ("C=CC=CCO", "penta-2,4-dien-1-ol", "戊-2,4-二烯-1-醇"),
    (r"CC/C=C\CC/C=C/CO", "(2E,6Z)-nona-2,6-dien-1-ol", "(2E,6Z)-壬-2,6-二烯-1-醇"),
    (
        r"CC/C=C\C/C=C\C/C=C\CCCCCCCCO",
        "(9Z,12Z,15Z)-octadeca-9,12,15-trien-1-ol",
        None,
    ),
    # polyalkenols beyond triene (multiplicity-generic mult-seg)
    ("C=CC=CC=CC=CCO", "nona-2,4,6,8-tetraen-1-ol", "壬-2,4,6,8-四烯-1-醇"),
    ("C=CC=CC=CC=CC=CCO", "undeca-2,4,6,8,10-pentaen-1-ol", "十一-2,4,6,8,10-五烯-1-醇"),
    # negatives
    ("CCCCO", "butan-1-ol", "丁-1-醇"),
    ("C1CCC(O)CC1", "cyclohexanol", "环己醇"),
    ("CCCCCCCCCCCC(=O)[O-]", "dodecanoate", "十二酸根"),
    ("C=CC=CC=CC", "hepta-1,3,5-triene", "庚-1,3,5-三烯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_alkenol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
