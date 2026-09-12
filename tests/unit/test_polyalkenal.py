# IUPAC: P-66.6.1 / P-31.1 / P-93.4
# Layer: L5
from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer
import pytest

CASES = [
    ("C=CC=CC=O", "penta-2,4-dienal", "戊-2,4-二烯醛"),
    ("CC=CC=CC=O", "hexa-2,4-dienal", "己-2,4-二烯醛"),
    ("CC(C)=CCC/C(C)=C/C=O", "(2E)-3,7-dimethylocta-2,6-dienal", "(2E)-3,7-二甲基辛-2,6-二烯醛"),
    # negatives
    ("CC=CC=O", "but-2-enal", "丁-2-烯醛"),
    ("CCCC=O", "butanal", "丁醛"),
    ("C=CC(=O)O", "prop-2-enoic acid", "丙-2-烯酸"),
]

@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_polyalkenal(smiles, en, zh):
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)
