# IUPAC: P-25.4
# Layer: L2,L4
"""Asymmetric fused arene locant numbering (P-25.4): quinoline/quinazoline/indole.

These fused heteroarenes have fixed IUPAC numbering (heteroatom = 1, then along
the ring with bridgehead letters 4a/7a). P-14.4 generic ring enumeration
mis-numbers them (10-chloroquinolin-5-ol instead of 2-chloroquinolin-6-ol).
The retained scaffold numbering facts (standard_path + template match) drive L4.
Symmetric templates (naphthalene / quinoxaline / benzimidazole) are excluded:
their GetSubstructMatch direction is not unique, so P-14.4's substituent
lowest-locant rule stays in control (fixed numbering would be unstable).
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # positive: quinoline locants (N=1; 2-chloro, 6-ol)
    ("ClC1=NC2=CC=C(C=C2C=C1)O", "2-chloroquinolin-6-ol", "2-氯喹啉-6-醇"),
    ("ClC1=NC2=CC=C(C=C2C(=C1)C)O", "2-chloro-4-methylquinolin-6-ol", "2-氯-4-甲基喹啉-6-醇"),
    # positive: quinoline CF3 (zh 括号文体与 gold 不同，只断言 EN locant)
    ("BrC1=CC=C2C=CC(=NC2=C1)C(F)(F)F", "7-bromo-2-(trifluoromethyl)quinoline", None),
    # positive: quinazoline locants (N1,N3; 2-chloro, 4-methyl, 7-methoxy)
    ("ClC1=NC2=CC(=CC=C2C(=N1)C)OC", "2-chloro-7-methoxy-4-methylquinazoline", "2-氯-7-甲氧基-4-甲基喹唑啉"),
    ("ClC1=NC(=NC2=CC=C(C=C12)F)C1=CC=C(C=C1)OC", "4-chloro-6-fluoro-2-(4-methoxyphenyl)quinazoline", "4-氯-6-氟-2-(4-甲氧基苯基)喹唑啉"),
    # negative near-misses: symmetric naphthalene/quinoxaline stay P-14.4; already-correct cases must not regress
    ("c1ccc2ccccc2c1", "naphthalene", "萘"),
    ("Cc1cccc2ccccc12", "1-methylnaphthalene", "1-甲基萘"),
    ("BrCC(=O)C1=NC2=CC=CC=C2N=C1C", "2-bromo-1-(3-methylquinoxalin-2-yl)ethanone", None),
    ("ClC1=NC2=CC=CC=C2C(=C1)C(F)(F)F", "2-chloro-4-(trifluoromethyl)quinoline", None),
    ("c1ccc2[nH]ccc2c1", "1H-indole", "1H-吲哚"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_fused_locant_numbering(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)
