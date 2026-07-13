# IUPAC: P-62.2.1 / P-63.1.4 / P-22.2.1
# Layer: L2,L4,L5
"""Retained pyridine with amino/hydroxy: pyridin-n-amine / pyridin-n-ol.

N fixed as locant 1; principal NH2 or OH on a ring carbon gets locant 2/3/4.
Ring may carry at most one simple extra sub: mono-halo (F/Cl/Br/I) or mono-methyl.
Orientation under N=1 chooses direction by lowest set of FG + ring sub locants.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

CASES = [
    # positive: unsubstituted pyridinamines
    ("Nc1ccccn1", "pyridin-2-amine", "吡啶-2-胺"),
    ("Nc1cccnc1", "pyridin-3-amine", "吡啶-3-胺"),
    ("Nc1ccncc1", "pyridin-4-amine", "吡啶-4-胺"),
    # positive: unsubstituted pyridinols
    ("Oc1ccccn1", "pyridin-2-ol", "吡啶-2-醇"),
    ("Oc1cccnc1", "pyridin-3-ol", "吡啶-3-醇"),
    ("Oc1ccncc1", "pyridin-4-ol", "吡啶-4-醇"),
    # positive: mono-methyl pyridinamines
    ("Cc1cccc(N)n1", "6-methylpyridin-2-amine", "6-甲基吡啶-2-胺"),
    ("Nc1c(C)cccn1", "3-methylpyridin-2-amine", "3-甲基吡啶-2-胺"),
    ("Nc1cc(C)ccn1", "4-methylpyridin-2-amine", "4-甲基吡啶-2-胺"),
    ("Nc1ccc(C)cn1", "5-methylpyridin-2-amine", "5-甲基吡啶-2-胺"),
    ("Nc1cccnc1C", "2-methylpyridin-3-amine", "2-甲基吡啶-3-胺"),
    ("Nc1cnccc1C", "4-methylpyridin-3-amine", "4-甲基吡啶-3-胺"),
    ("Nc1cncc(C)c1", "5-methylpyridin-3-amine", "5-甲基吡啶-3-胺"),
    ("Nc1cc(C)ncc1", "2-methylpyridin-4-amine", "2-甲基吡啶-4-胺"),
    # lowest-set orientation: 2-methyl + 5-amine preferred over 6-methyl + 3-amine
    ("Nc1ccc(C)nc1", "2-methylpyridin-5-amine", "2-甲基吡啶-5-胺"),
    # positive: mono-halo pyridinamines
    ("Clc1cccc(N)n1", "6-chloropyridin-2-amine", "6-氯吡啶-2-胺"),
    ("Clc1ccncc1N", "4-chloropyridin-3-amine", "4-氯吡啶-3-胺"),
    ("Nc1cc(Cl)ccn1", "4-chloropyridin-2-amine", "4-氯吡啶-2-胺"),
    ("Nc1cccnc1Cl", "2-chloropyridin-3-amine", "2-氯吡啶-3-胺"),
    ("Brc1ccncc1N", "4-bromopyridin-3-amine", "4-溴吡啶-3-胺"),
    ("Fc1cccc(N)n1", "6-fluoropyridin-2-amine", "6-氟吡啶-2-胺"),
    # positive: mono-methyl / mono-halo pyridinols
    ("Oc1cccc(C)n1", "6-methylpyridin-2-ol", "6-甲基吡啶-2-醇"),
    ("Cc1cccnc1O", "3-methylpyridin-2-ol", "3-甲基吡啶-2-醇"),
    ("Oc1cc(Cl)ccn1", "4-chloropyridin-2-ol", "4-氯吡啶-2-醇"),
    ("Clc1ncccc1O", "2-chloropyridin-3-ol", "2-氯吡啶-3-醇"),
    ("Ic1cccnc1O", "3-iodopyridin-2-ol", "3-碘吡啶-2-醇"),
    # negative: bare pyridine, aniline, phenol, chain amine, pyridinecarboxylic
    ("c1ccncc1", "pyridine", "吡啶"),
    ("Nc1ccccc1", "aniline", "苯胺"),
    ("Oc1ccccc1", "phenol", "苯酚"),
    ("CCN", "ethanamine", "乙胺"),
    ("OC(=O)c1ccccn1", "pyridine-2-carboxylic acid", "吡啶-2-甲酸"),
    # negative: mono-methyl / mono-halo retained pyridine (no amine/ol FG)
    ("Cc1ccccn1", "2-methylpyridine", "2-甲基吡啶"),
    ("Clc1ccccn1", "2-chloropyridine", "2-氯吡啶"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_pyridinamine_ol(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_amino_locant_is_2_not_6() -> None:
    """N fixed as 1; 2-aminopyridine must be pyridin-2-amine not 6-."""
    r = SMILESNNamer().name("Nc1ccccn1")
    assert r.success
    en = normalize_en(r.en)
    assert en == "pyridin-2-amine"
    assert "6-amine" not in en
    assert "pyridin-6" not in en


def test_multi_sub_not_pyridinamine() -> None:
    """Dimethyl / dichloro rings (>1 extra sub) must not stay as pyridinamine parent."""
    for s in ("Cc1cc(C)nc(N)c1", "Clc1cc(Cl)nc(N)c1"):
        r = SMILESNNamer().name(s)
        assert r.success
        assert "pyridin" not in normalize_en(r.en)


def test_mono_ethyl_pyridinamine() -> None:
    """Recursive alkyl sides (P-29.3.1): one ethyl side is now a valid substituent."""
    r = SMILESNNamer().name("CCc1cccc(N)n1")
    assert r.success
    assert normalize_en(r.en) == normalize_en("6-ethylpyridin-2-amine")
    assert normalize_zh(r.zh) == normalize_zh("6-乙基吡啶-2-胺")
