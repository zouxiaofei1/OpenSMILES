# IUPAC: P-29.3 / P-14.3.1 / P-22.2.1
# Layer: L2,L3
"""Recursive depth-1 aryl: multi-halo Ph/OPh, benzyl, benzyloxy; pyridine parent.

- Ph/OPh: up to 3 terminal ring-halos; lowest locants from attach=1 (2- not 6-).
- benzyl = parent–CH2–Ph; benzyloxy = parent–O–CH2–Ph (not ethoxy).
- pyridine (unfused C5N) may carry phenyl/phenoxy/benzyl/benzyloxy.
Negatives: anisole / ethoxybenzene / phenoxybenzene must stay correct.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# ("smiles", "expected_en", "expected_zh_or_None")
CASES = [
    # multi-halo phenyl as substituent on FG parent (attach=1 → lowest set)
    (
        "O=Cc1ccc(-c2ccc(Cl)c(Cl)c2)cc1",
        "4-(3,4-dichlorophenyl)benzaldehyde",
        "4-(3,4-二氯苯基)苯甲醛",
    ),
    # multi-halo phenoxy on phenol
    (
        "Oc1ccc(Oc2cc(Cl)cc(Cl)c2)cc1",
        "4-(3,5-dichlorophenoxy)phenol",
        "4-(3,5-二氯苯氧基)苯酚",
    ),
    # ortho-halo phenyl lowest locant is 2- (not 6-) when Ph is the arm
    (
        "O=Cc1ccc(-c2ccccc2Cl)cc1",
        "4-(2-chlorophenyl)benzaldehyde",
        "4-(2-氯苯基)苯甲醛",
    ),
    # benzylbenzene
    ("c1ccc(Cc2ccccc2)cc1", "benzylbenzene", "苄基苯"),
    # benzyloxybenzene (must not be phenoxyheptane / ethoxy)
    ("c1ccc(OCc2ccccc2)cc1", "benzyloxybenzene", "苄氧基苯"),
    # 4-benzyloxybenzaldehyde
    (
        "O=Cc1ccc(OCc2ccccc2)cc1",
        "4-benzyloxybenzaldehyde",
        "4-苄氧基苯甲醛",
    ),
    # 4-chlorobenzyloxy on phenol
    (
        "Oc1ccc(OCc2ccc(Cl)cc2)cc1",
        "4-(4-chlorobenzyloxy)phenol",
        "4-(4-氯苄氧基)苯酚",
    ),
    # pyridine + phenyl
    ("n1ccc(-c2ccccc2)cc1", "4-phenylpyridine", "4-苯基吡啶"),
    # pyridine + phenoxy
    ("n1ccc(Oc2ccccc2)cc1", "4-phenoxypyridine", "4-苯氧基吡啶"),
    # pyridine + benzyl
    ("n1ccc(Cc2ccccc2)cc1", "4-benzylpyridine", "4-苄基吡啶"),
    # negatives: retained alkoxy / phenoxybenzene
    ("COc1ccccc1", "anisole", "甲氧基苯"),
    ("CCOc1ccccc1", "ethoxybenzene", None),
    ("c1ccc(Oc2ccccc2)cc1", "phenoxybenzene", "苯氧基苯"),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_recursive_aryl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


def test_benzyloxy_not_ethoxy() -> None:
    r = SMILESNNamer().name("c1ccc(OCc2ccccc2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "ethoxy" not in en
    assert "heptane" not in en
    assert en == "benzyloxybenzene"


def test_phenylpyridine_not_nonane() -> None:
    r = SMILESNNamer().name("c1ccc(-c2ccncc2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "nonane" not in en
    assert "phenylpyridine" in en


def test_dichloro_n_phenyl_not_unsub() -> None:
    """Multi-halo N-Ph must not be treated as unsubstituted N-phenylamide."""
    r = SMILESNNamer().name("O=C(Nc1cc(Cl)cc(Cl)c1)C")
    assert r.success
    en = normalize_en(r.en)
    assert "N-phenylacetamide" not in en
    assert "N-phenyl" not in en


def test_mixed_halo_phenyl_alpha_order() -> None:
    """Mixed halos: cite by English name order (bromo before chloro)."""
    r = SMILESNNamer().name("O=Cc1ccc(-c2cc(Cl)cc(Br)c2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "5-bromo-3-chlorophenyl" in en
    assert "3-chloro-5-bromo" not in en
