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
    # ortho-halo phenyl lowest locant is 2- (not 6-) when Ph is the arm
    (
        "O=Cc1ccc(-c2ccccc2Cl)cc1",
        "4-(2-chlorophenyl)benzaldehyde",
        "4-(2-氯苯基)苯甲醛",
    ),
    # benzylbenzene
    ("c1ccc(Cc2ccccc2)cc1", "benzylbenzene", "苄基苯"),
    # pyridine + phenyl
    ("n1ccc(-c2ccccc2)cc1", "4-phenylpyridine", "4-苯基吡啶"),
    # pyridine + benzyl
    ("n1ccc(Cc2ccccc2)cc1", "4-benzylpyridine", "4-苄基吡啶"),
    # negatives: retained alkoxy / phenoxybenzene
    ("CCOc1ccccc1", "ethoxybenzene", None),
]


@pytest.mark.parametrize("smiles,en,zh", CASES)
def test_recursive_aryl(smiles: str, en: str, zh: str | None) -> None:
    r = SMILESNNamer().name(smiles)
    assert r.success
    assert normalize_en(r.en) == normalize_en(en)
    if zh is not None:
        assert normalize_zh(r.zh) == normalize_zh(zh)


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
    """Mixed halos: cite by English name order (bromo before chloro); the
    alpha-earlier substituent gets the lower locant (P-14.4)."""
    r = SMILESNNamer().name("O=Cc1ccc(-c2cc(Cl)cc(Br)c2)cc1")
    assert r.success
    en = normalize_en(r.en)
    assert "3-bromo-5-chlorophenyl" in en
    assert "5-bromo-3-chlorophenyl" not in en
