"""Verify general vs pin mode name switching for branched alkyls."""
from __future__ import annotations

import pytest

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en


CASES = [
    # (smiles, general_en_substring, pin_en_substring)
    # Must use a parent (benzene) so branched alkyl is extracted as substituent.
    ("CC(C)(C)c1ccccc1", "tert-butyl", "tert-butyl"),  # PIN level: same in both
]


@pytest.mark.parametrize("smiles,gen_token,pin_token", CASES)
def test_mode_switching(smiles, gen_token, pin_token):
    r_gen = SMILESNNamer().name(smiles)
    r_pin = SMILESNNamer(name_mode="pin").name(smiles)
    assert r_gen.success and r_pin.success
    assert gen_token in normalize_en(r_gen.en)
    assert pin_token in normalize_en(r_pin.en)


def test_default_is_general():
    """SMILESNNamer() == SMILESNNamer(name_mode='general')."""
    r_default = SMILESNNamer().name("CC(C)CC")
    r_general = SMILESNNamer(name_mode="general").name("CC(C)CC")
    assert normalize_en(r_default.en) == normalize_en(r_general.en)


def test_isopropyl_pin_uses_propan_2_yl():
    r = SMILESNNamer(name_mode="pin").name("CC(C)c1ccccc1")
    assert r.success
    assert "propan-2-yl" in normalize_en(r.en)
    assert "isopropyl" not in normalize_en(r.en)


def test_tert_butyl_unchanged_in_pin():
    """tert-butyl IS a PIN — no change in pin mode."""
    r_gen = SMILESNNamer().name("CC(C)(C)c1ccccc1")
    r_pin = SMILESNNamer(name_mode="pin").name("CC(C)(C)c1ccccc1")
    assert "tert-butyl" in normalize_en(r_gen.en)
    assert "tert-butyl" in normalize_en(r_pin.en)
