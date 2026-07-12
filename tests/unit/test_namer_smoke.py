from namepredict.namer import SMILESNNamer
from namepredict.constants import normalize_en, normalize_zh


def test_methane():
    r = SMILESNNamer().name("C")
    assert r.success
    assert normalize_en(r.en) == "methane"
    assert normalize_zh(r.zh) == "甲烷"


def test_ethanol():
    r = SMILESNNamer().name("CCO")
    assert r.success
    assert normalize_en(r.en) == "ethanol"
    assert normalize_zh(r.zh) == "乙醇"
