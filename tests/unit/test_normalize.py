from namepredict.tools.re import normalize_en, normalize_zh

def test_normalize_en_lower_and_space():
    assert normalize_en("  Ethanol  ") == "ethanol"
    assert normalize_en("propan-2-one") == "propan-2-one"
    assert normalize_en("A  B") == "a b"

def test_normalize_zh_strip_only():
    assert normalize_zh("  乙醇  ") == "乙醇"
    assert normalize_zh("乙醇") == "乙醇"
