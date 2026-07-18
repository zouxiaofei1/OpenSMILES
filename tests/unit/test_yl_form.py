from namepredict.layer3.yl_form import yl_form


def test_pyridine_to_yl():
    en, zh, paren = yl_form("pyridine", "吡啶", 2)
    assert en == "pyridin-2-yl"
    assert zh == "吡啶-2-基"
    assert paren is True


def test_thiazole_to_yl():
    en, zh, paren = yl_form("1,3-thiazole", "1,3-噻唑", 2)
    assert en == "1,3-thiazol-2-yl"
    assert "噻唑-2-基" in zh


def test_benzene_to_yl():
    en, zh, paren = yl_form("benzene", "苯", 1)
    # unsubstituted phenyl convention — locked phenyl
    assert en == "phenyl"
    assert zh == "苯基"
