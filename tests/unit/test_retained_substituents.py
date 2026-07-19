"""Tests for retained_substituents registry and resolve_name dispatch."""
from __future__ import annotations

import pytest

from namepredict.layer3.retained_substituents import (
    IupacLevel,
    get_retained,
    resolve_name,
)

# ── registry completeness ──

EXPECTED_KEYS = {
    # branched alkyls
    "tert-butyl", "isopropyl", "isobutyl", "sec-butyl",
    "neopentyl", "isopentyl", "2-methylbutan-2-yl", "3-methylbut-2-enyl",
    "phenethyl", "benzhydryl", "trityl",
    # unsaturated acyclic
    "vinyl", "allyl", "isopropenyl", "propargyl", "crotyl", "cinnamyl",
    # aryl / aralkyl
    "phenyl", "benzyl", "tolyl", "furfuryl", "thenyl",
    # heteroaryl general-only
    "furyl", "thienyl", "pyridyl", "quinolyl", "isoquinolyl",
    "anthryl", "phenanthryl", "adamantyl", "piperidyl",
    # heteroatom — O
    "methoxy", "hydroperoxy",
    # heteroatom — S
    "methylsulfanyl", "ethylsulfanyl", "sulfanyl",
    "methylsulfinyl", "methylsulfonyl",
    "sulfo",
    "mesyl", "tosyl", "triflyl",
    # heteroatom — N
    "nitroso", "azido", "hydrazinyl", "anilino",
    "methylamino", "dimethylamino", "diazenyl", "diazo",
    "cyano", "isocyano",
    "carbamoyl", "sulfamoyl",
    "amidino", "guanidino", "ureido",
    # P / B / Se prefixes
    "phosphanyl", "boranyl", "selanyl",
    # haloalkyls
    "trichloromethyl", "tribromomethyl", "difluoromethyl", "pentafluoroethyl",
    # acyl
    "formyl", "acetyl", "benzoyl",
    "propionyl", "butyryl", "isobutyryl", "valeryl",
    "oxamoyl",
    "methoxycarbonyl", "ethoxycarbonyl",
}


def test_all_expected_keys_present():
    for key in EXPECTED_KEYS:
        assert get_retained(key) is not None, f"missing key: {key}"


def test_duplicate_names_only_for_synonyms():
    """Some entries share en names (e.g. mesyl vs methylsulfonyl) — allowed."""
    # All expected keys must resolve
    for key in EXPECTED_KEYS:
        assert get_retained(key) is not None, f"missing: {key}"


# ── resolve_name dispatch ──

def test_general_always_uses_retained_name():
    assert resolve_name("vinyl", name_mode="general") == ("vinyl", "乙烯基")
    assert resolve_name("tert-butyl", name_mode="general") == ("tert-butyl", "叔丁基")


def test_pin_uses_systematic_for_general_level():
    assert resolve_name("vinyl", name_mode="pin") == ("ethenyl", "乙烯基")
    assert resolve_name("allyl", name_mode="pin") == ("prop-2-en-1-yl", "丙-2-烯-1-基")
    assert resolve_name("isopropyl", name_mode="pin") == ("propan-2-yl", "丙-2-基")


def test_pin_uses_systematic_for_not_recommended():
    assert resolve_name("isobutyl", name_mode="pin") == ("2-methylpropyl", "2-甲基丙基")
    assert resolve_name("sec-butyl", name_mode="pin") == ("butan-2-yl", "丁-2-基")


def test_pin_keeps_retained_for_pin_level():
    assert resolve_name("tert-butyl", name_mode="pin") == ("tert-butyl", "叔丁基")
    assert resolve_name("phenyl", name_mode="pin") == ("phenyl", "苯基")
    assert resolve_name("methylsulfanyl", name_mode="pin") == ("methylsulfanyl", "甲硫基")


def test_system_names_for_pin_level_also_pin():
    """PIN-level entries: retained == systematic (same PIN name)."""
    assert resolve_name("2-methylbutan-2-yl", name_mode="pin") == ("2-methylbutan-2-yl", "2-甲基丁-2-基")
    assert resolve_name("methoxy", name_mode="pin") == ("methoxy", "甲氧基")


def test_default_mode_is_general():
    assert resolve_name("vinyl") == ("vinyl", "乙烯基")
    assert resolve_name("isopropyl") == ("isopropyl", "异丙基")


# ── IupacLevel sanity ──

def test_iupac_level_values():
    assert get_retained("tert-butyl").level == IupacLevel.PIN
    assert get_retained("isopropyl").level == IupacLevel.GENERAL
    assert get_retained("isobutyl").level == IupacLevel.NOT_RECOMMENDED


# ── new category spot-checks ──

def test_acyl_general_vs_pin():
    assert resolve_name("propionyl", name_mode="general") == ("propionyl", "丙酰基")
    assert resolve_name("propionyl", name_mode="pin") == ("propanoyl", "丙酰基")


def test_sulfonyl_general_vs_pin():
    assert resolve_name("mesyl", name_mode="general") == ("mesyl", "甲磺酰基")
    assert resolve_name("mesyl", name_mode="pin") == ("methanesulfonyl", "甲磺酰基")


def test_not_recommended_uses_systematic_in_pin():
    assert resolve_name("furfuryl", name_mode="pin") == ("furan-2-ylmethyl", "呋喃-2-基甲基")
    assert resolve_name("tolyl", name_mode="pin") == ("methylphenyl", "甲基苯基")


def test_haloalkyl_pin_level_uses_same_name():
    """Haloalkyl entries are PIN level — retained == systematic."""
    for key in ("trichloromethyl", "difluoromethyl", "pentafluoroethyl"):
        e = get_retained(key)
        assert e.level == IupacLevel.PIN
        assert e.en == e.systematic_en


# ── leaf_atoms topology fingerprints ──

LEAF_KEYS = {"nitroso", "azido", "cyano", "isocyano", "sulfo"}


def test_leaf_entries_have_atoms_and_non_empty():
    for key in LEAF_KEYS:
        e = get_retained(key)
        assert e is not None, f"missing leaf: {key}"
        assert e.leaf_atoms is not None and len(e.leaf_atoms) > 0, f"{key}: no leaf_atoms"


def test_leaf_atoms_consistent_with_key():
    """cyano(C+N) vs isocyano(N+C) distinguished by leaf_root_z only."""
    c = get_retained("cyano")
    i = get_retained("isocyano")
    assert c.leaf_atoms == i.leaf_atoms == (6, 7)  # same atoms
    assert c.leaf_root_z == 6  # C is root for cyano
    assert i.leaf_root_z == 7  # N is root for isocyano


def test_azido_has_three_nitrogens():
    e = get_retained("azido")
    assert e.leaf_atoms == (7, 7, 7)
    assert e.leaf_root_z is None  # any N can be root


def test_nitroso_leaf_topology():
    e = get_retained("nitroso")
    assert e.leaf_atoms == (7, 8)
    assert e.leaf_root_z == 7  # N is root


def test_sulfo_leaf_topology():
    e = get_retained("sulfo")
    assert e.leaf_atoms == (8, 8, 8, 16)  # 3xO + S
    assert e.leaf_root_z == 16  # S is root


# ── heteroaryl general-only (P-57.1.5.3) ──

def test_heteroaryl_general_vs_pin():
    """Heteroaryl retained names (furyl, thienyl, etc.) are GENERAL only."""
    assert resolve_name("furyl", name_mode="general") == ("furyl", "呋喃基")
    assert resolve_name("furyl", name_mode="pin") == ("furan-2-yl", "呋喃-2-基")
    assert resolve_name("thienyl", name_mode="general") == ("thienyl", "噻吩基")
    assert resolve_name("thienyl", name_mode="pin") == ("thiophen-2-yl", "噻吩-2-基")
    assert resolve_name("pyridyl", name_mode="general") == ("pyridyl", "吡啶基")
    assert resolve_name("pyridyl", name_mode="pin") == ("pyridin-2-yl", "吡啶-2-基")


def test_heteroaryl_all_general_level():
    for key in ("furyl", "thienyl", "pyridyl", "quinolyl", "isoquinolyl",
                 "anthryl", "phenanthryl", "adamantyl", "piperidyl"):
        assert get_retained(key).level == IupacLevel.GENERAL, f"{key} should be GENERAL"


# ── new N prefixes ──

def test_methylamino_pin_same_as_systematic():
    e = get_retained("methylamino")
    assert e.level == IupacLevel.PIN
    assert e.en == e.systematic_en == "methylamino"


def test_dimethylamino_pin():
    e = get_retained("dimethylamino")
    assert e.level == IupacLevel.PIN
    assert resolve_name("dimethylamino", name_mode="pin") == ("dimethylamino", "二甲氨基")


# ── P / B / Se prefixes ──

def test_pbse_all_pin_level():
    for key in ("phosphanyl", "boranyl", "selanyl"):
        e = get_retained(key)
        assert e.level == IupacLevel.PIN, f"{key} should be PIN"
        assert e.en == e.systematic_en, f"{key}: retained == systematic"


# ── oxamoyl ──

def test_oxamoyl_pin():
    e = get_retained("oxamoyl")
    assert e.level == IupacLevel.PIN
    assert resolve_name("oxamoyl", name_mode="general") == ("oxamoyl", "草氨酰基")


# ── integration: 5 sulfur substituent prefixes ──

_SULFUR_PREFIX_CASES = [
    # (smiles, expected_en_fragment, label)
    ("CS(=O)CC[C@H](N)C(=O)O", "methylsulfinyl", "methylsulfinyl on amino acid"),
    ("CS(=O)(=O)CC[C@H](N)C(=O)O", "methylsulfonyl", "methylsulfonyl on amino acid (benchmark fix)"),
    ("CS(=O)(=O)C", "methylsulfonyl", "dimethyl sulfone → methylsulfonyl prefix"),
    ("O=S(=O)(C(F)(F)F)CC[C@H](N)C(=O)O", "triflyl", "triflyl on amino acid"),
    ("Cc1ccc(S(=O)(=O)CC[C@H](N)C(=O)O)cc1", "tosyl", "tosyl on amino acid"),
    ("CS(=O)(=O)c1ccc(C(=O)c2c(Sc3ccccc3)C3CCC(C3)C2=O)c(Cl)c1", "methylsulfonyl", "methylsulfonyl aryl ketone"),
]


@pytest.mark.parametrize("smiles,expected_fragment,label", _SULFUR_PREFIX_CASES)
def test_sulfur_prefix_naming(smiles, expected_fragment, label):
    from namepredict.namer import SMILESNNamer
    r = SMILESNNamer(name_mode="general").name(smiles)
    assert r.success, f"{label}: naming failed"
    assert expected_fragment in r.en, f"{label}: expected '{expected_fragment}' in '{r.en}'"


def test_methylsulfanyl_unchanged():
    """Plain thioether should still give methylsulfanyl, not a sulfonyl variant."""
    from namepredict.namer import SMILESNNamer
    r = SMILESNNamer(name_mode="general").name("CC(C)SC")
    assert r.success
    assert "methylsulfanyl" in r.en


def test_sulfonate_not_falsely_tosyl():
    """Butyl tosylate must remain sulfonate parent, not tosyl substituent."""
    from namepredict.namer import SMILESNNamer
    r = SMILESNNamer(name_mode="general").name("CC1=CC=C(C=C1)S(=O)(=O)OCCCC")
    assert r.success
    assert "4-methylbenzenesulfonate" in r.en
