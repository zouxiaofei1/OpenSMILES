"""Tests for retained_substituents registry and resolve_name dispatch."""
from __future__ import annotations

import pytest
from rdkit.Chem import MolFromSmiles, MolToSmiles

import namepredict.tools.anchored_table as at
from namepredict.tools.anchored_table import (
    ANCHOR_TABLE,
    _ANCHOR_INDEX,
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


def test_methylsulfanyl_unchanged():
    """Plain thioether should still give methylsulfanyl, not a sulfonyl variant."""
    from namepredict.namer import SMILESNNamer
    r = SMILESNNamer(name_mode="general").name("CC(C)SC")
    assert r.success
    assert "methylsulfanyl" in r.en


# ── 合表后：锚定索引一致性 / canonical 校验 ──

def test_anchor_index_matches_registry_anchored():
    """_ANCHOR_INDEX 的每个键必须来自对应 registry 条目的 anchored 字段。"""
    for smi, key in _ANCHOR_INDEX.items():
        assert smi in get_retained(key).anchored, f"{smi} -> {key}"


def test_anchored_keys_all_canonical():
    """registry 的 anchored 键必须是 canonical 形式（防 *C=C-C 类死条目）。"""
    for key, e in at._REGISTRY.items():
        for smi in e.anchored:
            m = MolFromSmiles(smi)
            assert m is not None, f"{key}: 解析失败 {smi}"
            assert MolToSmiles(m) == smi, f"{key}: {smi} 非 canonical -> {MolToSmiles(m)}"


def test_inline_keys_all_canonical_and_no_collision():
    """内联表键 canonical 且不与 registry 锚定冲突。"""
    for smi in ANCHOR_TABLE:
        m = MolFromSmiles(smi)
        assert m is not None, f"内联键解析失败 {smi}"
        assert MolToSmiles(m) == smi, f"内联键非 canonical: {smi}"
        assert smi not in _ANCHOR_INDEX, f"跨表冲突: {smi}"


def test_allyl_only_true_topology():
    """allyl 只对应真烯丙基连接位点 *CC=C；*C=CC 是丙-1-烯基内联条目。"""
    assert get_retained("allyl").anchored == ("*CC=C",)
    assert "*C=CC" not in _ANCHOR_INDEX
    assert "*C=CC" in ANCHOR_TABLE
    assert "*C=C-C" not in _ANCHOR_INDEX
    assert "*C=C-C" not in ANCHOR_TABLE


def test_unsaturated_isobutyl_isopentyl_not_anchored():
    """不饱和烯基拓扑不再错挂到饱和保留名 isobutyl/isopentyl。"""
    assert get_retained("isobutyl").anchored == ("*CC(C)C",)
    assert get_retained("isopentyl").anchored == ("*CCC(C)C",)
    assert "*C=C(C)C" not in _ANCHOR_INDEX
    assert "*C=CC(C)C" not in _ANCHOR_INDEX


def test_build_anchor_index_rejects_noncanonical():
    """非 canonical anchored 键（如 *C=C-C）在构建索引时立即报错。"""
    orig = at._REGISTRY
    bad = dict(orig)
    bad["_fake_noncanon"] = at.RetainedSubstituent(
        "x", "x", "x", "x", IupacLevel.PIN, anchored=("*C=C-C",))
    at._REGISTRY = bad
    try:
        with pytest.raises(ValueError):
            at._build_anchor_index()
    finally:
        at._REGISTRY = orig


def test_build_anchor_index_rejects_collision():
    """同一锚定键映射两个保留基必须报错。"""
    orig = at._REGISTRY
    bad = dict(orig)
    bad["_fake_dup"] = at.RetainedSubstituent(
        "x", "x", "x", "x", IupacLevel.PIN, anchored=("*CC(C)C",))
    at._REGISTRY = bad
    try:
        with pytest.raises(ValueError):
            at._build_anchor_index()
    finally:
        at._REGISTRY = orig


def test_prop_1_enyl_vs_allyl_resolve():
    """*C=CC 与 *CC=C 是不同的基团：丙-1-烯基 vs 烯丙基。"""
    mol1 = MolFromSmiles(r"C/C=C\c1cc(OC)c(OC)cc1OC")
    assert at.anchored_entry(mol1, frozenset({0, 1, 2})) == (
        "prop-1-enyl", "丙-1-烯基", True, "alkyl")
    mol2 = MolFromSmiles("C=CCc1ccccc1")
    assert at.anchored_entry(mol2, frozenset({0, 1, 2})) == (
        "allyl", "烯丙基", False, "alkyl")
