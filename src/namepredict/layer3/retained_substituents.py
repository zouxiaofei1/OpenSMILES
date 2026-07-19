"""Centralized retained substituent registry (IUPAC 2013 P-29/P-57).

Provides dual-mode name resolution: "general" (retained names, default)
and "pin" (Preferred IUPAC Names).
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class IupacLevel(Enum):
    PIN = "pin"
    GENERAL = "general"
    NOT_RECOMMENDED = "not_rec"


@dataclass(frozen=True)
class RetainedSubstituent:
    en: str
    zh: str
    systematic_en: str
    systematic_zh: str
    level: IupacLevel
    rule_ref: str  # e.g. "P-29.6.1"
    # Optional: topology fingerprint for _try_registry_leaf matching.
    # leaf_atoms = sorted tuple of atomic nums in the claim (e.g. (7,8) for N+O)
    # leaf_root_z = required atomic num of claim.root (None = don't care)
    # leaf_bond_order = for 2-atom claims, bond order between root and other atom
    leaf_atoms: tuple[int, ...] | None = None
    leaf_root_z: int | None = None
    leaf_bond_order: float | None = None


def _build() -> dict[str, RetainedSubstituent]:
    P, G, N = IupacLevel.PIN, IupacLevel.GENERAL, IupacLevel.NOT_RECOMMENDED
    return {
        # ── 1. branched alkyls (P-29.3.1 / P-57.1) ──
        "tert-butyl": RetainedSubstituent(
            "tert-butyl", "叔丁基",
            "1,1-dimethylethyl", "1,1-二甲基乙基", P, "P-57.1.2",
        ),
        "isopropyl": RetainedSubstituent(
            "isopropyl", "异丙基",
            "propan-2-yl", "丙-2-基", G, "P-29.6.2.2",
        ),
        "isobutyl": RetainedSubstituent(
            "isobutyl", "异丁基",
            "2-methylpropyl", "2-甲基丙基", N, "P-57.1.4",
        ),
        "sec-butyl": RetainedSubstituent(
            "sec-butyl", "仲丁基",
            "butan-2-yl", "丁-2-基", N, "P-57.1.4",
        ),
        "neopentyl": RetainedSubstituent(
            "neopentyl", "新戊基",
            "2,2-dimethylpropyl", "2,2-二甲基丙基", N, "P-57.1.4",
        ),
        "isopentyl": RetainedSubstituent(
            "isopentyl", "异戊基",
            "3-methylbutyl", "3-甲基丁基", N, "P-57.1.4",
        ),
        "2-methylbutan-2-yl": RetainedSubstituent(
            "2-methylbutan-2-yl", "2-甲基丁-2-基",
            "2-methylbutan-2-yl", "2-甲基丁-2-基", P, "P-57.1.2",
        ),
        "3-methylbut-2-enyl": RetainedSubstituent(
            "3-methylbut-2-enyl", "3-甲基丁-2-烯基",
            "3-methylbut-2-enyl", "3-甲基丁-2-烯基", P, "P-29.6.1",
        ),
        "phenethyl": RetainedSubstituent(
            "phenethyl", "2-苯乙基",
            "2-phenylethyl", "2-苯乙基", G, "P-29.6.2.2",
        ),
        "benzhydryl": RetainedSubstituent(
            "benzhydryl", "二苯甲基",
            "diphenylmethyl", "二苯甲基", G, "P-29.6.2.2",
        ),
        "trityl": RetainedSubstituent(
            "trityl", "三苯甲基",
            "triphenylmethyl", "三苯甲基", G, "P-29.6.2.2",
        ),
        # ── 2. unsaturated acyclic (P-29.3.2 / P-29.6.1) ──
        "vinyl": RetainedSubstituent(
            "vinyl", "乙烯基",
            "ethenyl", "乙烯基", G, "P-29.6.2.2",
        ),
        "allyl": RetainedSubstituent(
            "allyl", "烯丙基",
            "prop-2-en-1-yl", "丙-2-烯-1-基", G, "P-29.6.2.2",
        ),
        "isopropenyl": RetainedSubstituent(
            "isopropenyl", "异丙烯基",
            "prop-1-en-2-yl", "丙-1-烯-2-基", G, "P-29.6.2.2",
        ),
        "propargyl": RetainedSubstituent(
            "propargyl", "炔丙基",
            "prop-2-yn-1-yl", "丙-2-炔-1-基", G, "P-29.6.2.2",
        ),
        "crotyl": RetainedSubstituent(
            "crotyl", "巴豆基",
            "(E)-but-2-enyl", "(E)-丁-2-烯基", G, "P-29.6.1",
        ),
        "cinnamyl": RetainedSubstituent(
            "cinnamyl", "肉桂基",
            "3-phenylprop-2-en-1-yl", "3-苯基丙-2-烯-1-基", G, "P-29.6.1",
        ),
        # ── 3. aryl / aralkyl (P-29.3.2 / P-29.6 / P-57.1.2) ──
        "phenyl": RetainedSubstituent(
            "phenyl", "苯基",
            "phenyl", "苯基", P, "P-29.6",
        ),
        "benzyl": RetainedSubstituent(
            "benzyl", "苄基",
            "phenylmethyl", "苯甲基", P, "P-57.1.2",
        ),
        "tolyl": RetainedSubstituent(
            "tolyl", "甲苯基",
            "methylphenyl", "甲基苯基", N, "P-57.1.5.4",
        ),
        "furfuryl": RetainedSubstituent(
            "furfuryl", "糠基",
            "furan-2-ylmethyl", "呋喃-2-基甲基", N, "P-57.1.5.4",
        ),
        "thenyl": RetainedSubstituent(
            "thenyl", "噻吩甲基",
            "thiophen-2-ylmethyl", "噻吩-2-基甲基", N, "P-57.1.5.4",
        ),
        # ── 3b. heteroaryl general-only (P-57.1.5.3) ──
        "furyl": RetainedSubstituent(
            "furyl", "呋喃基",
            "furan-2-yl", "呋喃-2-基", G, "P-57.1.5.3",
        ),
        "thienyl": RetainedSubstituent(
            "thienyl", "噻吩基",
            "thiophen-2-yl", "噻吩-2-基", G, "P-57.1.5.3",
        ),
        "pyridyl": RetainedSubstituent(
            "pyridyl", "吡啶基",
            "pyridin-2-yl", "吡啶-2-基", G, "P-57.1.5.3",
        ),
        "quinolyl": RetainedSubstituent(
            "quinolyl", "喹啉基",
            "quinolin-2-yl", "喹啉-2-基", G, "P-57.1.5.3",
        ),
        "isoquinolyl": RetainedSubstituent(
            "isoquinolyl", "异喹啉基",
            "isoquinolin-1-yl", "异喹啉-1-基", G, "P-57.1.5.3",
        ),
        "anthryl": RetainedSubstituent(
            "anthryl", "蒽基",
            "anthracen-2-yl", "蒽-2-基", G, "P-57.1.5.3",
        ),
        "phenanthryl": RetainedSubstituent(
            "phenanthryl", "菲基",
            "phenanthren-9-yl", "菲-9-基", G, "P-57.1.5.3",
        ),
        "adamantyl": RetainedSubstituent(
            "adamantyl", "金刚烷基",
            "adamantan-2-yl", "金刚烷-2-基", G, "P-57.1.5.3",
        ),
        "piperidyl": RetainedSubstituent(
            "piperidyl", "哌啶基",
            "piperidin-2-yl", "哌啶-2-基", G, "P-57.1.5.3",
        ),
        # ── 4. heteroatom prefixes — oxygen (P-31.1 / P-63) ──
        "methoxy": RetainedSubstituent(
            "methoxy", "甲氧基",
            "methoxy", "甲氧基", P, "P-63.2.2",
        ),
        "hydroperoxy": RetainedSubstituent(
            "hydroperoxy", "氢过氧基",
            "hydroperoxy", "氢过氧基", P, "P-63.5",
        ),
        # ── 5. heteroatom prefixes — sulfur (P-31.1 / P-63) ──
        "methylsulfanyl": RetainedSubstituent(
            "methylsulfanyl", "甲硫基",
            "methylsulfanyl", "甲硫基", P, "P-63.2.1",
        ),
        "ethylsulfanyl": RetainedSubstituent(
            "ethylsulfanyl", "乙硫基",
            "ethylsulfanyl", "乙硫基", P, "P-63.2.1",
        ),
        "sulfanyl": RetainedSubstituent(
            "sulfanyl", "硫烷基",
            "sulfanyl", "硫烷基", P, "P-63.1.5",
        ),
        "methylsulfinyl": RetainedSubstituent(
            "methylsulfinyl", "甲亚磺酰基",
            "methanesulfinyl", "甲亚磺酰基", P, "P-65.3.1",
        ),
        "methylsulfonyl": RetainedSubstituent(
            "methylsulfonyl", "甲磺酰基",
            "methanesulfonyl", "甲磺酰基", P, "P-65.3.2",
        ),
        "mesyl": RetainedSubstituent(
            "mesyl", "甲磺酰基",
            "methanesulfonyl", "甲磺酰基", N, "P-65.3.2",
        ),
        "tosyl": RetainedSubstituent(
            "tosyl", "对甲苯磺酰基",
            "4-methylbenzenesulfonyl", "4-甲基苯磺酰基", N, "P-65.3.2",
        ),
        "triflyl": RetainedSubstituent(
            "triflyl", "三氟甲磺酰基",
            "trifluoromethanesulfonyl", "三氟甲磺酰基", N, "P-65.3.2",
        ),
        # ── 6. heteroatom prefixes — nitrogen (P-31.1 / P-62) ──
        "nitroso": RetainedSubstituent(
            "nitroso", "亚硝基",
            "nitroso", "亚硝基", P, "P-61.5.2",
            leaf_atoms=(7, 8), leaf_root_z=7, leaf_bond_order=2.0,
        ),
        "azido": RetainedSubstituent(
            "azido", "叠氮基",
            "azido", "叠氮基", P, "P-61.8",
            leaf_atoms=(7, 7, 7),
        ),
        "hydrazinyl": RetainedSubstituent(
            "hydrazinyl", "肼基",
            "hydrazinyl", "肼基", P, "P-68.3.1.1",
        ),
        "anilino": RetainedSubstituent(
            "anilino", "苯胺基",
            "phenylamino", "苯氨基", P, "P-62.2.1.1.1",
        ),
        "methylamino": RetainedSubstituent(
            "methylamino", "甲氨基",
            "methylamino", "甲氨基", P, "P-62.2.1",
        ),
        "dimethylamino": RetainedSubstituent(
            "dimethylamino", "二甲氨基",
            "dimethylamino", "二甲氨基", P, "P-62.2.1",
        ),
        "diazenyl": RetainedSubstituent(
            "diazenyl", "二氮烯基",
            "diazenyl", "二氮烯基", P, "P-68.3.1.2",
        ),
        "diazo": RetainedSubstituent(
            "diazo", "重氮基",
            "diazo", "重氮基", P, "P-68.3.1.3",
        ),
        "cyano": RetainedSubstituent(
            "cyano", "氰基",
            "cyano", "氰基", P, "P-66.5.1.1.1",
            leaf_atoms=(6, 7), leaf_root_z=6, leaf_bond_order=3.0,
        ),
        "isocyano": RetainedSubstituent(
            "isocyano", "异氰基",
            "isocyano", "异氰基", P, "P-61.9",
            leaf_atoms=(6, 7), leaf_root_z=7, leaf_bond_order=2.0,
        ),
        "carbamoyl": RetainedSubstituent(
            "carbamoyl", "氨基甲酰基",
            "aminocarbonyl", "氨基羰基", P, "P-66.1.6.1",
        ),
        "sulfamoyl": RetainedSubstituent(
            "sulfamoyl", "氨基磺酰基",
            "aminosulfonyl", "氨基磺酰基", P, "P-65.3.2.1",
        ),
        "amidino": RetainedSubstituent(
            "amidino", "脒基",
            "carbaminidoyl", "甲脒基", N, "P-66.1.6.1.2",
        ),
        "guanidino": RetainedSubstituent(
            "guanidino", "胍基",
            "carbamimidamido", "胍基", N, "P-66.1.6.1.3",
        ),
        "ureido": RetainedSubstituent(
            "ureido", "脲基",
            "carbamoylamino", "脲基", N, "P-66.1.6.1.1",
        ),
        # ── 6b. P / B / Se prefixes (P-68.3.1) ──
        "phosphanyl": RetainedSubstituent(
            "phosphanyl", "膦基",
            "phosphanyl", "膦基", P, "P-68.3.1.4",
        ),
        "boranyl": RetainedSubstituent(
            "boranyl", "硼烷基",
            "boranyl", "硼烷基", P, "P-68.3.1.5",
        ),
        "selanyl": RetainedSubstituent(
            "selanyl", "硒烷基",
            "selanyl", "硒烷基", P, "P-68.3.1.6",
        ),
        # ── 7. halogens / haloalkyls (P-31.1 / P-61) ──
        "trichloromethyl": RetainedSubstituent(
            "trichloromethyl", "三氯甲基",
            "trichloromethyl", "三氯甲基", P, "P-61.5.1",
        ),
        "tribromomethyl": RetainedSubstituent(
            "tribromomethyl", "三溴甲基",
            "tribromomethyl", "三溴甲基", P, "P-61.5.1",
        ),
        "difluoromethyl": RetainedSubstituent(
            "difluoromethyl", "二氟甲基",
            "difluoromethyl", "二氟甲基", P, "P-61.5.1",
        ),
        "pentafluoroethyl": RetainedSubstituent(
            "pentafluoroethyl", "五氟乙基",
            "pentafluoroethyl", "五氟乙基", P, "P-61.5.1",
        ),
        # ── 8. acyl prefixes (P-65/P-66) ──
        "formyl": RetainedSubstituent(
            "formyl", "甲酰基",
            "formyl", "甲酰基", P, "P-65.1.3.1",
        ),
        "acetyl": RetainedSubstituent(
            "acetyl", "乙酰基",
            "acetyl", "乙酰基", P, "P-65.1.3.1",
        ),
        "benzoyl": RetainedSubstituent(
            "benzoyl", "苯甲酰基",
            "benzoyl", "苯甲酰基", P, "P-65.1.3.1",
        ),
        "propionyl": RetainedSubstituent(
            "propionyl", "丙酰基",
            "propanoyl", "丙酰基", G, "P-65.1.3.1",
        ),
        "butyryl": RetainedSubstituent(
            "butyryl", "丁酰基",
            "butanoyl", "丁酰基", G, "P-65.1.3.1",
        ),
        "isobutyryl": RetainedSubstituent(
            "isobutyryl", "异丁酰基",
            "2-methylpropanoyl", "2-甲基丙酰基", G, "P-65.1.3.1",
        ),
        "valeryl": RetainedSubstituent(
            "valeryl", "戊酰基",
            "pentanoyl", "戊酰基", G, "P-65.1.3.1",
        ),
        "oxamoyl": RetainedSubstituent(
            "oxamoyl", "草氨酰基",
            "oxamoyl", "草氨酰基", P, "P-66.1.6.1",
        ),
        "methoxycarbonyl": RetainedSubstituent(
            "methoxycarbonyl", "甲氧羰基",
            "methoxycarbonyl", "甲氧羰基", P, "P-65.6.3.1",
        ),
        "ethoxycarbonyl": RetainedSubstituent(
            "ethoxycarbonyl", "乙氧羰基",
            "ethoxycarbonyl", "乙氧羰基", P, "P-65.6.3.1",
        ),
    }


_REGISTRY: dict[str, RetainedSubstituent] = _build()


def get_retained(key: str) -> RetainedSubstituent | None:
    """Look up a retained substituent by registry key."""
    return _REGISTRY.get(key)


def resolve_name(key: str, *, name_mode: str = "general") -> tuple[str, str]:
    """Return (en, zh) for a registry key under the given naming mode.

    - "general": always the retained/common name
    - "pin": systematic name unless the entry itself is PIN-level
    """
    entry = _REGISTRY[key]
    if name_mode == "pin" and entry.level != IupacLevel.PIN:
        return entry.systematic_en, entry.systematic_zh
    return entry.en, entry.zh
