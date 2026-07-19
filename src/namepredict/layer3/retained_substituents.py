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


def _build() -> dict[str, RetainedSubstituent]:
    P, G, N = IupacLevel.PIN, IupacLevel.GENERAL, IupacLevel.NOT_RECOMMENDED
    return {
        # ── branched alkyls (P-29.3.1 / P-57.1) ──
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
        # ── unsaturated acyclic (P-29.3.2 / P-29.6.1) ──
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
        # ── aryl (P-29.3.2 / P-29.6 / P-57.1.2) ──
        "phenyl": RetainedSubstituent(
            "phenyl", "苯基",
            "phenyl", "苯基", P, "P-29.6",
        ),
        "benzyl": RetainedSubstituent(
            "benzyl", "苄基",
            "phenylmethyl", "苯甲基", P, "P-57.1.2",
        ),
        # ── heteroatom prefixes (P-31.1 / P-63) ──
        "methoxy": RetainedSubstituent(
            "methoxy", "甲氧基",
            "methoxy", "甲氧基", P, "P-63.2.2",
        ),
        "methylsulfanyl": RetainedSubstituent(
            "methylsulfanyl", "甲硫基",
            "methylsulfanyl", "甲硫基", P, "P-63.2.1",
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
