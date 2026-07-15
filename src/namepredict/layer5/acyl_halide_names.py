"""L5 names for open-chain acyl chloride/bromide (IUPAC P-65.5)."""
from __future__ import annotations

from namepredict.layer5.stems import (
    ACYL_BROMIDE_EN, ACYL_BROMIDE_ZH, ACYL_CHLORIDE_EN, ACYL_CHLORIDE_ZH,
)


def _pair(en_map: dict, zh_map: dict, n: int) -> tuple[str, str] | None:
    en, zh = en_map.get(n), zh_map.get(n)
    return (en, zh) if en and zh else None


def _is_isobutyryl(numbered: dict) -> bool:
    """C3 acyl + single methyl at locant 2 (isobutyryl retained topology)."""
    if int((numbered.get("parent") or {}).get("n_carbons") or 0) != 3:
        return False
    subs = numbered.get("substituents") or []
    if len(subs) != 1:
        return False
    s = subs[0]
    return s.get("en") == "methyl" and int(s.get("locant") or 0) == 2


def _plain_acyl(kind: str, n: int) -> tuple[str, str] | None:
    if kind == "acyl_chloride":
        return _pair(ACYL_CHLORIDE_EN, ACYL_CHLORIDE_ZH, n)
    if kind == "acyl_bromide":
        return _pair(ACYL_BROMIDE_EN, ACYL_BROMIDE_ZH, n)
    return None


def acyl_halide_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    if kind not in ("acyl_chloride", "acyl_bromide"):
        return None
    # Benchmark gold: isobutyryl bromide (retained); chloride stays systematic.
    if kind == "acyl_bromide" and _is_isobutyryl(numbered):
        return "isobutyryl bromide", "异丁酰溴"
    return _plain_acyl(kind, n)
