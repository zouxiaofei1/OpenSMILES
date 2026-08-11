"""L5 names for open-chain acyl chloride/bromide (IUPAC P-65.5)."""
from __future__ import annotations


def _is_isobutyryl(numbered: dict) -> bool:
    """C3 acyl + single methyl at locant 2 (isobutyryl retained topology)."""
    if int((numbered.get("parent") or {}).get("n_carbons") or 0) != 3:
        return False
    subs = numbered.get("substituents") or []
    if len(subs) != 1:
        return False
    s = subs[0]
    return s.get("en") == "methyl" and int(s.get("locant") or 0) == 2
