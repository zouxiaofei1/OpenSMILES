"""Pure EN/ZH skeletal replacement stems for saturated monocycles.

Input: ring_size + pre-numbered hetero sites [(locant, Z), ...].
No RDKit / mol / namer. Unknown Z or size → None.
"""
from __future__ import annotations

from namepredict.layer2.scaffold.hetero_a_names import (
    MULT_EN,
    MULT_ZH,
    a_prefix_en,
    a_root_zh,
    cyclo_ane_en,
    cyclo_ane_zh,
)


def _sorted_sites(
    sites: list[tuple[int, int]],
) -> list[tuple[int, int]] | None:
    if not sites:
        return None
    return sorted(sites, key=lambda s: s[0])


def _all_known(sites: list[tuple[int, int]]) -> bool:
    return all(a_prefix_en(z) is not None for _, z in sites)


def _same_z(sites: list[tuple[int, int]]) -> bool:
    return len({z for _, z in sites}) == 1


def _locants(sites: list[tuple[int, int]]) -> str:
    return ",".join(str(loc) for loc, _ in sites)


def _en_same_chunk(sites: list[tuple[int, int]]) -> str | None:
    z = sites[0][1]
    pref = a_prefix_en(z)
    mult = MULT_EN.get(len(sites))
    if pref is None or mult is None:
        return None
    return f"{_locants(sites)}-{mult}{pref}"


def _en_mixed_chunk(sites: list[tuple[int, int]]) -> str | None:
    parts: list[str] = []
    for loc, z in sites:
        pref = a_prefix_en(z)
        if pref is None:
            return None
        parts.append(f"{loc}-{pref}")
    return "-".join(parts)


def _en_hetero_chunk(sites: list[tuple[int, int]]) -> str | None:
    if _same_z(sites):
        return _en_same_chunk(sites)
    return _en_mixed_chunk(sites)


def sat_hetero_stem_en(
    ring_size: int, sites: list[tuple[int, int]]
) -> str | None:
    ordered = _sorted_sites(sites)
    cyclo = cyclo_ane_en(ring_size)
    if ordered is None or cyclo is None or not _all_known(ordered):
        return None
    chunk = _en_hetero_chunk(ordered)
    return f"{chunk}{cyclo}" if chunk else None


def _zh_roots(sites: list[tuple[int, int]]) -> list[str] | None:
    roots = [a_root_zh(z) for _, z in sites]
    if any(r is None for r in roots):
        return None
    return [r for r in roots if r is not None]


def _zh_hetero_body(sites: list[tuple[int, int]]) -> str | None:
    roots = _zh_roots(sites)
    if roots is None:
        return None
    if _same_z(sites):
        mult = MULT_ZH.get(len(sites))
        return None if mult is None else f"{mult}{roots[0]}杂"
    return f"{''.join(roots)}杂"


def sat_hetero_stem_zh(
    ring_size: int, sites: list[tuple[int, int]]
) -> str | None:
    ordered = _sorted_sites(sites)
    cyclo = cyclo_ane_zh(ring_size)
    if ordered is None or cyclo is None or not _all_known(ordered):
        return None
    body = _zh_hetero_body(ordered)
    if body is None:
        return None
    if len(ordered) == 1:
        return f"{body}{cyclo}"
    return f"{_locants(ordered)}-{body}{cyclo}"
