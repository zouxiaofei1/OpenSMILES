"""Bridged (von Baeyer) parent hydride candidate (IUPAC P-23).

Round 1: saturated all-carbon bicyclic systems only.
Hetero / unsat / polycyclic-bridged deferred to future rounds.
"""
from __future__ import annotations


def _is_bridged_system(r: dict) -> bool:
    return (r.get("topology") == "bridged"
            and not r.get("hetero_atoms")
            and not r.get("is_aromatic_mancude")
            and r.get("n_rings") == 2)


def _is_simple_bridged(info: dict) -> dict | None:
    """Return the first saturated all-carbon bridged bicyclic system, or None."""
    for r in (info.get("ring_systems") or []):
        if _is_bridged_system(r):
            return r
    return None


def _bicyclo_stems(bridge_lengths: list[int]) -> tuple[str, str]:
    """Bicyclo prefix stems: ('bicyclo[2.2.1]', '双环[2.2.1]')."""
    s = ".".join(str(x) for x in sorted(bridge_lengths, reverse=True))
    return f"bicyclo[{s}]", f"双环[{s}]"


def _build_bridged_parent(br: dict) -> dict:
    bl = br.get("bridge_lengths") or []
    stem_en, stem_zh = _bicyclo_stems(bl)
    from namepredict.layer2.parent_core import _parent_dict
    return _parent_dict(br["atom_ids"], "bridged", stem_en=stem_en, stem_zh=stem_zh,
                        bridge_lengths=sorted(bl, reverse=True),
                        bridgeheads=br.get("bridgeheads"),
                        bridge_paths=br.get("bridge_paths"))


def _try_bridged_parent(info: dict) -> dict | None:
    """Recognize simple saturated all-carbon bicyclic bridged parent."""
    br = _is_simple_bridged(info)
    return None if br is None else _build_bridged_parent(br)
