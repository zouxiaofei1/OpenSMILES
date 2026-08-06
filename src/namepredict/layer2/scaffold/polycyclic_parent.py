"""Spiro + bridged (von Baeyer) parent hydride candidates (P-23 / P-24.2.1).

Round 1: saturated all-carbon systems with two monocyclic components only
(spiro) or bicyclic bridged systems only.  Hetero / unsat / polycyclic
components deferred to future rounds.
"""
from __future__ import annotations


# ── spiro (P-24.2.1) ──

def _is_simple_spiro(info: dict) -> dict | None:
    """Return the first all-carbon saturated mono+mono spiro system, or None."""
    for r in (info.get("ring_systems") or []):
        ok = (
            r.get("topology") == "spiro"
            and not r.get("hetero_atoms")
            and not r.get("is_aromatic_mancude")
            and r.get("n_rings") == 2
        )
        if ok:
            return r
    return None


def _spiro_stems(ring_sizes: list[int]) -> tuple[str, str]:
    """Spiro prefix stems: ('spiro[4.5]', '螺[4.5]')."""
    s = ".".join(str(x) for x in sorted(ring_sizes))
    return f"spiro[{s}]", f"螺[{s}]"


def _try_spiro_parent(info: dict) -> dict | None:
    """Recognize simple saturated all-carbon monospiro parent."""
    sp = _is_simple_spiro(info)
    if sp is None:
        return None
    from namepredict.layer2.parent_core import _parent_dict
    sizes = sp.get("ring_sizes") or []
    stem_en, stem_zh = _spiro_stems(sizes)
    return _parent_dict(
        sp["atom_ids"], "spiro",
        stem_en=stem_en, stem_zh=stem_zh, ring_sizes=sorted(sizes),
    )


# ── bridged / von Baeyer (P-23) ──

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
