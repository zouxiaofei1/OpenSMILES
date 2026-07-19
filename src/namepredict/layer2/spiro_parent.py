"""Spiro parent hydride candidate (IUPAC P-24.2.1).

Round 1: saturated all-carbon monospiro with two monocyclic components only.
Hetero / unsat / polycyclic-component spiros deferred to future rounds.
"""
from __future__ import annotations


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
