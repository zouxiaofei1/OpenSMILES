"""L4 locant-aligned relative-stereochemistry facts for cyclo polyacids."""
from __future__ import annotations


def _ordered(oriented: dict):
    stereo, plan = oriented.get("relative_stereo"), oriented.get("numbering")
    if stereo is None or plan is None or not stereo.faces: return []
    return sorted(stereo.faces, key=lambda item: int(plan.atom_to_label[item[0]]))


def _locants(ordered, plan) -> str:
    reference = ordered[0][1]
    return ",".join(f"{plan.atom_to_label[atom]}{'r' if i == 0 else ('c' if face == reference else 't')}" for i, (atom, face) in enumerate(ordered))


def relative_stereo_facts(oriented: dict) -> dict:
    """Materialize names only after NumberingPlan has selected locants."""
    if oriented.get("kind") != "cycloalkane_polycarboxylic": return {}
    ordered = _ordered(oriented)
    if len(ordered) == 2:
        return {"relative_stereo_prefix": "cis" if ordered[0][1] == ordered[1][1] else "trans"}
    if len(ordered) == 3: return {"relative_stereo_locants": _locants(ordered, oriented["numbering"])}
    return {}
