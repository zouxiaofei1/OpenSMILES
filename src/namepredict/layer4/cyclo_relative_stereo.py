"""L4 与位次对齐的环多酸相对立体化学事实。"""
from __future__ import annotations


def _ordered(oriented: dict):
    """按位次标签排序立体面列表；缺任一项则空列表。"""
    stereo, plan = oriented.get("relative_stereo"), oriented.get("numbering")
    if stereo is None or plan is None or not stereo.faces: return []
    return sorted(stereo.faces, key=lambda item: int(plan.atom_to_label[item[0]]))


def _locants(ordered, plan) -> str:
    """生成三位立体前缀（r/c/t）的 locant 字符串。"""
    reference = ordered[0][1]
    return ",".join(f"{plan.atom_to_label[atom]}-{'r' if i == 0 else ('c' if face == reference else 't')}" for i, (atom, face) in enumerate(ordered))


def relative_stereo_facts(oriented: dict) -> dict:
    """仅在 NumberingPlan 选定位次后生成名称。"""
    if oriented.get("kind") != "acid" or oriented.get("scaffold_id") != "carbocycle":
        return {}
    ordered = _ordered(oriented)
    if len(ordered) == 2:
        return {"relative_stereo_prefix": "cis" if ordered[0][1] == ordered[1][1] else "trans"}
    if len(ordered) == 3: return {"relative_stereo_locants": _locants(ordered, oriented["numbering"])}
    return {}
