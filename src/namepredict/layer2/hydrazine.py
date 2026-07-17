"""L2 simple hydrazine functional parent (P-68.3.1.2)."""
from __future__ import annotations


_HZ_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs
    return len(info.get("hydrazines") or []) == 1 and _no_fgs(info, _HZ_BAD)


def _arm_n(mol, c_idx: int, n_idx: int) -> int:
    from namepredict.layer2.parent_core import _longest_from
    return len(_longest_from(mol, c_idx, {n_idx}) or [c_idx])


def _is_plain_phenyl(mol, c_idx: int, n_idx: int) -> bool:
    from namepredict.layer2.aryl_sub import _count_side_leaves, _phenyl_at
    ph = _phenyl_at(mol, c_idx, n_idx)
    if ph is None:
        return False
    return _count_side_leaves(mol, ph, c_idx, n_idx) == 0


def _side_mono(mol, n_idx: int, c: int) -> dict:
    if mol.GetAtomWithIdx(c).GetIsAromatic():
        if _is_plain_phenyl(mol, c, n_idx):
            return {"kind": "phenyl", "n_idx": n_idx, "aryl_c": c}
        return {"kind": "complex", "n_idx": n_idx}
    n = _arm_n(mol, c, n_idx)
    if 1 <= n <= 4:
        return {"kind": "alkyl", "n_idx": n_idx, "ns": [n]}
    return {"kind": "complex", "n_idx": n_idx}


def _side_di(mol, n_idx: int, cs: list[int]) -> dict:
    if any(mol.GetAtomWithIdx(c).GetIsAromatic() for c in cs):
        return {"kind": "complex", "n_idx": n_idx}
    ns = [_arm_n(mol, c, n_idx) for c in cs]
    if len(cs) == 2 and all(1 <= x <= 4 for x in ns):
        return {"kind": "dialkyl", "n_idx": n_idx, "ns": ns}
    return {"kind": "complex", "n_idx": n_idx}


def _side(mol, n_idx: int, cs: list[int]) -> dict:
    if not cs:
        return {"kind": "h", "n_idx": n_idx}
    return _side_mono(mol, n_idx, cs[0]) if len(cs) == 1 else _side_di(mol, n_idx, cs)


def _simple_sides(s1: dict, s2: dict) -> bool:
    kinds = {s1["kind"], s2["kind"]}
    if kinds in ({"h"}, {"phenyl", "h"}):
        return True
    return kinds == {"dialkyl", "h"}


def _order_sides(s1: dict, s2: dict) -> tuple[dict, dict]:
    if s1["kind"] == "h" and s2["kind"] != "h":
        return s2, s1
    return s1, s2


def _sides_ok(info: dict, e: dict) -> tuple[dict, dict] | None:
    mol = info["mol"]
    s1 = _side(mol, e["n1_idx"], list(e.get("n1_c_idxs") or []))
    s2 = _side(mol, e["n2_idx"], list(e.get("n2_c_idxs") or []))
    return _order_sides(s1, s2) if _simple_sides(s1, s2) else None


def _pack_parent(info: dict, e: dict, n1: dict, n2: dict) -> dict:
    from namepredict.layer2.parent_core import _parent_dict
    # Empty carbon chain: N-subs live in n1/n2; avoid L3 side extraction.
    return _parent_dict(
        [], "hydrazine",
        n1_idx=e["n1_idx"], n2_idx=e["n2_idx"],
        mol=info["mol"], n1=n1, n2=n2,
    )


def _hydrazine_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e = info["hydrazines"][0]
    ordered = _sides_ok(info, e)
    return None if ordered is None else _pack_parent(info, e, *ordered)
