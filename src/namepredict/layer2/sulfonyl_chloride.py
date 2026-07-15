"""L2 simple mono-sulfonyl chloride functional parent (P-65.3)."""
from __future__ import annotations


_SC_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
    "has_sulfoxide", "has_isocyanate", "has_isothiocyanate",
    "has_sulfonamide", "has_sulfonate", "has_sulfonic_acid",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs
    return len(info.get("sulfonyl_chlorides") or []) == 1 and _no_fgs(info, _SC_BAD)


def _arm_n(mol, c_idx: int, forbid: int) -> int:
    from namepredict.layer2.chain_walk import _longest_from
    return len(_longest_from(mol, c_idx, {forbid}) or [c_idx])


def _arm_ok(mol, c_idx: int, forbid: int) -> bool:
    from namepredict.layer2.chain_walk import _longest_from
    from namepredict.layer2.parent_selector import _arm_ok as _ok
    arm = _longest_from(mol, c_idx, {forbid}) or [c_idx]
    return _ok(mol, arm, forbid)


def _aryl_of(mol, c_idx: int, parent: int) -> dict | None:
    atom = mol.GetAtomWithIdx(c_idx)
    if not (atom.GetIsAromatic() and atom.GetAtomicNum() == 6):
        return None
    from namepredict.layer2.aryl_sub import _phenyl_at
    ph = _phenyl_at(mol, c_idx, parent)
    return None if ph is None else {"c": c_idx, "ph": ph}


def _alkyl_s(mol, c: int, s: int) -> dict | None:
    if mol.GetAtomWithIdx(c).IsInRing() or mol.GetAtomWithIdx(c).GetIsAromatic():
        return None
    if not _arm_ok(mol, c, s):
        return None
    n = _arm_n(mol, c, s)
    return {"kind": "alkyl", "n": n} if 1 <= n <= 4 else None


def _s_side(mol, e: dict) -> dict | None:
    """S–C side: simple aryl (Ph ≤ leaves) or n-alkyl C1–C4 (+ terminal F/CF3)."""
    c, s = e["c_attach"], e["s_idx"]
    ar = _aryl_of(mol, c, s)
    return {"kind": "aryl", **ar} if ar is not None else _alkyl_s(mol, c, s)


def _s_chain(mol, e: dict, s: dict) -> list[int]:
    from namepredict.layer2.chain_walk import _longest_from
    if s["kind"] != "alkyl":
        return [e["c_attach"]]
    return _longest_from(mol, e["c_attach"], {e["s_idx"]}) or [e["c_attach"]]


def _pack(info: dict, e: dict, s: dict) -> dict:
    from namepredict.layer2.parent_selector import _parent_dict
    mol = info["mol"]
    return _parent_dict(
        _s_chain(mol, e, s), "sulfonyl_chloride",
        s_idx=e["s_idx"], c_attach=e["c_attach"], cl_idx=e["cl_idx"],
        s_side=s, mode=s["kind"], mol=mol,
    )


def _sulfonyl_chloride_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e, mol = info["sulfonyl_chlorides"][0], info["mol"]
    s = _s_side(mol, e)
    return None if s is None else _pack(info, e, s)
