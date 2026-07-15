"""L2 organic carbonate parent (P-65.6): dialkyl or alkyl–aryl carbonate."""
from __future__ import annotations


_CB_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
    "has_sulfoxide", "has_isocyanate", "has_isothiocyanate", "has_sulfonamide",
    "has_sulfonate",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs
    return len(info.get("carbonates") or []) == 1 and _no_fgs(info, _CB_BAD)


def _arm_ok(mol, c_idx: int, o_idx: int) -> bool:
    from namepredict.layer2.chain_walk import _longest_from
    from namepredict.layer2.parent_selector import _arm_ok as _ok
    arm = _longest_from(mol, c_idx, {o_idx}) or [c_idx]
    return _ok(mol, arm, o_idx)


def _alkyl_side(mol, o_idx: int, c_idx: int) -> dict | None:
    """n-alkyl C1–C4 only (first cut)."""
    from namepredict.layer2.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(mol, o_idx, c_idx)
    n = side.get("alkoxy_n")
    if n is None or not (1 <= int(n) <= 4):
        return None
    if not _arm_ok(mol, c_idx, o_idx):
        return None
    return {"kind": "alkyl", "n": int(n), "o_idx": o_idx, "c_idx": c_idx}


def _aryl_side(mol, o_idx: int, c_idx: int) -> dict | None:
    """Unfused Ph with ≤3 Me/halo leaves (aryl_sub gate)."""
    atom = mol.GetAtomWithIdx(c_idx)
    if not (atom.GetIsAromatic() and atom.GetAtomicNum() == 6):
        return None
    from namepredict.layer2.aryl_sub import _phenyl_at
    ph = _phenyl_at(mol, c_idx, o_idx)
    if ph is None:
        return None
    return {"kind": "aryl", "o_idx": o_idx, "c_idx": c_idx, "ph": ph}


def _side(mol, o_idx: int, c_idx: int) -> dict | None:
    return _aryl_side(mol, o_idx, c_idx) or _alkyl_side(mol, o_idx, c_idx)


def _sides(mol, e: dict) -> tuple[dict, dict] | None:
    s1 = _side(mol, e["o1_idx"], e["alkoxy1_c_idx"])
    s2 = _side(mol, e["o2_idx"], e["alkoxy2_c_idx"])
    if s1 is None or s2 is None:
        return None
    return s1, s2


def _mode_of(s1: dict, s2: dict) -> str | None:
    kinds = {s1["kind"], s2["kind"]}
    if kinds == {"alkyl"} and s1["n"] == s2["n"]:
        return "sym_dialkyl"
    if kinds == {"alkyl", "aryl"}:
        return "alkyl_aryl"
    return None


def _chain_of(mol, e: dict, s1: dict, s2: dict) -> list[int]:
    from namepredict.layer2.chain_walk import _longest_from
    alk = s1 if s1["kind"] == "alkyl" else s2
    return _longest_from(mol, alk["c_idx"], {alk["o_idx"]}) or [e["c_idx"]]


def _pack(info: dict, e: dict, s1: dict, s2: dict, mode: str) -> dict:
    from namepredict.layer2.parent_selector import _parent_dict
    mol = info["mol"]
    return _parent_dict(
        _chain_of(mol, e, s1, s2), "carbonate",
        c_idx=e["c_idx"], side1=s1, side2=s2, mode=mode, mol=mol,
    )


def _carbonate_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e, mol = info["carbonates"][0], info["mol"]
    got = _sides(mol, e)
    if got is None:
        return None
    s1, s2 = got
    mode = _mode_of(s1, s2)
    return None if mode is None else _pack(info, e, s1, s2, mode)
