"""L2 simple mono-sulfonamide functional parent (P-65.3)."""
from __future__ import annotations


_SA_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
    "has_sulfoxide", "has_isocyanate", "has_isothiocyanate", "has_sulfonate",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs
    return len(info.get("sulfonamides") or []) == 1 and _no_fgs(info, _SA_BAD)


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


def _cyclo_of(mol, c_idx: int, parent: int) -> dict | None:
    from namepredict.layer2.side_cycloalkyl import _is_monocycloalkyl
    ring = _is_monocycloalkyl(mol, c_idx, {parent})
    if ring is None or not (3 <= len(ring) <= 6):
        return None
    return {"c": c_idx, "size": len(ring)}


def _s_side(mol, e: dict) -> dict | None:
    c, s = e["c_attach"], e["s_idx"]
    ar = _aryl_of(mol, c, s)
    if ar is not None:
        return {"kind": "aryl", **ar}
    if mol.GetAtomWithIdx(c).IsInRing() or mol.GetAtomWithIdx(c).GetIsAromatic():
        return None
    if not _arm_ok(mol, c, s):
        return None
    n = _arm_n(mol, c, s)
    return {"kind": "alkyl", "n": n} if 1 <= n <= 4 else None


def _n_c_side(mol, n_idx: int, c: int) -> dict | None:
    ar = _aryl_of(mol, c, n_idx)
    if ar is not None:
        return {"kind": "aryl", **ar}
    cy = _cyclo_of(mol, c, n_idx)
    return {"kind": "cyclo", **cy} if cy is not None else None


def _n_side(mol, e: dict) -> dict | None:
    n_idx, cs = e["n_idx"], list(e.get("n_c_idxs") or [])
    if not cs:
        return {"kind": "h"}
    return _n_c_side(mol, n_idx, cs[0]) if len(cs) == 1 else None


_MODE_MAP = {
    ("alkyl", "h"): "alkyl", ("aryl", "h"): "aryl",
    ("alkyl", "aryl"): "n_aryl_alkyl", ("aryl", "cyclo"): "n_cyclo_aryl",
}


def _mode(s: dict, n: dict) -> str | None:
    return _MODE_MAP.get((s["kind"], n["kind"]))


def _s_chain(mol, e: dict, s: dict) -> list[int]:
    """Alkyl: full open arm; aryl: only S-attach (L5 builds ring stem)."""
    from namepredict.layer2.chain_walk import _longest_from
    if s["kind"] != "alkyl":
        return [e["c_attach"]]
    return _longest_from(mol, e["c_attach"], {e["s_idx"]}) or [e["c_attach"]]


def _pack(info: dict, e: dict, s: dict, n: dict, mode: str) -> dict:
    from namepredict.layer2.parent_selector import _parent_dict
    mol = info["mol"]
    return _parent_dict(
        _s_chain(mol, e, s), "sulfonamide",
        s_idx=e["s_idx"], c_attach=e["c_attach"], n_idx=e["n_idx"],
        s_side=s, n_side=n, mode=mode, mol=mol,
    )


def _sulfonamide_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e, mol = info["sulfonamides"][0], info["mol"]
    s, n = _s_side(mol, e), _n_side(mol, e)
    if s is None or n is None:
        return None
    mode = _mode(s, n)
    return None if mode is None else _pack(info, e, s, n, mode)
