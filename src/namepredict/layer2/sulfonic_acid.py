"""L2 simple mono free sulfonic acid / sulfonate parent (P-65.3)."""
from __future__ import annotations


_SA_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
    "has_sulfoxide", "has_isocyanate", "has_isothiocyanate",
    "has_sulfonamide", "has_sulfonate", "has_sulfonyl_chloride",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs
    return len(info.get("sulfonic_acids") or []) == 1 and _no_fgs(info, _SA_BAD)


def _arm_n(mol, c_idx: int, forbid: int) -> int:
    from namepredict.layer2.chain_walk import _longest_from
    return len(_longest_from(mol, c_idx, {forbid}) or [c_idx])


def _arm_ok(mol, c_idx: int, forbid: int) -> bool:
    from namepredict.layer2.chain_walk import _longest_from
    from namepredict.layer2.parent_core import _arm_ok as _ok
    arm = _longest_from(mol, c_idx, {forbid}) or [c_idx]
    return _ok(mol, arm, forbid)


def _is_me_or_halo(mol, nb, ring_i: int) -> bool:
    from namepredict.layer2.aryl_sub import _is_terminal_halo, _is_terminal_me_leaf
    if _is_terminal_halo(nb):
        return True
    return nb.GetAtomicNum() == 6 and _is_terminal_me_leaf(mol, nb.GetIdx(), ring_i)


def _one_leaf(mol, nb, i: int, attach: int, parent: int) -> int | None:
    if i == attach and nb.GetIdx() == parent:
        return 0
    return 1 if _is_me_or_halo(mol, nb, i) else None


def _leaf_n(mol, ph, attach: int, parent: int) -> int | None:
    from namepredict.layer2.aryl_sub import _nb_outside
    n = 0
    for i in ph:
        for nb in _nb_outside(mol, i, ph):
            k = _one_leaf(mol, nb, i, attach, parent)
            if k is None:
                return None
            n += k
    return n


def _ph_leaf_ok(mol, ph, attach: int, parent: int) -> bool:
    n = _leaf_n(mol, ph, attach, parent)
    return n is not None and n <= 2


def _aryl_of(mol, c_idx: int, parent: int) -> dict | None:
    """Unfused Ph with only terminal Me/halo leaves, count ≤2 (first cut)."""
    atom = mol.GetAtomWithIdx(c_idx)
    if not (atom.GetIsAromatic() and atom.GetAtomicNum() == 6):
        return None
    from namepredict.layer2.aryl_sub import _phenyl_at
    from namepredict.layer2.aryl_stem import with_aryl_names
    ph = _phenyl_at(mol, c_idx, parent)
    if ph is None or not _ph_leaf_ok(mol, ph, c_idx, parent):
        return None
    return with_aryl_names({"c": c_idx, "ph": ph}, mol, c_idx, parent)


def _alkyl_s(mol, c: int, s: int) -> dict | None:
    if mol.GetAtomWithIdx(c).IsInRing() or mol.GetAtomWithIdx(c).GetIsAromatic():
        return None
    if not _arm_ok(mol, c, s):
        return None
    n = _arm_n(mol, c, s)
    return {"kind": "alkyl", "n": n} if 1 <= n <= 4 else None


def _s_side(mol, e: dict) -> dict | None:
    """S–C side: simple aryl (Ph ≤ leaves) or n-alkyl C1–C4."""
    c, s = e["c_attach"], e["s_idx"]
    ar = _aryl_of(mol, c, s)
    return {"kind": "aryl", **ar} if ar is not None else _alkyl_s(mol, c, s)


def _s_chain(mol, e: dict, s: dict) -> list[int]:
    from namepredict.layer2.chain_walk import _longest_from
    if s["kind"] != "alkyl":
        return [e["c_attach"]]
    return _longest_from(mol, e["c_attach"], {e["s_idx"]}) or [e["c_attach"]]


def _pack(info: dict, e: dict, s: dict) -> dict:
    from namepredict.layer2.parent_core import _parent_dict
    mol = info["mol"]
    return _parent_dict(
        _s_chain(mol, e, s), "sulfonic_acid",
        s_idx=e["s_idx"], c_attach=e["c_attach"], o_idx=e["o_idx"],
        s_side=s, mode=s["kind"], mol=mol, anion=bool(e.get("anion")),
    )


def _sulfonic_acid_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e, mol = info["sulfonic_acids"][0], info["mol"]
    s = _s_side(mol, e)
    return None if s is None else _pack(info, e, s)
