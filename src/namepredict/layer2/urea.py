"""L2 simple mono-urea functional parent (P-66.1.6.1.1)."""
from __future__ import annotations


_UREA_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs
    return len(info.get("ureas") or []) == 1 and _no_fgs(info, _UREA_BAD)


def _arm_n(mol, c_idx: int, n_idx: int) -> int:
    from namepredict.layer2.parent_core import _longest_from
    return len(_longest_from(mol, c_idx, {n_idx}) or [c_idx])


def _aryl_of(mol, n_idx: int, c_idx: int) -> dict | None:
    atom = mol.GetAtomWithIdx(c_idx)
    if not (atom.GetIsAromatic() and atom.GetAtomicNum() == 6):
        return None
    from namepredict.layer2.aryl_sub import _phenyl_at
    from namepredict.layer2.aryl_stem import with_aryl_names
    ph = _phenyl_at(mol, c_idx, n_idx)
    if ph is None:
        return None
    base = {"aryl_c": c_idx, "n_idx": n_idx, "ph": ph}
    return with_aryl_names(base, mol, c_idx, n_idx)


def _alkyl_ns(mol, n_idx: int, cs: list[int]) -> list[int] | None:
    if not cs or any(mol.GetAtomWithIdx(c).GetIsAromatic() for c in cs):
        return None
    return [_arm_n(mol, c, n_idx) for c in cs]


def _side_mono(mol, n_idx: int, c: int) -> dict:
    ar = _aryl_of(mol, n_idx, c)
    if ar is not None:
        return {"kind": "aryl", **ar}
    ns = _alkyl_ns(mol, n_idx, [c])
    return {"kind": "alkyl", "n_idx": n_idx, "ns": ns or []}


def _side_multi(mol, n_idx: int, cs: list[int]) -> dict:
    ns = _alkyl_ns(mol, n_idx, cs)
    if ns is None:
        return {"kind": "complex", "n_idx": n_idx}
    return {"kind": "dialkyl", "n_idx": n_idx, "ns": ns}


def _side(mol, n_idx: int, cs: list[int]) -> dict:
    """Classify one urea N side: free NH2 / mono-aryl / dialkyl."""
    if not cs:
        return {"kind": "h", "n_idx": n_idx}
    return _side_mono(mol, n_idx, cs[0]) if len(cs) == 1 else _side_multi(mol, n_idx, cs)


def _is_me2(s: dict) -> bool:
    return s["kind"] == "dialkyl" and list(s.get("ns") or []) == [1, 1]


def _simple_sides(s1: dict, s2: dict) -> bool:
    """Accept: H+H, aryl+H, Me2+aryl. Reject complex / di-aryl."""
    kinds = {s1["kind"], s2["kind"]}
    if kinds in ({"h"}, {"aryl", "h"}):
        return True
    return kinds == {"dialkyl", "aryl"} and (_is_me2(s1) or _is_me2(s2))


def _order_sides(s1: dict, s2: dict) -> tuple[dict, dict]:
    """N1 = dialkyl/H preferred as locant-1; N3 = aryl or other."""
    if s1["kind"] == "dialkyl" or (s1["kind"] == "h" and s2["kind"] == "aryl"):
        return s1, s2
    if s2["kind"] == "dialkyl" or (s2["kind"] == "h" and s1["kind"] == "aryl"):
        return s2, s1
    return s1, s2


def _sides_ok(info: dict, e: dict) -> tuple[dict, dict] | None:
    mol = info["mol"]
    s1 = _side(mol, e["n1_idx"], list(e.get("n1_c_idxs") or []))
    s2 = _side(mol, e["n2_idx"], list(e.get("n2_c_idxs") or []))
    return _order_sides(s1, s2) if _simple_sides(s1, s2) else None


def _pack_parent(info: dict, e: dict, n1: dict, n3: dict) -> dict:
    from namepredict.layer2.parent_core import _parent_dict
    return _parent_dict(
        [e["c_idx"]], "urea",
        c_idx=e["c_idx"], o_idx=e["o_idx"], mol=info["mol"],
        enol=e.get("enol"), n1=n1, n3=n3,
    )


def _urea_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e = info["ureas"][0]
    ordered = _sides_ok(info, e)
    return None if ordered is None else _pack_parent(info, e, *ordered)
