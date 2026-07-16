"""L2 simple mono-guanidine functional parent (P-66.4.1.2.1)."""
from __future__ import annotations


_GU_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic", "has_carbamate", "has_urea",
    "has_sulfoxide", "has_isocyanate", "has_isothiocyanate", "has_sulfonate",
    "has_sulfonyl_chloride", "has_sulfonamide",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs
    return len(info.get("guanidines") or []) == 1 and _no_fgs(info, _GU_BAD)


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


def _is_so2_s(mol, s_idx: int, n_idx: int) -> bool:
    """True if S is SO2 single-bonded to N (arylsulfonyl bridge)."""
    from namepredict.layer1.sulfonamide import _dbl_o_nbs, _sgl_nbs
    s = mol.GetAtomWithIdx(s_idx)
    if s.GetAtomicNum() != 16 or s.GetTotalDegree() != 4 or len(_dbl_o_nbs(s)) != 2:
        return False
    ns, cs = _sgl_nbs(s, 7), _sgl_nbs(s, 6)
    return any(n.GetIdx() == n_idx for n in ns) and len(cs) == 1


def _pack_as(n_idx: int, s_idx: int, ar: dict) -> dict:
    out = {
        "kind": "arylsulfonyl", "n_idx": n_idx, "s_idx": s_idx,
        "aryl_c": ar["aryl_c"], "ph": ar["ph"],
    }
    if "en" in ar:
        out["en"], out["zh"] = ar["en"], ar.get("zh") or ar["en"]
    return out


def _arylsulfonyl_of(mol, n_idx: int, s_idx: int) -> dict | None:
    if not _is_so2_s(mol, s_idx, n_idx):
        return None
    from namepredict.layer1.sulfonamide import _sgl_nbs
    c = _sgl_nbs(mol.GetAtomWithIdx(s_idx), 6)[0]
    ar = _aryl_of(mol, s_idx, c.GetIdx())
    return None if ar is None else _pack_as(n_idx, s_idx, ar)


def _side_from_rest(mol, n_idx: int, rest_idxs: list[int], rest_zs: list[int]) -> dict:
    """Classify one guanidine N: free H / mono-aryl / mono-arylsulfonyl."""
    if len(rest_idxs) != 1:
        return {"kind": "h" if not rest_idxs else "complex", "n_idx": n_idx}
    z, r = rest_zs[0], rest_idxs[0]
    if z == 16:
        return _arylsulfonyl_of(mol, n_idx, r) or {"kind": "complex", "n_idx": n_idx}
    ar = _aryl_of(mol, n_idx, r)
    return {"kind": "aryl", **ar} if ar is not None else {"kind": "complex", "n_idx": n_idx}


def _sides(mol, e: dict) -> list[dict]:
    out = []
    for s in e.get("n_sides") or []:
        out.append(_side_from_rest(
            mol, s["n_idx"], list(s.get("rest_idxs") or []),
            list(s.get("rest_zs") or []),
        ))
    return out


def _simple_sides(sides: list[dict]) -> bool:
    """Accept: all H; one aryl + two H; one arylsulfonyl + two H."""
    if len(sides) != 3 or any(s["kind"] == "complex" for s in sides):
        return False
    kinds = [s["kind"] for s in sides]
    n_h, n_ar, n_as = kinds.count("h"), kinds.count("aryl"), kinds.count("arylsulfonyl")
    return (n_h == 3) or (n_h == 2 and n_ar == 1) or (n_h == 2 and n_as == 1)


def _pick_sub(sides: list[dict]) -> dict | None:
    """The single non-H side, or None if unsubstituted."""
    subs = [s for s in sides if s["kind"] != "h"]
    return subs[0] if len(subs) == 1 else None


def _pack_parent(info: dict, e: dict, sides: list[dict]) -> dict:
    from namepredict.layer2.parent_selector import _parent_dict
    sub = _pick_sub(sides)
    return _parent_dict(
        [e["c_idx"]], "guanidine",
        c_idx=e["c_idx"], mol=info["mol"],
        sides=sides, sub=sub,
        mode="unsub" if sub is None else sub["kind"],
    )


def _guanidine_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e, mol = info["guanidines"][0], info["mol"]
    sides = _sides(mol, e)
    return None if not _simple_sides(sides) else _pack_parent(info, e, sides)
