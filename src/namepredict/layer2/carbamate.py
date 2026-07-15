"""Simple mono carbamate parent (P-65): alkyl N-carbamate."""
from __future__ import annotations


_CB_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_selector import _no_fgs
    return len(info.get("carbamates") or []) == 1 and _no_fgs(info, _CB_BAD)


def _methyl_count(c, o_idx: int, alkoxy_c: int) -> int:
    n = 0
    for nb in c.GetNeighbors():
        if nb.GetIdx() == o_idx or nb.GetAtomicNum() != 6:
            continue
        heavy = [x for x in nb.GetNeighbors() if x.GetAtomicNum() != 1]
        if len(heavy) == 1 and heavy[0].GetIdx() == alkoxy_c:
            n += 1
    return n


def _is_tert_butyl_o(mol, o_idx: int, alkoxy_c: int) -> bool:
    c = mol.GetAtomWithIdx(alkoxy_c)
    nbs = [n for n in c.GetNeighbors() if n.GetAtomicNum() != 1]
    return c.GetAtomicNum() == 6 and len(nbs) == 4 and _methyl_count(c, o_idx, alkoxy_c) == 3


def _alkoxy_label(info: dict, e: dict) -> tuple[str, str, int | None]:
    mol, o_idx, ac = info["mol"], e["o_idx"], e["alkoxy_c_idx"]
    if _is_tert_butyl_o(mol, o_idx, ac):
        return "tert-butyl", "叔丁", None
    from namepredict.layer2.parent_selector import _longest_from
    chain = _longest_from(mol, ac)
    return "", "", len(chain) if chain else 1


def _arm_n(mol, c_idx: int, n_idx: int) -> int:
    from namepredict.layer2.parent_selector import _longest_from
    return len(_longest_from(mol, c_idx, {n_idx}) or [c_idx])


def _one_n_sub(mol, n_idx: int, c_idx: int) -> dict:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetIsAromatic() and atom.GetAtomicNum() == 6:
        return {"n_phenyl": True, "n_aryl_c": c_idx}
    from namepredict.layer2.aryl_sub import _ch2_ph_at
    if _ch2_ph_at(mol, c_idx, n_idx) is not None:
        return {"n_benzyl": True, "n_benzyl_ch2": c_idx}
    return {"n_alkyl_n": _arm_n(mol, c_idx, n_idx)}


def _n_meta(info: dict, e: dict) -> dict:
    mol, n_cs = info["mol"], list(e.get("n_c_idxs") or [])
    out: dict = {"n_c_idxs": n_cs, "n_idx": e["n_idx"]}
    if len(n_cs) == 1:
        out.update(_one_n_sub(mol, e["n_idx"], n_cs[0]))
    elif len(n_cs) == 2:
        out["n_alkyl_ns"] = [_arm_n(mol, c, e["n_idx"]) for c in n_cs]
    return out


def _parent_meta(info: dict, e: dict) -> dict:
    en_a, zh_a, n = _alkoxy_label(info, e)
    return {
        "c_idx": e["c_idx"], "o_idx": e["o_idx"], "alkoxy_c_idx": e["alkoxy_c_idx"],
        "alkoxy_en": en_a, "alkoxy_zh": zh_a, "alkoxy_n": n, "mol": info["mol"],
        **_n_meta(info, e),
    }


def _parent_chain(info: dict, e: dict, n: int | None) -> list[int]:
    if not n:
        return [e["c_idx"]]
    from namepredict.layer2.parent_selector import _longest_from
    return _longest_from(info["mol"], e["alkoxy_c_idx"]) or [e["c_idx"]]


def _carbamate_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e = info["carbamates"][0]
    meta = _parent_meta(info, e)
    from namepredict.layer2.parent_selector import _parent_dict
    return _parent_dict(_parent_chain(info, e, meta.get("alkoxy_n")), "carbamate", **meta)
