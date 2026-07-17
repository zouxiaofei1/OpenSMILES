"""Simple mono carbamate parent (P-65): alkyl N-carbamate."""
from __future__ import annotations


_CB_BAD = (
    "has_acid", "has_ester", "has_amide", "has_nitrile", "has_acyl_chloride",
    "has_aldehyde", "has_ketone", "has_anhydride", "has_thiol",
    "has_phosphate", "has_phosphonic",
)


def _mono_ok(info: dict) -> bool:
    from namepredict.layer2.parent_core import _no_fgs
    return len(info.get("carbamates") or []) == 1 and _no_fgs(info, _CB_BAD)


def _alkoxy_label(info: dict, e: dict) -> tuple[str, str, int | None]:
    from namepredict.layer2.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(info["mol"], e["o_idx"], e["alkoxy_c_idx"])
    return side["alkoxy_en"], side["alkoxy_zh"], side["alkoxy_n"]


def _arm_n(mol, c_idx: int, n_idx: int) -> int:
    from namepredict.layer2.parent_core import _longest_from
    return len(_longest_from(mol, c_idx, {n_idx}) or [c_idx])


def _n_aryl_fields(mol, n_idx: int, c_idx: int) -> dict:
    from namepredict.layer2.aryl_stem import aryl_en_zh
    en, zh = aryl_en_zh(mol, c_idx, n_idx)
    return {"n_phenyl": True, "n_aryl_c": c_idx, "n_aryl_en": en, "n_aryl_zh": zh}


def _one_n_sub(mol, n_idx: int, c_idx: int) -> dict:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetIsAromatic() and atom.GetAtomicNum() == 6:
        return _n_aryl_fields(mol, n_idx, c_idx)
    from namepredict.layer2.aryl_sub import _ch2_ph_at
    if _ch2_ph_at(mol, c_idx, n_idx) is not None:
        return {"n_benzyl": True, "n_benzyl_ch2": c_idx}
    return {"n_alkyl_n": _arm_n(mol, c_idx, n_idx)}


def _n_arm_label(mol, n_idx: int, c_idx: int) -> dict:
    """One N-sub label: kind + en/zh (+n for alkyl)."""
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetIsAromatic() and atom.GetAtomicNum() == 6:
        from namepredict.layer2.aryl_stem import aryl_en_zh
        en, zh = aryl_en_zh(mol, c_idx, n_idx)
        return {"kind": "aryl", "en": en, "zh": zh, "c": c_idx}
    n = _arm_n(mol, c_idx, n_idx)
    return {"kind": "alk", "n": n, "en": None, "zh": None, "c": c_idx}


_N_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _one_alk_fill(lab: dict) -> dict:
    if lab["kind"] != "alk":
        return lab
    n = int(lab.get("n") or 0)
    return {**lab, "en": _N_EN.get(n, "alkyl"), "zh": _N_ZH.get(n, "烷基")}


def _fill_alk_labels(labs: list[dict]) -> list[dict]:
    return [_one_alk_fill(lab) for lab in labs]


def _n_meta(info: dict, e: dict) -> dict:
    mol, n_cs = info["mol"], list(e.get("n_c_idxs") or [])
    out: dict = {"n_c_idxs": n_cs, "n_idx": e["n_idx"]}
    if len(n_cs) == 1:
        out.update(_one_n_sub(mol, e["n_idx"], n_cs[0]))
    elif len(n_cs) == 2:
        out["n_alkyl_ns"] = [_arm_n(mol, c, e["n_idx"]) for c in n_cs]
        labs = [_n_arm_label(mol, e["n_idx"], c) for c in n_cs]
        out["n_labels"] = _fill_alk_labels(labs)
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
    from namepredict.layer2.parent_core import _longest_from
    return _longest_from(info["mol"], e["alkoxy_c_idx"]) or [e["c_idx"]]


def _carbamate_parent(info: dict) -> dict | None:
    if not _mono_ok(info):
        return None
    e = info["carbamates"][0]
    meta = _parent_meta(info, e)
    from namepredict.layer2.parent_core import _parent_dict
    return _parent_dict(_parent_chain(info, e, meta.get("alkoxy_n")), "carbamate", **meta)
