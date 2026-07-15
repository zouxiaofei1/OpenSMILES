"""L5 assembly for simple carbamate (P-65)."""
from __future__ import annotations

from namepredict.layer5.stems import ESTER_ALKYL_EN, ESTER_ALKYL_ZH


_N_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _alkyl_pair(n: int | None) -> tuple[str, str] | None:
    if n is None:
        return None
    en, zh = ESTER_ALKYL_EN.get(n), ESTER_ALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _alkoxy_names(parent: dict) -> tuple[str, str] | None:
    if parent.get("alkoxy_en"):
        return parent["alkoxy_en"], parent.get("alkoxy_zh") or parent["alkoxy_en"]
    return _alkyl_pair(parent.get("alkoxy_n"))


def _aryl_en_zh(mol, n_idx: int, c_idx: int) -> tuple[str, str]:
    from namepredict.layer2.aryl_sub import _phenyl_at, _phenyl_name
    ph = _phenyl_at(mol, c_idx, n_idx)
    if ph is None:
        return "phenyl", "苯基"
    en, zh, _ = _phenyl_name(mol, ph, c_idx)
    return en, zh


def _arm_len(mol, c_idx: int, n_idx: int | None) -> int:
    from namepredict.layer2.parent_selector import _longest_from
    forbid = {n_idx} if n_idx is not None else set()
    return len(_longest_from(mol, c_idx, forbid) or [c_idx])


def _n_aryl_piece(parent: dict) -> tuple[str, str]:
    mol, c, n_idx = parent.get("mol"), parent.get("n_aryl_c"), parent.get("n_idx")
    if mol is not None and c is not None and n_idx is not None:
        en, zh = _aryl_en_zh(mol, n_idx, c)
    else:
        en, zh = "phenyl", "苯基"
    return (f"({en})", f"({zh})") if en != "phenyl" else ("N-phenyl", "N-苯基")


def _n_alkyl_piece(parent: dict) -> tuple[str, str]:
    n = parent.get("n_alkyl_n")
    en, zh = _N_ALKYL_EN.get(n), _N_ALKYL_ZH.get(n)
    return (f"N-{en}", f"N-{zh}") if en else ("", "")


def _n_piece(parent: dict) -> tuple[str, str]:
    mol, n_cs = parent.get("mol"), list(parent.get("n_c_idxs") or [])
    if len(n_cs) == 2 and mol is not None:
        return _mixed_n(parent, mol, n_cs)
    if parent.get("n_phenyl") or parent.get("n_aryl_c") is not None:
        return _n_aryl_piece(parent)
    if parent.get("n_benzyl"):
        return "N-benzyl", "N-苄基"
    return _n_alkyl_piece(parent)


def _label_one(mol, c, n_idx) -> tuple[str, str, str]:
    if mol.GetAtomWithIdx(c).GetIsAromatic():
        en, zh = _aryl_en_zh(mol, n_idx, c)
        return en, zh, "aryl"
    n = _arm_len(mol, c, n_idx)
    return _N_ALKYL_EN.get(n, "alkyl"), _N_ALKYL_ZH.get(n, "烷基"), "alk"


def _mixed_n(parent: dict, mol, n_cs: list[int]) -> tuple[str, str]:
    labels = [_label_one(mol, c, parent.get("n_idx")) for c in n_cs]
    labels.sort(key=lambda t: 0 if t[2] == "alk" else 1)
    a, b = labels[0], labels[1]
    return f"{a[0]}({b[0]})", f"{a[1]}({b[1]})"


def _en_join(a_en: str, n_en: str) -> str:
    if n_en.startswith("N-") or n_en:
        return f"{a_en} {n_en}carbamate"
    return f"{a_en} carbamate"


def _zh_join(a_zh: str, n_zh: str) -> str:
    if n_zh.startswith("N-") or n_zh.startswith("(") or n_zh:
        return f"{n_zh}氨基甲酸{a_zh}酯"
    return f"氨基甲酸{a_zh}酯"


def carbamate_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    alk = _alkoxy_names(parent)
    if not alk:
        return None
    a_en, a_zh = alk
    n_en, n_zh = _n_piece(parent)
    return _en_join(a_en, n_en), _zh_join(a_zh, n_zh)
