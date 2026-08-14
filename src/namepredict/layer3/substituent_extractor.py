from __future__ import annotations

from rdkit.Chem import Mol

import namepredict.layer3.side_facts as side_facts
from namepredict.cache.common_names import CommonNameCache
from namepredict.layer3.amino_side import (
    _extract_aminos as _extract_aminos_impl,
    _principal_attachments,
)

def _strip_ital_prefix(stem: str) -> str:
    if stem.startswith("tert-") or stem.startswith("sec-"):
        return stem[stem.index("-") + 1 :]
    return stem

def _strip_n_prefix(stem: str) -> str:
    if stem.startswith("N,"):
        return stem.split("-")[-1] if "-" in stem else stem
    return stem[2:] if stem.startswith("N-") else stem

def _strip_lead_locant(stem: str) -> str:
    """去掉一个前导位次集：'4-'、'1,3-'、'1,1,1-'、'1H-'（P-14.5）。"""
    i = 0
    n = len(stem)
    while i < n and stem[i].isdigit():
        i += 1
        while i < n and stem[i] == ",":
            i += 1
            while i < n and stem[i].isdigit():
                i += 1
    # 指示氢前缀：1H-、2H-、3H-（P-14.5 / P-65.3.2.5）
    if i and i + 1 < n and stem[i] == "H" and stem[i + 1] == "-":
        i += 1
    return stem[i + 1 :] if i and i < n and stem[i] == "-" else stem

def _strip_outer_parens(stem: str) -> str:
    if len(stem) >= 2 and stem[0] == "(" and stem[-1] == ")":
        return stem[1:-1]
    return stem

def alkyl_alpha_key(stem: str) -> str:
    """字母数字序键：忽略 sec-/tert-/N-/括号/前导位次（P-14.5）。"""
    s = _strip_n_prefix(_strip_ital_prefix(stem))
    s = _strip_outer_parens(s)
    return _strip_lead_locant(s)

HALO_EN = {9: "fluoro", 17: "chloro", 35: "bromo", 53: "iodo"}
HALO_ZH = {9: "氟", 17: "氯", 35: "溴", 53: "碘"}

def _side_starts(mol: Mol, chain: list[int]) -> list[tuple[int, int]]:
    cs = set(chain)
    return [(c, n) for c in chain for n in side_facts.carbon_neighbors(mol, c) if n not in cs]








def _one_anchored_alkyl(mol: Mol, attach: int, start: int, chain_set: set[int], *, name_mode: str) -> dict | None:
    """尝试用锚定 canonical-SMILES 表解析从 `start` 开始的侧链。

    使用 side_atoms（完整非母体连通组分）作为原子集合，选取取代基一侧的
    连接原子，再在 name_mode 下查找锚定键。这里只 claim "alkyl" kind 条目
    （纯碳侧链）；杂原子叶子（cyano/nitroso/...）仍由 FG/claim 提取器处理。
    未命中则继续向下。
    """
    from namepredict.tools.anchored_table import anchored_entry
    from namepredict.tools.block_cut import side_atoms

    atoms = side_atoms(mol, frozenset(chain_set), attach, frozenset({start}))
    if not atoms:
        return None
    entry = anchored_entry(mol, atoms, name_mode=name_mode)
    if entry is None or entry[3] != "alkyl":
        return None
    en, zh, paren, _kind = entry
    return {
        "kind": "alkyl", "n_carbons": len(atoms), "attach_idx": attach,
        "atoms": sorted(atoms), "en": en, "zh": zh, "paren": paren,
    }


def _one_alkyl(mol: Mol, attach: int, start: int, chain_set: set[int], *, name_mode: str = "general") -> dict | None:
    anchored = _one_anchored_alkyl(mol, attach, start, chain_set, name_mode=name_mode)
    return anchored


def _make_halo(attach: int, halo_idx: int, z: int) -> dict:
    return {
        "kind": "halo", "attach_idx": attach, "atoms": [halo_idx],
        "en": HALO_EN[z], "zh": HALO_ZH[z],
    }

def _halo_on_carbon(mol: Mol, c_idx: int) -> list[dict]:
    return [
        _make_halo(c_idx, n.GetIdx(), n.GetAtomicNum())
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() in HALO_EN
    ]

def _extract_halos(mol: Mol, chain: list[int]) -> list[dict]:
    return [h for c in chain for h in _halo_on_carbon(mol, c)]

def _filter_fg_halos(halos: list, parent: dict) -> list:
    # 官能团类醚臂已编码 F（如 HFIP）；不要重复加前缀。
    if parent.get("kind") == "ether" and parent.get("ether_arms"):
        return []
    if parent.get("kind") not in ("acyl_chloride", "acyl_bromide"):
        return halos
    cl = parent.get("cl_idx") or parent.get("hal_idx")
    return [h for h in halos if cl not in (h.get("atoms") or [])]

_PARENT_OH_KINDS = frozenset({"alcohol", "phenol", "benzenediol"})
_PARENT_NH2_KINDS = frozenset({"amine", "sec_amine", "tert_amine", "aniline"})
_PARENT_OXO_KINDS = frozenset({"ketone", "dione"})

def _make_hydroxy(attach: int, o_idx: int) -> dict:
    return {
        "kind": "hydroxy", "attach_idx": attach, "atoms": [o_idx],
        "en": "hydroxy", "zh": "羟基",
    }

def _make_oxo(attach: int) -> dict:
    return {
        "kind": "oxo", "attach_idx": attach, "atoms": [attach],
        "en": "oxo", "zh": "氧代",
    }

def _extract_hydroxys(info: dict, parent: dict) -> list[dict]:
    principal = _principal_attachments(parent, "alcohol")
    if parent.get("kind") in _PARENT_OH_KINDS and not principal:
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_hydroxy(h["c_idx"], h["o_idx"])
        for h in info.get("hydroxyls") or []
        if h["c_idx"] in chain and h["c_idx"] not in principal
    ]

def _extract_aminos(info: dict, parent: dict) -> list[dict]:
    return _extract_aminos_impl(info, parent, _PARENT_NH2_KINDS)

def _extract_oxos(info: dict, parent: dict) -> list[dict]:
    if parent.get("kind") in _PARENT_OXO_KINDS:
        return []
    chain = set(parent.get("chain") or [])
    return [
        _make_oxo(k["c_idx"]) for k in info.get("ketones") or [] if k["c_idx"] in chain
    ]

def _extract_alkyls_no_aryl(mol: Mol, chain: list[int], *, name_mode: str = "general") -> list[dict]:
    cs, out = set(chain), []
    for attach, start in _side_starts(mol, chain):
        one = _one_alkyl(mol, attach, start, cs, name_mode=name_mode)
        if one is not None:
            out.append(one)
    return out

def _extract_core_subs(info: dict, parent: dict) -> list:
    mol, chain = info["mol"], parent.get("chain") or []
    halo = _filter_fg_halos(_extract_halos(mol, chain), parent)
    return (
        halo + _extract_hydroxys(info, parent) + _extract_aminos(info, parent)
        + _extract_oxos(info, parent)
    )

def _with_full_atoms(mol, owned, s: dict) -> dict:
    """用取代基的完整非母体连通组分（side_atoms）替换其原子集，
    使叶子（halo/alkoxy/嵌套环）只计数一次。"""
    from namepredict.tools.block_cut import side_atoms

    seed = frozenset(s.get("atoms") or [])
    attach = s.get("attach_idx")
    if not seed or attach is None or attach not in owned:
        return s
    full = side_atoms(mol, owned, attach, seed)
    return s if not full else {**s, "atoms": sorted(full)}


def extract_substituents(info: dict, parent: dict, *, name_mode: str = "general", cache: CommonNameCache | None = None) -> list:
    from namepredict.layer3.claim_extract import extract_claimed_sides

    mol, chain = info["mol"], parent.get("chain") or []
    base = (
       _extract_core_subs(info, parent)
       +  _extract_alkyls_no_aryl(mol, chain, name_mode=name_mode)
    )
    owned = parent.get("owned_atoms")
    if owned:
        base = [_with_full_atoms(mol, owned, s) for s in base]
    return base + extract_claimed_sides(info, parent, base, name_mode=name_mode, cache=cache)
