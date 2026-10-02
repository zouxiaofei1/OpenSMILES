"""L5 立体前缀：E/Z（P-91.2）与 CIP R/S（P-92）。"""
from __future__ import annotations

import re

from rdkit.Chem import BondStereo, BondType, Mol

from namepredict.layer4.locant_calc import _atom_locant, locant_key
from namepredict.layer4.numbering_engine import assign_cip


def _split_stereo_lead(name: str) -> tuple[str, str]:
    """从名称/词干中切分前导 '(…)-' 立体块。"""
    if not name.startswith("("):
        return "", name
    close = name.find(")-")
    if close < 0:
        return "", name
    return name[: close + 2], name[close + 2 :]


# --- E/Z 立体描述符 ------------------------

_STEREO_TAG = {BondStereo.STEREOE: "(E)-", BondStereo.STEREOZ: "(Z)-"}  # RDKit 立体键枚举 → E/Z 前缀标记


def _bond_stereo(mol: Mol | None, double_bond) -> str:
    """查指定双键的立体标签（无键/无立体返回空串）。"""
    if mol is None or not double_bond:
        return ""
    c1, c2 = double_bond
    bond = mol.GetBondBetweenAtoms(int(c1), int(c2))
    return _STEREO_TAG.get(bond.GetStereo(), "") if bond is not None else ""


def _ez_prefix(numbered: dict) -> str:
    """单双键母体的 E/Z 前缀（带位次，如 '(2E)-'）。"""
    parent = numbered.get("parent") or {}
    mol = parent.get("mol")
    bond = parent.get("double_bond")
    tag = _bond_stereo(mol, bond)
    if not tag:
        return ""
    letter = tag[1]  # 'E' or 'Z'
    chain = parent.get("chain") or []
    loc = _bond_min_loc(chain, bond)
    if loc is None:
        return tag
    return f"({loc}{letter})-"


def _bond_min_loc(chain: list[int], pair) -> int | None:
    """求双键两端在母体链中的较小位次（不在链上返回 None）。"""
    if not pair or pair[0] not in chain or pair[1] not in chain:
        return None
    return min(chain.index(pair[0]) + 1, chain.index(pair[1]) + 1)


def _ez_bond_part(mol, chain: list[int], bond) -> tuple[int, str] | None:
    """单条立体双键的 (位次, 字母) 部件。"""
    loc, tag = _bond_min_loc(chain, bond), _bond_stereo(mol, bond)
    return (loc, tag[1]) if loc is not None and tag else None


def _ez_parts(mol, chain: list[int], bonds) -> list[tuple[int, str]]:
    """只收集定义了立体化学的键的 (loc, letter)（部分立体也可）。"""
    parts = [_ez_bond_part(mol, chain, b) for b in bonds]
    return sorted((p for p in parts if p is not None), key=lambda x: x[0])


def _ez_multi_prefix(numbered: dict) -> str:
    """由有立体化学的键生成多烯前缀：(2E,6Z)- 或 (14Z)-。"""
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), parent.get("chain") or []
    bonds = list(parent.get("double_bonds") or [])
    if mol is None or not bonds or not chain:
        return ""
    parts = _ez_parts(mol, chain, bonds)
    if not parts:
        return ""
    return f"({','.join(f'{loc}{let}' for loc, let in parts)})-"


def ez_for_parent(numbered: dict) -> str:
    """母体 E/Z 前缀：有 double_bonds 时多烯，否则单键。"""
    parent = numbered.get("parent") or {}
    if parent.get("double_bonds"):
        return _ez_multi_prefix(numbered)
    return _ez_prefix(numbered)


def _parent_locants(parent: dict) -> dict[int, int | str]:
    """母体原子索引 → 位次（整体编号标签优先，否则链序号）。"""
    chain = parent.get("chain") or []
    facts = parent.get("numbering_scaffold") or {}
    return {int(idx): _atom_locant(chain, idx, facts) for idx in chain}


def _exo_ez_parts(numbered: dict) -> list[tuple[int | str, str]]:
    """母体外挂立体双键的 (位次, 字母)：仅一端在母体内，位次取母体侧（P-91.2）。"""
    parent = numbered.get("parent") or {}
    mol, locants = parent.get("mol"), _parent_locants(parent)
    if mol is None or not locants:
        return []
    parts: list[tuple[int | str, str]] = []
    for bond in mol.GetBonds():
        tag = _STEREO_TAG.get(bond.GetStereo(), "")
        if not tag or bond.GetBondType() != BondType.DOUBLE:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in locants and b in locants:
            continue  # 两端都在母体内：已由母体名承载
        idx = a if a in locants else b
        if idx in locants:
            parts.append((locants[idx], tag[1]))
    return parts


def join_ez_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """对外 E/Z 入口：补母体外挂双键的 E/Z 部件（母体名已有的位次不重复）。"""
    extra = _exo_ez_parts(numbered)
    if not extra:
        return en, zh
    have = {loc for loc, _ in _parse_stereo(_split_stereo_lead(en)[0])}
    extra = [(loc, let) for loc, let in extra if loc not in have]
    if not extra:
        return en, zh
    return _apply_rs(en, zh, extra, numbered)


# --- CIP R/S 立体描述符 --------------------

def _cip_on_chain(mol: Mol, chain: list[int]) -> list[tuple[int, str]]:
    """返回母体链上手性中心的 (链序号, R/S) 列表。"""
    assign_cip(mol)
    out: list[tuple[int, str]] = []
    for loc, idx in enumerate(chain, 1):
        atom = mol.GetAtomWithIdx(int(idx))
        code = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else ""
        if code in ("R", "S"):
            out.append((loc, code))
    return out


def _chain_locant(parent: dict, pos: int) -> int | str:
    """链上第 pos 位（1 起）→ locant，退回链序号。"""
    chain = parent.get("chain") or []
    if not 1 <= pos <= len(chain):
        return pos
    return _atom_locant(chain, chain[pos - 1], parent.get("numbering_scaffold") or {}) or pos


def _collapsed_parent(parent: dict) -> bool:
    """当母体 C1 来自更大/环分子折叠时跳过 R/S。"""
    n = int(parent.get("n_carbons") or 0)
    if n > 1:
        return False
    mol = parent.get("mol")
    if mol is None:
        return False
    if mol.GetRingInfo().NumRings() > 0:
        return True
    return mol.GetNumHeavyAtoms() > n + 2


def _rs_parts(numbered: dict) -> list[tuple[int | str, str]]:
    """取母体手性中心的 (位次, R/S) 列表（P-92）。

    不做 kind 门控：parent 的 kind 混用两套词汇表，含氧酸母体存的是
    L1 的 oxo_kind（sulfonic/sulfonate/…），与 FG 类别枚举对不上。
    """
    parent = numbered.get("parent") or {}
    mol, chain = parent.get("mol"), parent.get("chain") or []
    if mol is None or not chain:
        return []
    # 碳/硫单原子母体省略位次号（同 as_substituent 约定）；磷碎片位次由外部骨架给出，不适用
    if len(chain) == 1 and mol.GetAtomWithIdx(int(chain[0])).GetAtomicNum() != 15:
        return [(None, code) for _, code in _cip_on_chain(mol, chain)]
    if _collapsed_parent(parent):
        return []
    return [(_chain_locant(parent, pos), code) for pos, code in _cip_on_chain(mol, chain)]


def _parse_token(tok: str) -> tuple[int | str | None, str] | None:
    """解析单个 token：'E'→(None,'E')；'8R'→(8,'R')。"""
    if tok in ("E", "Z", "R", "S"):
        return None, tok
    m = re.match(r"^(\d+)([a-z]*)([EZRS])$", tok)
    if not m:
        return None
    num, letter, code = m.groups()
    return (int(num) if not letter else f"{num}{letter}"), code


def _parse_stereo(tag: str) -> list[tuple[int | str | None, str]]:
    """将 '(E)-' / '(2E,6Z)-' 解析为部件列表。"""
    if not tag.startswith("(") or not tag.endswith(")-"):
        return []
    raw = tag[1:-2]
    out: list[tuple[int | None, str]] = []
    for tok in raw.split(","):
        p = _parse_token(tok.strip())
        if p:
            out.append(p)
    return out


def _part_key(part: tuple[int | str | None, str]) -> tuple[bool, tuple[int, str]]:
    """立体部件排序键：无位次者排前，余按 locant_key。"""
    loc = part[0]
    return (loc is not None, locant_key(loc) if loc is not None else (0, ""))


def _format_stereo(parts: list[tuple[int | str | None, str]]) -> str:
    """将立体部件列表拼成 '(…)-' 前缀（无部件返回空串）。"""
    if not parts:
        return ""
    ordered = sorted(parts, key=_part_key)
    body = ",".join((let if loc is None else f"{loc}{let}") for loc, let in ordered)
    return f"({body})-"


def _with_rs(name: str, rs: list[tuple[int | str | None, str]]) -> str:
    """将 R/S 部件并入名称已有的立体前缀（保留非 R/S 立体，按位次重排）。"""
    if not rs:
        return name
    tag, stem = _split_stereo_lead(name)
    if tag and not _parse_stereo(tag):  # 前导括号非立体块（如 (4-chlorophenyl)-）不应被当立体标签吞掉
        tag, stem = "", name
    parts = [(loc, let) for loc, let in _parse_stereo(tag) if let not in ("R", "S")] + list(rs)
    return f"{_format_stereo(parts)}{stem}"


def _ester_en_rs(en: str, rs: list[tuple[int | str | None, str]]) -> str:
    """在酯名烷基词后插入 R/S 部件。"""
    if " " not in en:
        return _with_rs(en, rs)
    alkyl, acyl = en.split(" ", 1)
    return f"{alkyl} {_with_rs(acyl, rs)}"


def _apply_rs(en: str, zh: str, parts: list, numbered: dict) -> tuple[str, str]:
    """立体部件应用到名称：酯在烷基词后插入，其余直接前缀（中文恒为前缀）。"""
    if (numbered.get("parent") or {}).get("kind") == "ester":
        return _ester_en_rs(en, parts), _with_rs(zh, parts)
    return _with_rs(en, parts), _with_rs(zh, parts)


def join_rs_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """对外 R/S 入口：算手性部件并应用到中英文名称（酯特殊插入）。"""
    rs = _rs_parts(numbered)
    if not rs:
        return en, zh
    return _apply_rs(en, zh, rs, numbered)
