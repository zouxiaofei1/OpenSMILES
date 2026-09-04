"""L5 立体描述符前缀：E/Z（P-91.2/P-93.4）+ CIP R/S（P-92/P-93）；由 stereo_ez.py、stereo_rs.py 与 _stereo_common.py 合并而来，_split_stereo_lead 为共享立体块拆分器。"""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import BondStereo, ChiralType, Mol

from namepredict.layer1 import fg_registry as _fg_reg


def _split_stereo_lead(name: str) -> tuple[str, str]:
    """从名称/词干中切分前导 '(…)-' 立体块。"""
    if not name.startswith("("):
        return "", name
    close = name.find(")-")
    if close < 0:
        return "", name
    return name[: close + 2], name[close + 2 :]


# --- E/Z 立体描述符 -------------------------------------------------

def _stereo_tag(st) -> str:
    """将 RDKit 立体键枚举转成 E/Z 前缀标记。"""
    if st == BondStereo.STEREOE:
        return "(E)-"
    if st == BondStereo.STEREOZ:
        return "(Z)-"
    return ""


def _bond_stereo(mol: Mol | None, double_bond) -> str:
    """查指定双键的立体标签（无键/无立体返回空串）。"""
    if mol is None or not double_bond:
        return ""
    c1, c2 = double_bond
    bond = mol.GetBondBetweenAtoms(int(c1), int(c2))
    return _stereo_tag(bond.GetStereo()) if bond is not None else ""


def _ez_prefix(numbered: dict) -> str:
    """单双键母体的 E/Z 前缀（带位次，如 '(2E)-'；双键不在母体链上时退化为裸 '(E)-'）。"""
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


def _ez_letter(tag: str) -> str:
    """'(E)-' → 'E'；空 → ''。"""
    return tag[1] if len(tag) >= 3 and tag[0] == "(" else ""


def _ez_bond_part(mol, chain: list[int], bond) -> tuple[int, str] | None:
    """计算单条立体双键的 (位次, 字母) 部件。"""
    loc = _bond_min_loc(chain, bond)
    letter = _ez_letter(_bond_stereo(mol, bond))
    return (loc, letter) if loc is not None and letter else None


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


# --- CIP R/S 立体描述符 ----------------------------------------------

_RS_KINDS = _fg_reg.rs_fgs() | frozenset({"radical"})

def _assign_cip(mol: Mol) -> None:
    """强制重算分子立体化学（CIP 分配）；隐式 H 的 [C@]/[C@@] 手性碳先补显式 H 再赋。"""
    r = Chem.RWMol(mol)
    for a in list(r.GetAtoms()):
        if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED and a.GetTotalNumHs() == 0 and a.GetDegree() < 4:
            r.AddBond(a.GetIdx(), r.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
    m = r.GetMol()
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)
    for a in mol.GetAtoms():
        b = m.GetAtomWithIdx(a.GetIdx())
        if b.HasProp("_CIPCode"):
            a.SetProp("_CIPCode", b.GetProp("_CIPCode"))


def _cip_code(atom) -> str | None:
    """取原子的 CIP 代码，仅 R/S 有效时返回。"""
    if not atom.HasProp("_CIPCode"):
        return None
    code = atom.GetProp("_CIPCode")
    return code if code in ("R", "S") else None


def _cip_on_chain(mol: Mol, chain: list[int]) -> list[tuple[int, str]]:
    """返回母体链上手性中心的 (位次, R/S)。"""
    _assign_cip(mol)
    out: list[tuple[int, str]] = []
    for loc, idx in enumerate(chain, 1):
        code = _cip_code(mol.GetAtomWithIdx(int(idx)))
        if code:
            out.append((loc, code))
    return out


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


def _rs_parts(numbered: dict) -> list[tuple[int, str]]:
    """取母体上手性中心的 (位次, R/S) 列表（kind 不支持、折叠环或空链时为空）。

    链式主官能团母体按 kind ∈ _RS_KINDS 放行；环/稠合骨架母体以 scaffold_id 识别
    （其 chain 是 L4 定向编号的整环 walk，环上 sp3 手性中心可被 _cip_on_chain 扫到）。
    折叠环仍由 _collapsed_parent 跳过。
    """
    parent = numbered.get("parent") or {}
    kind = parent.get("kind")
    is_ring_parent = bool(parent.get("scaffold_id"))
    if (kind not in _RS_KINDS and not is_ring_parent) or _collapsed_parent(parent):
        return []
    mol, chain = parent.get("mol"), parent.get("chain") or []
    if mol is None or not chain:
        return []
    return _cip_on_chain(mol, chain)


def _parse_token(tok: str) -> tuple[int | None, str] | None:
    """解析单个 token：'E'→(None,'E')；'8R'→(8,'R')；'2E'→(2,'E')。"""
    if tok in ("E", "Z", "R", "S"):
        return None, tok
    i = 0
    while i < len(tok) and tok[i].isdigit():
        i += 1
    if i and tok[i:] in ("E", "Z", "R", "S"):
        return int(tok[:i]), tok[i:]
    return None


def _parse_stereo(tag: str) -> list[tuple[int | None, str]]:
    """将 '(E)-' / '(2E,6Z)-' / '(E,8R)-' 解析为 (loc, letter) 列表。"""
    if not tag.startswith("(") or not tag.endswith(")-"):
        return []
    raw = tag[1:-2]
    out: list[tuple[int | None, str]] = []
    for tok in raw.split(","):
        p = _parse_token(tok.strip())
        if p:
            out.append(p)
    return out


def _fmt_part(loc: int | None, letter: str) -> str:
    """单部件格式化：有位次拼 loc+字母，无位次只给字母。"""
    return letter if loc is None else f"{loc}{letter}"


def _format_stereo(parts: list[tuple[int | None, str]]) -> str:
    """将立体部件列表拼成 '(…)-' 前缀（无部件返回空串）。"""
    if not parts:
        return ""
    ordered = sorted(parts, key=lambda x: (x[0] is not None, x[0] or 0))
    body = ",".join(_fmt_part(loc, let) for loc, let in ordered)
    return f"({body})-"


def _merge_parts(
    old: list[tuple[int | None, str]], rs: list[tuple[int | None, str]],
) -> list[tuple[int | None, str]]:
    """保留非 R/S 立体；按位次添加 R/S。"""
    keep = [(loc, let) for loc, let in old if let not in ("R", "S")]
    return keep + list(rs)


def _with_rs(name: str, rs: list[tuple[int | None, str]]) -> str:
    """将 R/S 部件并入名称已有的立体前缀。"""
    if not rs:
        return name
    tag, stem = _split_stereo_lead(name)
    parts = _merge_parts(_parse_stereo(tag), rs)
    return f"{_format_stereo(parts)}{stem}"


def _ester_en_rs(en: str, rs: list[tuple[int | None, str]]) -> str:
    """在烷基词后插入 R/S：'methyl X' → 'methyl (2S)-X'。"""
    if not rs or " " not in en:
        return _with_rs(en, rs)
    alkyl, acyl = en.split(" ", 1)
    return f"{alkyl} {_with_rs(acyl, rs)}"


def apply_rs_prefix(numbered: dict, en: str, zh: str) -> tuple[str, str]:
    """对外 R/S 入口：算手性部件并应用到中英文名称（酯特殊插入）。"""
    rs_raw = _rs_parts(numbered)
    if not rs_raw:
        return en, zh
    kind = (numbered.get("parent") or {}).get("kind")
    rs = rs_raw
    if kind == "ester":
        return _ester_en_rs(en, rs), _with_rs(zh, rs)
    return _with_rs(en, rs), _with_rs(zh, rs)
