"""P-24 螺环母体名：von Baeyer 螺描述符 + 组分式螺环各组分名与连接位次。

单环组分出 spiro[..] 描述符 + 烷词干；组分式螺环按 P-24.5.1 拼
`spiro[组分1-位次,位次′-组分2-…]`，多螺用 dispiro/trispiro。L5 不得 import L2，
node 只按鸭子类型读 descriptor / components / links。
"""
from __future__ import annotations

from dataclasses import replace

from rdkit import Chem

from namepredict.constants import BIS_EN, BIS_ZH, C, MULT_EN, MULT_ZH
from namepredict.layer1.ring_systems import kekulized
from namepredict.layer4.locant_calc import locant_key
from namepredict.layer5.skeleton_replacement import prefix_from_chain, skeleton_replacement_prefix
from namepredict.layer5.stems import alkane_en, alkane_zh, stem_forms

APOSTROPHE = "'"


def spiro_multiplier(n_spiro: int) -> tuple[str, str] | None:
    """螺原子数词头：1 用 spiro/螺，2 起用 dispiro/二螺 型计数词（不是 bispiro）。"""
    if n_spiro == 1:
        return "spiro", "螺"
    en, zh = MULT_EN.get(n_spiro), MULT_ZH.get(n_spiro)
    return (f"{en}spiro", f"{zh}螺") if en and zh else None


def spiro_descriptor_str(descriptor: tuple[int, ...], superscripts: tuple[int, ...]) -> str:
    """螺描述符串：重访螺原子的位次紧接段长之后（判分口径无 ^{}）。"""
    parts: list[str] = []
    for i, d in enumerate(descriptor):
        sup = superscripts[i] if i < len(superscripts) else 0
        parts.append(f"{d}{sup}" if sup else str(d))
    return ".".join(parts)


def spiro_parent_names(mol, node, chain: list[int]):
    """螺环名 → ((完整英, 完整中), (裸词干英, 裸词干中))；不可组装返回 None。"""
    if not chain or not node.descriptor:
        return None
    mult = spiro_multiplier(len(node.free_spiro_atoms))
    if mult is None:
        return None
    a_en, a_zh = prefix_from_chain(mol, chain)
    if a_en is None or a_zh is None:
        return None
    stem_en, stem_zh = alkane_en(len(chain)), alkane_zh(len(chain))
    if not stem_en or not stem_zh:
        return None
    desc = spiro_descriptor_str(node.descriptor, node.descriptor_superscripts)
    body_en, body_zh = f"{mult[0]}[{desc}]", f"{mult[1]}[{desc}]"
    return stem_forms(f"{a_en}{body_en}", f"{a_zh}{body_zh}", stem_en, stem_zh)


# ── P-24.5~24.7 组分式螺环母体名 ────────────────────────

def fbs_parent_name(mol, node) -> tuple[str, str] | None:
    """组分式螺环母体名 (en, zh)；不可组装返回 None。"""
    comps = list(getattr(node, "components", ()) or ())
    links = list(getattr(node, "links", ()) or ())
    cites = list(getattr(node, "cites", ()) or ())
    if mol is None or len(comps) < 2 or not cites or len(links) != len(comps) - 1:
        return None
    mult = spiro_multiplier(len(comps) - 1)
    if mult is None:
        return None
    kek = kekulized(mol) or mol
    names: list[tuple[str, str]] = []
    a_en: list[str] = []
    a_zh: list[str] = []
    for comp in comps:  # 组分名与 'a' 前缀都按本组分自身编号渲染
        num = _numbering(comp)
        if num is None:
            return None
        name = _component_name(comp, num, kek)
        pre = _a_prefix(mol, comp, num)
        if name is None or pre is None:
            return None
        names.append(name)
        if pre[0]:
            a_en.append(pre[0])
            a_zh.append(pre[1])
    head = _cite_names(cites[0], names)
    if head is None:
        return None
    body_en, body_zh = head
    link_of = {slot: (parent, spiro) for spiro, parent, slot in links}
    seen = set(cites[0])
    for k in range(1, len(cites)):  # 相邻引用项之间写螺位次对（组内多对用冒号连接）
        items = []
        for slot in cites[k]:  # 父槽位须已列出（支链的第三个端环接在中心组上）
            parent, spiro = link_of.get(slot, (None, None))
            if parent not in seen or spiro is None:
                return None
            items.append((parent, slot, spiro))
        pairs = []
        for parent, slot, spiro in sorted(items):  # 每对 = (前列组分位次, 本次组分位次)
            la, lb = _locant(comps[parent], spiro), _locant(comps[slot], spiro)
            if la is None or lb is None:
                return None
            pairs.append(f"{la}{APOSTROPHE * parent},{lb}{APOSTROPHE * slot}")
        label = _cite_names(cites[k], names)
        if label is None:
            return None
        body_en = f"{body_en}-{':'.join(pairs)}-{label[0]}"
        body_zh = f"{body_zh}-{':'.join(pairs)}-{label[1]}"
        seen.update(cites[k])
    return (f"{'-'.join(a_en)}{mult[0]}[{body_en}]",
            f"{'-'.join(a_zh)}{mult[1]}[{body_zh}]")


def _cite_names(cite, names) -> tuple[str, str] | None:
    """引用项的双语名：多项时为 bis(...)/tris(...) 倍增词组（P-24.7.1/24.7.2）。"""
    if len(cite) == 1:
        return names[cite[0]]
    en, zh = BIS_EN.get(len(cite)), BIS_ZH.get(len(cite))
    if en is None or zh is None:
        return None
    return (f"{en}({names[cite[0]][0]})", f"{zh}({names[cite[0]][1]})")


def _numbering(comp):
    """L4 已选定的组分编号（无候选返回 None）。"""
    return comp.numberings[0] if comp.numberings else None


def _locant(comp, atom) -> str | None:
    """组分内某原子的本位次标签（不含撇号）。"""
    num = _numbering(comp)
    return None if num is None else num.locants.get(atom)


def _component_name(comp, num, kek) -> tuple[str, str] | None:
    """组分全名：保留名/稠合名原样；可接裸词干者按本组分环内不饱和补词尾。"""
    if comp.bare_en is None:
        return (comp.base_en, comp.base_zh)
    ene, yne = _unsat_locants(comp, num, kek)
    if not ene and not yne:
        return (comp.base_en, comp.base_zh)
    from namepredict.layer5.chain_engine import _KIND_TABLE, _chain_names
    spec = replace(_KIND_TABLE["alkane"], stem=(comp.bare_en, comp.bare_zh), coda="")
    numbered = {"parent": {"scaffold_id": None},
                "ene_locants": ene or None, "yne_locants": yne or None}
    return _chain_names(spec, len(comp.atom_ids), numbered)


def _unsat_locants(comp, num, kek) -> tuple[list, list]:
    """组分内环上双/三键位次（本位次，取编号方向较低端）。"""
    loc = num.locants
    order = {a: i for i, a in enumerate(num.chain)}
    kinds = {Chem.BondType.DOUBLE: [], Chem.BondType.TRIPLE: []}
    for bond in kek.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i not in loc or j not in loc or not bond.IsInRing() or bond.GetBondType() not in kinds:
            continue
        kinds[bond.GetBondType()].append(loc[num.chain[min(order[i], order[j])]])
    return (sorted(kinds[Chem.BondType.DOUBLE], key=locant_key),
            sorted(kinds[Chem.BondType.TRIPLE], key=locant_key))


def _a_prefix(mol, comp, num) -> tuple[str, str] | None:
    """组分杂原子的 'a' 前缀；P-24.5.2 要求挂在 spiro 之前，位次带本组分撇号。"""
    if comp.kind != "bridged_ring":
        return "", ""  # 保留名/稠合名已含杂原子，无需 'a' 前缀
    by_z: dict[int, list[int]] = {}
    for a in comp.atom_ids:
        z = mol.GetAtomWithIdx(a).GetAtomicNum()
        if z != C:
            by_z.setdefault(z, []).append(int(num.locants[a]))
    if not by_z:
        return "", ""
    pre = skeleton_replacement_prefix(by_z, APOSTROPHE * comp.index)
    return None if pre is None else pre
