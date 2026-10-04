"""L3 取代基提取器：claim 枚举 → 递归命名 → 取代基 dict。
claim 枚举是本层唯一提取路径，未命名成功的 claim 静默跳过。
"""
from __future__ import annotations

from rdkit.Chem import BondType

from opensmiles.tools.common_names import CommonNameCache
from opensmiles.constants import (
    C, CLAIM_KIND, ESTER_O_SIDE_KINDS, N, NAME_KIND, N_PREFIX_KINDS,
)


def _claim_kind(slot_value: str) -> str:
    """将槽位值映射为取代基 kind。"""
    return CLAIM_KIND.get(slot_value, "side")


def _amidine_n_owned(mol, parent: dict | None) -> frozenset[int]:
    """脒/胍母体（单碳 C 带 =N 且另有单键 N）的单键氮并入所有权边界。"""
    if parent is None:
        return frozenset()
    chain = list(parent.get("chain") or ())
    if len(chain) != 1:
        return frozenset()
    atom = mol.GetAtomWithIdx(int(chain[0]))
    if atom.GetAtomicNum() != C:
        return frozenset()
    sgl_n: list[int] = []
    has_dbl_n = False
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() != N:
            continue
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
        if bond.GetBondType() == BondType.DOUBLE:
            has_dbl_n = True
        elif bond.GetBondType() == BondType.SINGLE:
            sgl_n.append(nb.GetIdx())
    if not has_dbl_n or not sgl_n:  # 非脒/胍中心：不动边界
        return frozenset()
    return frozenset(sgl_n)  # 亚胺氮属母体：其上的臂须从 N 外侧键起切


def _hydrazide_n_owned(mol, parent: dict | None) -> frozenset[int]:
    """酰肼母体（P-66.3.0）：远端 N 并入所有权边界并改记 N′。"""
    idx = (parent or {}).get("hydrazide_n_idx")
    return frozenset({int(idx)}) if idx is not None else frozenset()


def sub_from_named(named, mol, parent: dict | None = None) -> dict:
    """将命名结果封装为取代基字典。"""
    claim = named.claim
    n_carbons = sum(1 for i in claim.atoms
                    if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
    kind = NAME_KIND.get(named.en, _claim_kind(named.claim.slot.value))  # 优先特殊保留名
    attach = mol.GetAtomWithIdx(int(claim.attach_parent))
    in_parent = parent is not None and int(claim.attach_parent) in set(parent.get("chain") or ())
    if kind in N_PREFIX_KINDS and (attach.IsInRing() or in_parent):  # 附着原子为环员或母体骨架成员：有数字位次可用，不写 N- 前缀
        kind = _claim_kind("ring_c")  # 环氮与均一氮链改用位次定位而非 N- 前缀
    return {
        "kind": kind, "n_carbons": n_carbons,
        "attach_idx": claim.attach_parent, "atoms": sorted(claim.atoms),
        "en": named.en, "zh": named.zh, "paren": named.requires_parentheses,
    }


def _append_named(mol, claim, namer, out: list[dict], *, o_side: bool = False, side_z: int = 8,
                  parent: dict | None = None) -> None:
    """为单个 claim 命名并追加到输出（可标记 O 侧；硫代酯的侧臂元素为 S）。"""
    named = namer.name(mol, claim)
    if named is None:
        return
    s = sub_from_named(named, mol, parent)
    if o_side and mol.GetAtomWithIdx(claim.attach_parent).GetAtomicNum() == side_z:
        s["o_side"] = True
    out.append(s)


def _side_arm_claim(mol, claim, chain: frozenset[int]):
    """把切在所有权内非链杂原子上的侧臂并回该杂原子、改挂到链上原子。"""
    from dataclasses import replace

    a = int(claim.attach_parent)
    if a in chain:
        return claim
    atom = mol.GetAtomWithIdx(a)
    if atom.GetAtomicNum() not in (8, 16) or atom.GetDegree() < 2:
        return claim
    linked = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in chain]
    if len(linked) != 1:
        return claim
    return replace(claim, attach_parent=linked[0], root=a, atoms=frozenset(set(claim.atoms) | {a}))


def extract_substituents(info: dict, parent: dict, *, cache: CommonNameCache | None = None) -> list[dict]:
    """L3 入口：为所有权边界内的每个 claim 命名并封装为取代基。"""
    from opensmiles.layer3.claimable_block import iter_claims
    from opensmiles.layer3.substituent_namer import SubstituentNamer

    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    mol = info["mol"]
    o_side = parent.get("kind") in ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None
    side_z = 16 if parent.get("thio_side") else 8  # 硫代酯的酯侧臂元素为 S（P-65.6.3.3.7.1）
    chain = frozenset(parent.get("chain") or ())
    namer, out = SubstituentNamer(cache=cache, root_ctx=info.get("root_ctx")), []
    cut_owned = (frozenset(owned) | _amidine_n_owned(mol, parent)  # 胍 N 上的臂从 N 外侧键起切
                 | _hydrazide_n_owned(mol, parent))  # 酰肼远端 N 同理（P-66.3.0）
    for claim in iter_claims(mol, cut_owned):
        if not (o_side and mol.GetAtomWithIdx(claim.attach_parent).GetAtomicNum() == side_z):
            claim = _side_arm_claim(mol, claim, chain)  # 侧臂切在桥杂原子上会丢臂：并入桥原子改挂链上
        _append_named(mol, claim, namer, out, o_side=o_side, side_z=side_z, parent=parent)
    return _mark_hydrazide_primes(out, parent)


def _mark_hydrazide_primes(subs: list[dict], parent: dict | None) -> list[dict]:
    """酰肼（P-66.3.3.1）两端 N 取代基定死撇号：羰基侧 N，远端 N′。"""
    hn = (parent or {}).get("hydrazide_n_idx")
    if hn is None:
        return subs
    hnear = (parent or {}).get("hydrazide_near_n_idx")
    hi, lo = int(hn), (int(hnear) if hnear is not None else None)
    return [{**s, "n_prime": 1} if s.get("attach_idx") == hi
            else ({**s, "n_prime": 0} if s.get("attach_idx") == lo else s)
            for s in subs]
