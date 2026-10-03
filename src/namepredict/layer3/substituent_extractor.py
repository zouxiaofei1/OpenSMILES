"""L3 取代基提取器：claim 枚举 → 递归命名 → 取代基 dict。
claim 枚举是本层唯一提取路径，未命名成功的 claim 静默跳过。
"""
from __future__ import annotations

from namepredict.tools.common_names import CommonNameCache
from namepredict.constants import (
    CLAIM_KIND, ESTER_O_SIDE_KINDS, NAME_KIND, N_PREFIX_KINDS,
)


def _claim_kind(slot_value: str) -> str:
    """将槽位值映射为取代基 kind。"""
    return CLAIM_KIND.get(slot_value, "side")


def sub_from_named(named, mol, parent: dict | None = None) -> dict:
    """将命名结果封装为取代基字典。"""
    claim = named.claim
    n_carbons = sum(1 for i in claim.atoms
                    if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
    kind = NAME_KIND.get(named.en, _claim_kind(named.claim.slot.value))  # 优先特殊保留名
    attach = mol.GetAtomWithIdx(int(claim.attach_parent))
    in_parent = parent is not None and int(claim.attach_parent) in set(parent.get("chain") or ())
    if kind in N_PREFIX_KINDS and (attach.IsInRing() or in_parent):  # 附着原子为环员或母体骨架成员：有数字位次可用，不写 N- 前缀
        kind = _claim_kind("ring_c")  # 环氮（内酰胺/环胺母体的环员 N）与均一氮链（1-methyltriazane）改用位次定位而非 N- 前缀
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


def extract_substituents(info: dict, parent: dict, *, cache: CommonNameCache | None = None) -> list[dict]:
    """L3 入口：为所有权边界内的每个 claim 命名并封装为取代基。"""
    from namepredict.layer3.claimable_block import iter_claims
    from namepredict.layer3.substituent_namer import SubstituentNamer

    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    mol = info["mol"]
    # benzoate（苯 base + ester FG）靠 o_idx 字段识别 O-side；链状 ester 走 kind 表。
    o_side = parent.get("kind") in ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None
    side_z = 16 if parent.get("thio_side") else 8  # 硫代酯的酯侧臂元素为 S（P-65.6.3.3.7.1）
    namer, out = SubstituentNamer(cache=cache, root_ctx=info.get("root_ctx")), []
    for claim in iter_claims(mol, owned):
        _append_named(mol, claim, namer, out, o_side=o_side, side_z=side_z, parent=parent)
    return out
