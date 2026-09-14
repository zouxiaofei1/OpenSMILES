"""基于 claim 的侧链提取，经 SubstituentNamer 填充。"""
from __future__ import annotations

from namepredict.tools.common_names import CommonNameCache
from namepredict.constants import (
    CLAIM_KIND, ESTER_O_SIDE_KINDS, NAME_KIND, N_PREFIX_KINDS,
)


def _claim_kind(slot_value: str) -> str:
    """将槽位值映射为取代基 kind。"""
    return CLAIM_KIND.get(slot_value, "side")


def _kind_for_named(named) -> str:
    """由命名结果确定取代基 kind（优先特殊保留名）。"""
    return NAME_KIND.get(named.en, _claim_kind(named.claim.slot.value))


def _is_ring_attach(mol, idx: int) -> bool:
    """附着原子是否为环员（环 N 用环上位次定位，不写 N- 前缀）。"""
    try:
        return mol.GetAtomWithIdx(int(idx)).IsInRing()
    except Exception:  # noqa: BLE001
        return False


def sub_from_named(named, mol=None) -> dict:
    """将命名结果封装为取代基字典。"""
    claim = named.claim
    if mol is not None:
        n_carbons = sum(1 for i in claim.atoms
                        if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
    else:
        n_carbons = len(claim.atoms)
    kind = _kind_for_named(named)
    if kind in N_PREFIX_KINDS and mol is not None and _is_ring_attach(mol, claim.attach_parent):
        kind = _claim_kind("ring_c")  # 环氮（内酰胺/环胺母体的环员 N）有环上位次可用，改用位次定位而非 N- 前缀
    return {
        "kind": kind, "n_carbons": n_carbons,
        "attach_idx": claim.attach_parent, "atoms": sorted(claim.atoms),
        "en": named.en, "zh": named.zh, "paren": named.requires_parentheses,
    }


def _covered_atoms(subs: list[dict]) -> set[int]:
    """收集已覆盖的原子索引集合。"""
    out: set[int] = set()
    for s in subs:
        out.update(s.get("atoms") or [])
    return out

def _append_named(mol, claim, namer, covered: set[int], out: list[dict], *, o_side: bool = False) -> None:
    """为单个 claim 命名并追加到输出（可标记 O 侧）。"""
    named = namer.name(mol, claim)
    if named is None:
        return
    s = sub_from_named(named, mol)
    if o_side and mol.GetAtomWithIdx(claim.attach_parent).GetAtomicNum() == 8:
        s["o_side"] = True
    out.append(s)
    covered |= set(named.claim.atoms)


def _named_new_sides(mol, owned, covered: set[int], *, cache: CommonNameCache | None = None, o_side: bool = False, root_ctx: tuple | None = None) -> list[dict]:
    """为所有权边界的所有 claim 生成命名侧链。"""
    from namepredict.layer3.claimable_block import iter_claims
    from namepredict.layer3.substituent_namer import SubstituentNamer

    namer, out = SubstituentNamer(cache=cache, root_ctx=root_ctx), []
    for claim in iter_claims(mol, owned):
        _append_named(mol, claim, namer, covered, out, o_side=o_side)
    return out


def extract_claimed_sides(info: dict, parent: dict, existing: list[dict], *, cache: CommonNameCache | None = None) -> list[dict]:
    """ claim 命名。"""
    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    o_side = parent.get("kind") in ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None  # benzoate（苯 base + ester FG）靠 o_idx 字段识别 O-side；链状 ester 走 kind 表。
    return _named_new_sides(info["mol"], owned, _covered_atoms(existing), cache=cache, o_side=o_side, root_ctx=info.get("root_ctx"))
