"""基于 claim 的侧链提取，经由 SubstituentNamer（覆盖填充）。"""
from __future__ import annotations

from namepredict.cache.common_names import CommonNameCache


def _claim_kind(slot_value: str) -> str:
    """将槽位值映射为取代基 kind。"""
    return {
        "ether_o": "alkoxy",
        "amide_n": "n_block",
        "ring_c": "alkyl",
        "chain_c": "alkyl",
    }.get(slot_value, "side")


# 锚定表/保留叶子的名称暗含非烷基 kind，使 L5 以 `kind` 键（iso 稠合、聚茴香醚、卤代苯）触发。
_NAME_KIND = {
    "fluoro": "halo", "chloro": "halo", "bromo": "halo", "iodo": "halo",
    "nitro": "nitro",
    "isocyanato": "isocyanato", "isothiocyanato": "isothiocyanato",
}


def _kind_for_named(named) -> str:
    """由命名结果确定取代基 kind（优先特殊保留名）。"""
    if named.en in ("methoxy", "ethoxy", "propoxy", "butoxy"):
        return "alkoxy"
    return _NAME_KIND.get(named.en, _claim_kind(named.claim.slot.value))


def sub_from_named(named, mol=None) -> dict:
    """将命名结果封装为取代基字典。"""
    claim = named.claim
    if mol is not None:
        n_carbons = sum(1 for i in claim.atoms
                        if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
    else:
        n_carbons = len(claim.atoms)
    return {
        "kind": _kind_for_named(named), "n_carbons": n_carbons,
        "attach_idx": claim.attach_parent, "atoms": sorted(claim.atoms),
        "en": named.en, "zh": named.zh, "paren": named.requires_parentheses,
        "backend": named.backend,
    }


def _covered_atoms(subs: list[dict]) -> set[int]:
    """收集已覆盖的原子索引集合。"""
    out: set[int] = set()
    for s in subs:
        out.update(s.get("atoms") or [])
    return out


def _should_skip(claim, covered: set[int]) -> bool:
    """判断 claim 是否应跳过（酰胺 N 或原子已被覆盖）。"""
    from namepredict.layer3.claimable_block import SideSlot
    return claim.slot == SideSlot.AMIDE_N or bool(set(claim.atoms) & covered)


# 酯酸侧 O（烷氧基臂）：连在酯 parent 的 O 原子上的侧链是 O 侧烷基，由 L5 join_kind_name 消费。
_ESTER_O_SIDE_KINDS = frozenset({"ester"})


def _append_named(mol, claim, namer, covered: set[int], out: list[dict], *, o_side: bool = False) -> None:
    """为单个 claim 命名并追加到输出（可标记 O 侧）。"""
    if _should_skip(claim, covered):
        return
    named = namer.name(mol, claim, depth=0)
    if named is None:
        return
    s = sub_from_named(named, mol)
    if o_side and mol.GetAtomWithIdx(claim.attach_parent).GetAtomicNum() == 8:
        s["o_side"] = True
    out.append(s)
    covered |= set(named.claim.atoms)


def _named_new_sides(mol, owned, covered: set[int], *, name_mode: str = "general", cache: CommonNameCache | None = None, o_side: bool = False) -> list[dict]:
    """为所有权边界的所有 claim 生成命名侧链。"""
    from namepredict.layer3.claimable_block import iter_claims
    from namepredict.layer3.substituent_namer import SubstituentNamer

    namer, out = SubstituentNamer(name_mode=name_mode, cache=cache), []
    for claim in iter_claims(mol, owned):
        _append_named(mol, claim, namer, covered, out, o_side=o_side)
    return out


def extract_claimed_sides(info: dict, parent: dict, existing: list[dict], *, name_mode: str = "general", cache: CommonNameCache | None = None) -> list[dict]:
    """为尚未被旧提取器覆盖的所有权边界 claim 命名。"""
    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    # benzoate（苯 base + ester FG）靠 o_idx 字段识别 O-side；链状 ester 走 kind 表。
    o_side = parent.get("kind") in _ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None
    return _named_new_sides(info["mol"], owned, _covered_atoms(existing), name_mode=name_mode, cache=cache, o_side=o_side)
