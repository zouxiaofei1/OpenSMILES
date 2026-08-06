"""Claim-based side extraction via SubstituentNamer (coverage fill-in)."""
from __future__ import annotations

from namepredict.cache.common_names import CommonNameCache


def _claim_kind(slot_value: str) -> str:
    return {
        "ether_o": "alkoxy",
        "amide_n": "n_block",
        "ring_c": "alkyl",
        "chain_c": "alkyl",
    }.get(slot_value, "side")


def _kind_for_named(named) -> str:
    if named.en in ("methoxy", "ethoxy", "propoxy", "butoxy"):
        return "alkoxy"
    return _claim_kind(named.claim.slot.value)


def sub_from_named(named) -> dict:
    claim = named.claim
    return {
        "kind": _kind_for_named(named), "n_carbons": len(claim.atoms),
        "attach_idx": claim.attach_parent, "atoms": sorted(claim.atoms),
        "en": named.en, "zh": named.zh, "paren": named.requires_parentheses,
        "backend": named.backend,
    }


def _covered_atoms(subs: list[dict]) -> set[int]:
    out: set[int] = set()
    for s in subs:
        out.update(s.get("atoms") or [])
    return out


def _should_skip(claim, covered: set[int]) -> bool:
    from namepredict.layer3.claimable_block import SideSlot
    return claim.slot == SideSlot.AMIDE_N or bool(set(claim.atoms) & covered)


def _append_named(mol, claim, namer, covered: set[int], out: list[dict]) -> None:
    if _should_skip(claim, covered):
        return
    named = namer.name(mol, claim, depth=0)
    if named is None:
        return
    out.append(sub_from_named(named))
    covered |= set(named.claim.atoms)


def _named_new_sides(mol, owned, covered: set[int], *, name_mode: str = "general", cache: CommonNameCache | None = None) -> list[dict]:
    from namepredict.layer3.claimable_block import iter_claims
    from namepredict.layer3.substituent_namer import SubstituentNamer

    namer, out = SubstituentNamer(name_mode=name_mode, cache=cache), []
    for claim in iter_claims(mol, owned):
        _append_named(mol, claim, namer, covered, out)
    return out


def extract_claimed_sides(info: dict, parent: dict, existing: list[dict], *, name_mode: str = "general", cache: CommonNameCache | None = None) -> list[dict]:
    """Name ownership-boundary claims not already covered by legacy extractors."""
    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    return _named_new_sides(info["mol"], owned, _covered_atoms(existing), name_mode=name_mode, cache=cache)
