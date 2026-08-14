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


# Anchored-table / retained leaves whose name implies a non-alkyl kind, so the
# L5 layers keying off `kind` (iso fusion, poly-anisole, halo benzene) fire.
_NAME_KIND = {
    "fluoro": "halo", "chloro": "halo", "bromo": "halo", "iodo": "halo",
    "nitro": "nitro",
    "isocyanato": "isocyanato", "isothiocyanato": "isothiocyanato",
}


def _kind_for_named(named) -> str:
    if named.en in ("methoxy", "ethoxy", "propoxy", "butoxy"):
        return "alkoxy"
    return _NAME_KIND.get(named.en, _claim_kind(named.claim.slot.value))


def sub_from_named(named, mol=None) -> dict:
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
    out: set[int] = set()
    for s in subs:
        out.update(s.get("atoms") or [])
    return out


def _should_skip(claim, covered: set[int]) -> bool:
    from namepredict.layer3.claimable_block import SideSlot
    return claim.slot == SideSlot.AMIDE_N or bool(set(claim.atoms) & covered)


# Ester acid-side O (alkoxy arm): a claimed side attached to an O atom of an
# ester/benzoate parent is the O-side alkyl, consumed by L5 join_kind_name.
_ESTER_O_SIDE_KINDS = frozenset({"ester", "diester", "benzoate"})


def _append_named(mol, claim, namer, covered: set[int], out: list[dict], *, o_side: bool = False) -> None:
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
    from namepredict.layer3.claimable_block import iter_claims
    from namepredict.layer3.substituent_namer import SubstituentNamer

    namer, out = SubstituentNamer(name_mode=name_mode, cache=cache), []
    for claim in iter_claims(mol, owned):
        _append_named(mol, claim, namer, covered, out, o_side=o_side)
    return out


def extract_claimed_sides(info: dict, parent: dict, existing: list[dict], *, name_mode: str = "general", cache: CommonNameCache | None = None) -> list[dict]:
    """Name ownership-boundary claims not already covered by legacy extractors."""
    owned = parent.get("owned_atoms")
    if owned is None:
        return []
    # benzoate（kind='benzene' + ester FG）靠 o_idx 字段识别 O-side；链状 ester 走 kind 表。
    o_side = parent.get("kind") in _ESTER_O_SIDE_KINDS or parent.get("o_idx") is not None
    return _named_new_sides(info["mol"], owned, _covered_atoms(existing), name_mode=name_mode, cache=cache, o_side=o_side)
