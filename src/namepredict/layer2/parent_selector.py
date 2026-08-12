from __future__ import annotations


def _rank_candidates(info: dict, cands: list[dict]) -> list[dict]:
    from namepredict.layer2.scoring import _score_parent
    return sorted(cands, key=lambda c: _score_parent(info, c), reverse=True)

def _finalize_ranked(info: dict, cands: list[dict]) -> list[dict]:
    from namepredict.layer2.kind_registry import pack_parent_stem
    from namepredict.layer2.parent_candidate import with_principal_group_contract
    from namepredict.layer2.parent_ownership import finalize_parent_ownership
    mol = info.get("mol")
    return [
        finalize_parent_ownership(
            pack_parent_stem(with_principal_group_contract(c), mol), mol,
        )
        for c in _rank_candidates(info, cands)
    ]

def select_parent(info: dict, *, all_candidates: bool = False) -> dict | list[dict] | None:
    """Ranked parent candidates, each finalized with immutable owned_atoms.
    """
    from namepredict.layer2.candidates import _collect_candidates
   
    cands = _finalize_ranked(info, _collect_candidates(info))
    if all_candidates:
        return cands
    return next(iter(cands), None)
