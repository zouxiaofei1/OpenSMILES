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

def iter_parent_candidates(info: dict) -> list[dict]:
    """Ranked parent candidates, each finalized with immutable owned_atoms."""
    from namepredict.layer2.candidates import _alkane_fallback, _collect_candidates
    cands = _collect_candidates(info) or [_alkane_fallback(info)]
    return _finalize_ranked(info, cands)

def select_parent(info: dict) -> dict:
    """Legacy: first finalized ranked candidate (or alkane fallback)."""
    return iter_parent_candidates(info)[0]
