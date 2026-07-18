"""Typed L3 substituent naming contract (types only for Task 3)."""
from __future__ import annotations

from dataclasses import dataclass

from namepredict.layer2.claimable_block import ClaimedBlock


@dataclass(frozen=True)
class SubstituentName:
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained", "rooted_tree", or "recursive"
