"""L2 parent assembly / chain-walk / gate helpers (single authority).

Producers import helpers from here. parent_selector keeps FG try + select_parent
and may thin-re-export for compatibility.
"""
from __future__ import annotations

from namepredict.layer2.chain_walk import _longest_chain, _longest_from

# Re-export chain / gate primitives used by producers.
__all__ = [
    "_parent_dict",
    "_longest_from",
    "_longest_chain",
]


def _parent_dict(chain: list[int], kind: str, **kw) -> dict:
    return {"chain": chain, "n_carbons": len(chain), "kind": kind, **kw}
