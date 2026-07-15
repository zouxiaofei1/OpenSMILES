"""Ring numbering candidate generation (L4 pure)."""
from __future__ import annotations


def _rotations(chain: list[int]) -> list[list[int]]:
    return [chain[i:] + chain[:i] for i in range(len(chain))]


def ring_candidates(chain: list[int]) -> list[list[int]]:
    """All rotation × reverse candidates (6×2=12 for n=6)."""
    out: list[list[int]] = []
    for base in (list(chain), list(reversed(chain))):
        out.extend(_rotations(base))
    return out


def labels_for(n: int) -> tuple[str, ...]:
    """Sequential locant labels '1'..'n'."""
    return tuple(str(i) for i in range(1, n + 1))


def candidates_for(mode: str, chain: list[int]) -> list[list[int]]:
    """Dispatch candidate generation by NumberingPolicy mode."""
    if mode in ("carbocycle_free", "poly_unsat", "multi_hetero"):
        return ring_candidates(chain)
    return [list(chain)]
