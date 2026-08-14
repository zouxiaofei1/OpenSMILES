"""可组合的母体候选评分（简化的 IUPAC P-44 优先规则）：每个候选母体 dict 被评分为元组（越大越优）(has_principal_fg, fg_class_rank, sides_ok, is_hetero_ring, is_carbo_ring, n_rings, ring_size, retained_bonus, n_unsat, n_carbons, -n_unhandled_side)。"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr
from namepredict.layer2.parent_candidate import principal_key, with_principal_group_contract

def _n_unhandled(parent: dict) -> int:
    return int(parent.get("n_unhandled") or 0)

def _sides_ok(parent: dict) -> int:
    return 0 if _n_unhandled(parent) else 1

def _is_hetero_ring(kind: str) -> int:
    return _kr.is_hetero_ring(kind)

def _is_carbo_ring(kind: str) -> int:
    return _kr.is_carbo_ring(kind)

def _n_rings(kind: str) -> int:
    return _kr.n_rings_of(kind)

def _ring_size(parent: dict, kind: str) -> int:
    if not _n_rings(kind):
        return 0
    return len(parent.get("chain") or [])

def _n_unsat(parent: dict) -> int:
    dbs = parent.get("double_bonds")
    if dbs:
        return len(dbs)
    return 1 if parent.get("double_bond") or parent.get("triple_bond") else 0

def _n_carbons(parent: dict) -> int:
    return int(parent.get("n_carbons") or len(parent.get("chain") or []))

def _retained_bonus(kind: str) -> int:
    return _kr.retained_bonus(kind)

def _p44_1_1(parent: dict) -> tuple[int, int]:
    facts = principal_key(with_principal_group_contract(parent))
    return facts.principal_group_class, facts.principal_group_count

def _later_score(parent: dict, kind: str) -> tuple:
    return (_sides_ok(parent), _is_hetero_ring(kind), _is_carbo_ring(kind),
            _n_rings(kind), _ring_size(parent, kind), _retained_bonus(kind),
            _n_unsat(parent), _n_carbons(parent), -_n_unhandled(parent))

def _score_parent(info: dict, parent: dict) -> tuple:
    kind = parent.get("kind") or ""
    return (*_p44_1_1(parent), *_later_score(parent, kind))

