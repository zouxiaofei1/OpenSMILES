"""可组合的母体候选评分（简化的 IUPAC P-44 优先规则）：每个候选母体 dict 被评分为元组（越大越优）(has_principal_fg, fg_class_rank, sides_ok, is_hetero_ring, is_carbo_ring, n_rings, ring_size, retained_bonus, n_unsat, n_carbons, -n_unhandled_side)。"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr
from namepredict.layer2.parent_candidate import principal_key, with_principal_group_contract

def _n_unhandled(parent: dict) -> int:
    """取候选的未处理侧链数（缺省 0）。"""
    return int(parent.get("n_unhandled") or 0)

def _sides_ok(parent: dict) -> int:
    """无未处理侧链为 1，否则 0。"""
    return 0 if _n_unhandled(parent) else 1

def _is_hetero_ring(kind: str) -> int:
    """kind 是否为杂环（1/0）。"""
    return _kr.is_hetero_ring(kind)

def _is_carbo_ring(kind: str) -> int:
    """kind 是否为碳环（1/0）。"""
    return _kr.is_carbo_ring(kind)

def _n_rings(kind: str) -> int:
    """取 kind 的环数。"""
    return _kr.n_rings_of(kind)

def _ring_size(parent: dict, kind: str) -> int:
    """环骨架的环大小（链长兜底）。"""
    if not _n_rings(kind):
        return 0
    return len(parent.get("chain") or [])

def _n_unsat(parent: dict) -> int:
    """统计候选的不饱和键数。"""
    dbs = parent.get("double_bonds")
    if dbs:
        return len(dbs)
    return 1 if parent.get("double_bond") or parent.get("triple_bond") else 0

def _n_carbons(parent: dict) -> int:
    """取候选碳数（缺省按链长）。"""
    return int(parent.get("n_carbons") or len(parent.get("chain") or []))

def _retained_bonus(kind: str) -> int:
    """kind 的保留名加分（1/0）。"""
    return _kr.retained_bonus(kind)

def _p44_1_1(parent: dict) -> tuple[int, int]:
    """P-44.1.1 键：(FG 等级, 基团个数)。"""
    facts = principal_key(with_principal_group_contract(parent))
    return facts.principal_group_class, facts.principal_group_count

def _later_score(parent: dict, kind: str) -> tuple:
    """P-44.1.1 之后的比较键（侧链/环性/不饱和等）。"""
    return (_sides_ok(parent), _is_hetero_ring(kind), _is_carbo_ring(kind),
            _n_rings(kind), _ring_size(parent, kind), _retained_bonus(kind),
            _n_unsat(parent), _n_carbons(parent), -_n_unhandled(parent))

def _score_parent(info: dict, parent: dict) -> tuple:
    """组合全部打分项，返回候选母体评分元组。"""
    kind = parent.get("kind") or ""
    return (*_p44_1_1(parent), *_later_score(parent, kind))

