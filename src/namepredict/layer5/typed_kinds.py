from __future__ import annotations

# 稠环/杂环 scaffold：+FG 时返回 FG 类别，词干由 chain_engine 注入（正交化）。
_RING_FG_SCAFFOLDS = frozenset({
    "naphthalene", "indole", "pyridine", "quinoline",
})


def _scaffold(numbered: dict) -> str | None:
    """Parent scaffold id (None for open chain). L2 kind 的环系语义不再被 L5 依赖。"""
    return (numbered.get("parent") or {}).get("scaffold_id")


def _typed_acid_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "acid":
        return kind
    # 环酸（carbocycle/稠环）exocyclic：kind 收敛为 FG 类别 acid，"carboxylic acid" 后缀由 chain_engine 组装（正交化，不再枚举 cycloalkanecarboxylic）。
    if facts.relation.value == "exocyclic" and (
        _scaffold(numbered) == "carbocycle" or _scaffold(numbered) in _RING_FG_SCAFFOLDS
    ):
        return "acid"
    if facts.relation.value == "exocyclic":
        return kind
    return "acid"  # 链式:数量统一由 multiplicity 承载,kind 恒为 acid


def _typed_ester_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ester":
        return kind
    return kind


def _typed_aldehyde_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "aldehyde":
        return kind
    return kind


def _typed_nitrile_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "nitrile":
        return kind
    return kind


def _typed_amide_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amide":
        return kind
    return kind

def _typed_ketone_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ketone":
        return kind
    # 数量由 multiplicity 承载，kind 恒为 ketone；dione 后缀由 chain_engine variant 切换（对齐 alcohol）。
    return "ketone"

def _typed_ring_alcohol_kind(kind: str, numbered: dict, facts) -> str | None:
    sid = _scaffold(numbered)
    if sid == "carbocycle" or sid in _RING_FG_SCAFFOLDS:
        return "alcohol"
    return None

def _typed_alcohol_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "alcohol":
        return kind
    ring_kind = _typed_ring_alcohol_kind(kind, numbered, facts)
    if ring_kind:
        return ring_kind
    if kind not in {"alcohol", "diol", "triol"}:
        return kind
    return "alcohol"  # 链式:数量统一由 multiplicity 承载,kind 恒为 alcohol

def _typed_amine_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amine":
        return kind
    # 环胺：kind 收敛为 FG 类别 amine；cyclo 前缀 / 稠环词干由 chain_engine 动态加。
    if _scaffold(numbered) == "carbocycle" or _scaffold(numbered) in _RING_FG_SCAFFOLDS:
        return "amine"
    if kind not in {"amine", "diamine", "triamine", "tetraamine"}:
        return kind
    return "amine"  # 链式:数量统一由 multiplicity 承载,kind 恒为 amine

def _typed_expression_kind(kind: str, numbered: dict) -> str:
    kind = _typed_acid_kind(kind, numbered)
    kind = _typed_ester_kind(kind, numbered)
    kind = _typed_aldehyde_kind(kind, numbered)
    kind = _typed_nitrile_kind(kind, numbered)
    kind = _typed_amide_kind(kind, numbered)
    kind = _typed_ketone_kind(kind, numbered)
    kind = _typed_alcohol_kind(kind, numbered)
    return _typed_amine_kind(kind, numbered)
