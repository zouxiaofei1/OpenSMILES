from __future__ import annotations

# 苯环 + 单主官能团 → P-22.1.3 保留名（原 L2 _RETAINED_RING_KINDS 迁此）。
# L2 只产生结构 kind='benzene'；保留名决策完全在 L5。
_BENZENE_RETAINED = {
    "acid": "benzoic",
    "ester": "benzoate",
    "aldehyde": "benzaldehyde",
    "nitrile": "benzonitrile",
    "amide": "benzamide",
    "amine": "aniline",
    "alcohol": "phenol",
}


def _benzene_retained(kind: str, numbered: dict) -> str | None:
    """Benzene ring + exactly one principal FG → retained ring kind, else None."""
    if kind != "benzene":
        return None
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.multiplicity != 1:
        return None
    return _BENZENE_RETAINED.get(facts.group_class.value)


def _typed_acid_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "acid":
        return kind
    retained = _benzene_retained(kind, numbered)
    if retained:
        return retained
    if facts.relation.value == "exocyclic" and kind == "cycloalkane":
        return "cycloalkanecarboxylic" if facts.multiplicity == 1 else "cycloalkane_polycarboxylic"
    if facts.relation.value == "exocyclic":
        return kind
    return "acid"  # 链式:数量统一由 multiplicity 承载,kind 恒为 acid


def _typed_ester_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ester":
        return kind
    return _benzene_retained(kind, numbered) or kind


def _typed_aldehyde_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "aldehyde":
        return kind
    return _benzene_retained(kind, numbered) or kind


def _typed_nitrile_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "nitrile":
        return kind
    return _benzene_retained(kind, numbered) or kind


def _typed_amide_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amide":
        return kind
    return _benzene_retained(kind, numbered) or kind

def _typed_ketone_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "ketone":
        return kind
    if kind in {"cycloalkane", "cycloketone", "cycloalkanedione"}:
        return "cycloketone" if facts.multiplicity == 1 else "cycloalkanedione"
    if kind not in {"ketone", "dione"}:
        return kind
    return "ketone" if facts.multiplicity == 1 else "dione"

def _typed_ring_alcohol_kind(kind: str, numbered: dict, facts) -> str | None:
    if kind in {"cycloalkane", "cycloalcohol", "cycloalkanediol"}:
        return "cycloalcohol" if facts.multiplicity == 1 else "cycloalkanediol"
    if kind in {"benzene", "phenol", "benzenediol"}:
        return "phenol" if facts.multiplicity == 1 else "benzenediol"
    parent = numbered.get("parent") or {}
    scaffold = parent.get("scaffold_identity")
    if parent.get("typed_ring_expression_supported") and scaffold and scaffold.id == "naphthalene":
        return "naphthalenol" if facts.multiplicity == 1 else "naphthalenediol"
    return None

def _typed_alcohol_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "alcohol":
        return kind
    retained = _benzene_retained(kind, numbered)
    if retained:
        return retained
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
    retained = _benzene_retained(kind, numbered)
    if retained:
        return retained
    if kind in {"cycloalkane", "cycloamine"} and facts.multiplicity == 1:
        return "cycloamine"
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
