from __future__ import annotations

def _typed_acid_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "acid":
        return kind
    if facts.relation.value == "exocyclic" and kind == "cycloalkane":
        return "cycloalkanecarboxylic" if facts.multiplicity == 1 else "cycloalkane_polycarboxylic"
    if facts.relation.value == "exocyclic":
        return kind
    return "acid" if facts.multiplicity == 1 else "diacid" if facts.multiplicity == 2 else "polycarboxylic"

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
    ring_kind = _typed_ring_alcohol_kind(kind, numbered, facts)
    if ring_kind:
        return ring_kind
    if kind not in {"alcohol", "diol", "triol"}:
        return kind
    return "alcohol" if facts.multiplicity == 1 else "diol" if facts.multiplicity == 2 else "triol"

def _typed_amine_kind(kind: str, numbered: dict) -> str:
    facts = (numbered.get("parent") or {}).get("principal_expression_facts")
    if not facts or facts.group_class.value != "amine":
        return kind
    if kind in {"cycloalkane", "cycloamine"} and facts.multiplicity == 1:
        return "cycloamine"
    if kind not in {"amine", "diamine", "triamine", "tetraamine"}:
        return kind
    return {1: "amine", 2: "diamine", 3: "triamine", 4: "tetraamine"}.get(facts.multiplicity, kind)

def _typed_expression_kind(kind: str, numbered: dict) -> str:
    kind = _typed_acid_kind(kind, numbered)
    kind = _typed_ketone_kind(kind, numbered)
    kind = _typed_alcohol_kind(kind, numbered)
    return _typed_amine_kind(kind, numbered)
