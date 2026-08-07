"""Arene retained parents: benzoic, benzaldehyde, acetophenone, benzoate, etc."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.candidate_gate import CandidateGate, GateScope, GateStatus, pass_gate, scoped_reject

from namepredict.layer3.aryl_sub import (
    _exocyclic_fg_ring,
)
from namepredict.layer2.scaffold.ring_parent import (
    _is_benzene_core,
)


def _ring_c_neighbors(mol: Mol, c_idx: int, ring_set: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() in ring_set
    ]


def _cooh_oxygen_idxs(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8}


def _arene_fg_conflict(
    info: dict, *extra: str, allow: frozenset[str] | set[str] = frozenset(),
) -> bool:
    bad = (
        "has_ester", "has_amide", "has_acyl_chloride", "has_anhydride",
        "has_nitrile", *extra,
    )
    if any(info.get(k) for k in bad if k not in allow):
        return True
    return bool(info.get("thiols"))


def _is_alkoxy_c(mol: Mol, cur: int, prev: int) -> bool:
    atom = mol.GetAtomWithIdx(cur)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    for n in atom.GetNeighbors():
        z = n.GetAtomicNum()
        if z == 1 or n.GetIdx() == prev:
            continue
        if z != 6:
            return False
    return True


def _alkoxy_next(mol: Mol, cur: int, prev: int) -> int | None:
    free = [
        n.GetIdx() for n in mol.GetAtomWithIdx(cur).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() != prev
    ]
    return free[0] if len(free) == 1 else None if not free else -1


def _simple_alkoxy_n(mol: Mol, start: int, o_idx: int) -> int | None:
    """Count linear n-alkyl ester alcohol side (C1–C16); reject branch/hetero."""
    n, cur, prev = 0, start, o_idx
    while cur is not None and n < 17:
        if not _is_alkoxy_c(mol, cur, prev) or (nxt := _alkoxy_next(mol, cur, prev)) == -1:
            return None
        n, prev, cur = n + 1, cur, nxt
    return n if 1 <= n <= 16 else None


def _benzoate_alkoxy(mol: Mol, o_idx: int, ac: int) -> dict | None:
    """Linear C1–C16 or P-65.6 special alkoxy (benzyl/tBu/iPr/Ph)."""
    from namepredict.tools.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(mol, o_idx, ac)
    if side.get("alkoxy_en"):
        return side
    n = _simple_alkoxy_n(mol, ac, o_idx)
    return {**side, "alkoxy_n": n} if n is not None else None


def _benzoate_side_fields(side: dict) -> dict:
    return {
        "alkoxy_n": side.get("alkoxy_n"),
        "alkoxy_en": side.get("alkoxy_en") or "",
        "alkoxy_zh": side.get("alkoxy_zh") or "",
    }


def _benzene_polyacid_attach(info: dict, ring: set[int], acid: dict) -> int | None:
    nbs = _ring_c_neighbors(info["mol"], acid["c_idx"], ring)
    return nbs[0] if len(nbs) == 1 else None


def _benzene_polyacid_ring(info: dict) -> set[int] | None:
    acids = info.get("carboxyls") or []
    if len(acids) not in (2, 3) or not _is_benzene_core(info):
        return None
    ring = _exocyclic_fg_ring(info["mol"], acids[0]["c_idx"])
    return ring if ring and all(_benzene_polyacid_attach(info, ring, acid) is not None for acid in acids) else None


def _is_benzene_polycarboxylic(info: dict) -> bool:
    if _arene_fg_conflict(info, "has_aldehyde", "has_ketone") or any(
        acid.get("anion") for acid in info["carboxyls"]
    ):
        return False
    return _benzene_polyacid_ring(info) is not None


def _benzene_polyacid_like(info: dict) -> bool:
    acids = info.get("carboxyls") or []
    derivatives = sum(len(info.get(key) or []) for key in ("esters", "amides", "acyl_chlorides"))
    return len(acids) >= 2 or bool(acids and derivatives)


def _has_benzene_polyacid_attachment(info: dict) -> bool:
    for acid in info.get("carboxyls") or []:
        ring = _exocyclic_fg_ring(info["mol"], acid["c_idx"])
        if ring and _benzene_polyacid_attach(info, ring, acid) is not None:
            return True
    return False


def benzene_polycarboxylic_gate(info: dict) -> CandidateGate:
    """Gate only direct multiacid attachments on an unfused benzene core."""
    scope = GateScope.BENZENE_POLYCARBOXYLIC
    if not _is_benzene_core(info) or not _benzene_polyacid_like(info):
        return pass_gate(scope)
    if not _has_benzene_polyacid_attachment(info):
        return pass_gate(scope)
    if _is_benzene_polycarboxylic(info):
        return pass_gate(scope)
    return scoped_reject(scope, "unsupported_benzene_polyacid")


def benzene_polycarboxylic_eligibility(info: dict) -> bool | None:
    """Compatibility projection of the typed benzene gate."""
    gate = benzene_polycarboxylic_gate(info)
    in_scope = _is_benzene_core(info) and _benzene_polyacid_like(info)
    return None if not in_scope or not _has_benzene_polyacid_attachment(info) else gate.status is GateStatus.PASS
