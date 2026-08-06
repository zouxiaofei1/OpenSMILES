"""Arene retained parents: benzoic, benzaldehyde, acetophenone, benzoate, etc."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.candidate_gate import CandidateGate, GateScope, GateStatus, pass_gate, scoped_reject

from namepredict.layer2.aryl_sub import (
    _aryl_atoms,
    _aryl_exclude,
    _aryl_sub_n,
    _exocyclic_fg_ring,
    _ring_benzyloxys,
    _ring_phenoxys,
)
from namepredict.layer2.scaffold.ring_parent import (
    _arene_alkoxy,
    _arene_fg_subs_ok,
    _dbl_o_idx,
    _is_benzene_core,
    _phenol_amines_ok,
    _ring_alkoxy_ethers,
    _ring_nitro_atoms,
    _ring_nitro_n,
    _ring_primary_amines,
)


def _ring_c_neighbors(mol: Mol, c_idx: int, ring_set: set[int]) -> list[int]:
    return [
        n.GetIdx()
        for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() in ring_set
    ]


def _fg_ring_c(info: dict, ring_set: set[int], ekey: str) -> int | None:
    entries = info.get(ekey) or []
    if len(entries) != 1:
        return None
    nbs = _ring_c_neighbors(info["mol"], entries[0]["c_idx"], ring_set)
    return nbs[0] if len(nbs) == 1 else None


def _carboxyl_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "carboxyls")


def _aldehyde_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "aldehydes")


def _ketone_ring_c(info: dict, ring_set: set[int]) -> int | None:
    return _fg_ring_c(info, ring_set, "ketones")


def _cooh_oxygen_idxs(mol: Mol, c_idx: int) -> set[int]:
    atom = mol.GetAtomWithIdx(c_idx)
    return {n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 8}


def _ring_phenol_ohs(info: dict, ring_set: set[int]) -> list[dict]:
    return [h for h in (info.get("hydroxyls") or []) if h["c_idx"] in ring_set]


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


def _arene_ethers_ok(info: dict, ring_set: set[int]) -> bool:
    n_ether = len(info.get("ethers") or [])
    n_ok = (
        len(_ring_alkoxy_ethers(info, ring_set))
        + len(_ring_phenoxys(info, ring_set))
        + len(_ring_benzyloxys(info, ring_set))
    )
    return n_ether == n_ok


def _arene_ring_prefix_ok(info: dict, ring_set: set[int]) -> bool:
    if _phenol_amines_ok(info, ring_set) is None:
        return False
    return _arene_ethers_ok(info, ring_set)


def _arene_prefix_atoms(info: dict, ring_set: set[int]) -> set[int]:
    ams = _ring_primary_amines(info, ring_set)
    ohs = _ring_phenol_ohs(info, ring_set)
    alk, _ = _arene_alkoxy(info, ring_set)
    return (
        {a["n_idx"] for a in ams}
        | {h["o_idx"] for h in ohs}
        | _ring_nitro_atoms(info, ring_set)
        | alk
        | _aryl_atoms(info, ring_set)
    )


def _arene_sub_ctx(info: dict, ring_set: set[int], exclude: set[int], allowed: set[int]):
    alk, n_alk = _arene_alkoxy(info, ring_set)
    ams = _ring_primary_amines(info, ring_set)
    n_oh = len(_ring_phenol_ohs(info, ring_set))
    full = allowed | _arene_prefix_atoms(info, ring_set)
    excl = exclude | alk | _aryl_exclude(info, ring_set)
    return full, excl, n_alk, len(ams) + n_oh


def _arene_subs_ok(
    info: dict, mol: Mol, ring_set: set[int], exclude: set[int], allowed: set[int],
) -> bool:
    if not _arene_ring_prefix_ok(info, ring_set):
        return False
    full, excl, n_alk, n_am = _arene_sub_ctx(info, ring_set, exclude, allowed)
    return _arene_fg_subs_ok(
        mol, ring_set, full, _ring_nitro_n(info, ring_set), n_am, n_alk, excl,
        _aryl_sub_n(info, ring_set),
    )


def _pick_fg_ring(info: dict, ekey: str) -> set[int] | None:
    entries = info.get(ekey) or []
    if len(entries) != 1:
        return None
    return _exocyclic_fg_ring(info["mol"], entries[0]["c_idx"])


def _is_simple_benzoic(info: dict) -> bool:
    if not _is_benzene_core(info) or _arene_fg_conflict(info, "has_aldehyde", "has_ketone"):
        return False
    ring = _pick_fg_ring(info, "carboxyls")
    if ring is None or _carboxyl_ring_c(info, ring) is None:
        return False
    mol, fg_c = info["mol"], info["carboxyls"][0]["c_idx"]
    return _arene_subs_ok(info, mol, ring, {fg_c}, _cooh_oxygen_idxs(mol, fg_c))


def _benzoic_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "carboxyls") or set()
    cooh_c = info["carboxyls"][0]["c_idx"]
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "benzoic",
        "cooh_c_idx": cooh_c, "ring_attach_idx": _carboxyl_ring_c(info, ring),
    }


def _try_benzoic_parent(info: dict) -> dict | None:
    return _benzoic_parent(info) if _is_simple_benzoic(info) else None


def _oxo_fg_ok(info: dict, ekey: str, *conflict: str) -> tuple | None:
    if not _is_benzene_core(info) or _arene_fg_conflict(info, *conflict):
        return None
    ring = _pick_fg_ring(info, ekey)
    if ring is None or _fg_ring_c(info, ring, ekey) is None:
        return None
    mol, fg_c = info["mol"], info[ekey][0]["c_idx"]
    o_idx = _dbl_o_idx(mol, fg_c)
    return (mol, ring, fg_c, o_idx) if o_idx is not None else None


def _is_simple_benzaldehyde(info: dict) -> bool:
    got = _oxo_fg_ok(info, "aldehydes", "has_acid", "has_ketone")
    if got is None:
        return False
    mol, ring, fg_c, o_idx = got
    return _arene_subs_ok(info, mol, ring, {fg_c}, {o_idx})


def _benzaldehyde_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "aldehydes") or set()
    ald_c = info["aldehydes"][0]["c_idx"]
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "benzaldehyde",
        "aldehyde_c_idx": ald_c, "ring_attach_idx": _aldehyde_ring_c(info, ring),
    }


def _try_benzaldehyde_parent(info: dict) -> dict | None:
    return _benzaldehyde_parent(info) if _is_simple_benzaldehyde(info) else None


def _is_methyl_carbon(mol: Mol, c_idx: int, only_nb: int) -> bool:
    atom = mol.GetAtomWithIdx(c_idx)
    if atom.GetAtomicNum() != 6 or atom.IsInRing():
        return False
    heavies = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == only_nb


def _acetyl_methyl_c(mol: Mol, ket_c: int, ring_set: set[int]) -> int | None:
    atom = mol.GetAtomWithIdx(ket_c)
    cands = [
        n.GetIdx() for n in atom.GetNeighbors()
        if n.GetAtomicNum() == 6 and n.GetIdx() not in ring_set
    ]
    if len(cands) != 1:
        return None
    me = cands[0]
    return me if _is_methyl_carbon(mol, me, ket_c) else None


def _is_simple_acetophenone(info: dict) -> bool:
    got = _oxo_fg_ok(info, "ketones", "has_acid", "has_aldehyde")
    if got is None:
        return False
    mol, ring, ket_c, o_idx = got
    me = _acetyl_methyl_c(mol, ket_c, ring)
    if me is None:
        return False
    return _arene_subs_ok(info, mol, ring, {ket_c, me}, {o_idx})


def _acetophenone_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "ketones") or set()
    ket_c = info["ketones"][0]["c_idx"]
    me = _acetyl_methyl_c(info["mol"], ket_c, ring)
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "acetophenone",
        "ketone_c_idx": ket_c, "acetyl_methyl_idx": me,
        "ring_attach_idx": _ketone_ring_c(info, ring),
    }


def _try_acetophenone_parent(info: dict) -> dict | None:
    return _acetophenone_parent(info) if _is_simple_acetophenone(info) else None


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


def _ester_alkoxy_arm(mol: Mol, start: int, o_idx: int) -> set[int]:
    arm, cur, prev = set(), start, o_idx
    while cur is not None and len(arm) < 17:
        if not _is_alkoxy_c(mol, cur, prev):
            break
        arm.add(cur)
        nxt = _alkoxy_next(mol, cur, prev)
        if nxt is None or nxt == -1:
            break
        prev, cur = cur, nxt
    return arm


def _bfs_side(mol: Mol, start: int, ban: set[int]) -> set[int]:
    """Heavy-atom component from start, not crossing ban."""
    from collections import deque
    seen: set[int] = set()
    q: deque[int] = deque([start])
    while q:
        i = q.popleft()
        if i in seen or i in ban:
            continue
        seen.add(i)
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetAtomicNum() != 1:
                q.append(n.GetIdx())
    return seen


def _ester_o_side(mol: Mol, e: dict) -> set[int]:
    """Full alcohol-side component beyond ester O (linear or polycyclic)."""
    linear = _ester_alkoxy_arm(mol, e["alkoxy_c_idx"], e["o_idx"])
    if linear and _simple_alkoxy_n(mol, e["alkoxy_c_idx"], e["o_idx"]) is not None:
        return linear
    return _bfs_side(mol, e["alkoxy_c_idx"], {e["o_idx"]})


def _ester_exclude(mol: Mol, e: dict) -> tuple[set[int], set[int]]:
    fg_c, o_idx, ac = e["c_idx"], e["o_idx"], e["alkoxy_c_idx"]
    o_dbl = _dbl_o_idx(mol, fg_c)
    allowed = {o_idx, ac} | ({o_dbl} if o_dbl is not None else set())
    exclude = {fg_c, o_idx} | _ester_o_side(mol, e) | ({o_dbl} if o_dbl is not None else set())
    return exclude, allowed | exclude


def _arene_fg_ctx(
    info: dict, allow_flag: str, ekey: str,
) -> tuple[Mol, set[int]] | None:
    allow = frozenset({allow_flag})
    if not _is_benzene_core(info) or _arene_fg_conflict(
        info, "has_acid", "has_aldehyde", "has_ketone", allow=allow,
    ):
        return None
    ring = _pick_fg_ring(info, ekey)
    if ring is None or _fg_ring_c(info, ring, ekey) is None:
        return None
    return info["mol"], ring


def _benzoate_alkoxy(mol: Mol, o_idx: int, ac: int) -> dict | None:
    """Linear C1–C16 or P-65.6 special alkoxy (benzyl/tBu/iPr/Ph)."""
    from namepredict.layer2.alkoxy_side import classify_alkoxy
    side = classify_alkoxy(mol, o_idx, ac)
    if side.get("alkoxy_en"):
        return side
    n = _simple_alkoxy_n(mol, ac, o_idx)
    return {**side, "alkoxy_n": n} if n is not None else None


def _ph_ring_plain(mol: Mol, ring: set[int], excl: set[int]) -> bool:
    """Benzene ring has no unclaimed outside heavy attachments."""
    from namepredict.layer2.scaffold.ring_parent import _outside_ok
    return _outside_ok(mol, ring, excl)


def _complex_benzoate_ok(info: dict) -> bool:
    """Ph–C(=O)–O–R: clean Ph, any alcohol-side R (not simple linear/special)."""
    ctx = _arene_fg_ctx(info, "has_ester", "esters")
    if ctx is None:
        return False
    mol, ring = ctx
    e = info["esters"][0]
    if _benzoate_alkoxy(mol, e["o_idx"], e["alkoxy_c_idx"]) is not None:
        return False  # simple path owns these
    excl, allowed = _ester_exclude(mol, e)
    return _ph_ring_plain(mol, ring, allowed | excl)


def _is_simple_benzoate(info: dict) -> bool:
    ctx = _arene_fg_ctx(info, "has_ester", "esters")
    if ctx is None:
        return False
    mol, ring = ctx
    e = info["esters"][0]
    if _benzoate_alkoxy(mol, e["o_idx"], e["alkoxy_c_idx"]) is None:
        return False
    excl, allowed = _ester_exclude(mol, e)
    return _arene_subs_ok(info, mol, ring, excl, allowed)


def _benzoate_side_fields(side: dict) -> dict:
    return {
        "alkoxy_n": side.get("alkoxy_n"),
        "alkoxy_en": side.get("alkoxy_en") or "",
        "alkoxy_zh": side.get("alkoxy_zh") or "",
    }


def _benzoate_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "esters") or set()
    e, mol = info["esters"][0], info["mol"]
    side = _benzoate_alkoxy(mol, e["o_idx"], e["alkoxy_c_idx"]) or {}
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "benzoate",
        "ester_c_idx": e["c_idx"], "o_idx": e["o_idx"],
        "alkoxy_c_idx": e["alkoxy_c_idx"], **_benzoate_side_fields(side),
        "ring_attach_idx": _fg_ring_c(info, ring, "esters"),
    }


def _complex_benzoate_parent(info: dict) -> dict:
    """Benzoate with unresolved O-alkyl; L5 emits bare benzoate until named."""
    p = _benzoate_parent(info)
    p["alkoxy_complex"] = True
    return p


def _try_benzoate_parent(info: dict) -> dict | None:
    if _is_simple_benzoate(info):
        return _benzoate_parent(info)
    return _complex_benzoate_parent(info) if _complex_benzoate_ok(info) else None


def _nitrile_n_idx(mol: Mol, c_idx: int) -> int | None:
    atom = mol.GetAtomWithIdx(c_idx)
    for n in atom.GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _is_simple_benzonitrile(info: dict) -> bool:
    ctx = _arene_fg_ctx(info, "has_nitrile", "nitriles")
    if ctx is None:
        return False
    mol, ring = ctx
    fg_c = info["nitriles"][0]["c_idx"]
    n_idx = _nitrile_n_idx(mol, fg_c)
    return n_idx is not None and _arene_subs_ok(info, mol, ring, {fg_c}, {n_idx})


def _benzonitrile_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "nitriles") or set()
    c = info["nitriles"][0]["c_idx"]
    return {
        "chain": list(ring), "n_carbons": 6, "kind": "benzonitrile",
        "nitrile_c_idx": c, "ring_attach_idx": _fg_ring_c(info, ring, "nitriles"),
    }


def _try_benzonitrile_parent(info: dict) -> dict | None:
    return _benzonitrile_parent(info) if _is_simple_benzonitrile(info) else None


def _acyl_allowed(mol: Mol, e: dict) -> set[int]:
    o_dbl = _dbl_o_idx(mol, e["c_idx"])
    allowed = {e.get("hal_idx", e["cl_idx"])}
    return allowed | ({o_dbl} if o_dbl is not None else set())


def _is_simple_benzoyl_halide(info: dict) -> bool:
    ctx = _arene_fg_ctx(info, "has_acyl_chloride", "acyl_chlorides")
    if ctx is None:
        return False
    mol, ring = ctx
    e = info["acyl_chlorides"][0]
    return _arene_subs_ok(info, mol, ring, {e["c_idx"]}, _acyl_allowed(mol, e))


def _benzoyl_hal_kind(e: dict) -> str:
    return "benzoyl_bromide" if int(e.get("hal_z") or 17) == 35 else "benzoyl_chloride"


def _benzoyl_chloride_parent(info: dict) -> dict:
    ring = _pick_fg_ring(info, "acyl_chlorides") or set()
    e = info["acyl_chlorides"][0]
    return {
        "chain": list(ring), "n_carbons": 6, "kind": _benzoyl_hal_kind(e),
        "acyl_c_idx": e["c_idx"], "cl_idx": e["cl_idx"],
        "hal_idx": e.get("hal_idx", e["cl_idx"]), "hal_z": e.get("hal_z", 17),
        "ring_attach_idx": _fg_ring_c(info, ring, "acyl_chlorides"),
    }


def _try_benzoyl_chloride_parent(info: dict) -> dict | None:
    return _benzoyl_chloride_parent(info) if _is_simple_benzoyl_halide(info) else None


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


def _benzene_polycarboxylic_parent(info: dict) -> dict:
    ring = _benzene_polyacid_ring(info) or set()
    attaches = [_benzene_polyacid_attach(info, ring, acid) for acid in info["carboxyls"]]
    return {"chain": list(ring), "n_carbons": 6, "kind": "benzene_polycarboxylic",
            "cooh_c_idxs": attaches, "acid_count": len(attaches)}


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


def _try_benzene_polycarboxylic(info: dict) -> dict | None:
    return _benzene_polycarboxylic_parent(info) if _is_benzene_polycarboxylic(info) else None


def _try_benzamide_parent(info: dict) -> dict | None:
    """Retained Ph–C(=O)–N parent (implementation in benzamide module)."""
    from namepredict.layer2.benzamide import _try_benzamide_parent as _try
    return _try(info)


def _try_arene_other_fg(info: dict) -> dict | None:
    for fn in (
        _try_benzoyl_chloride_parent, _try_benzoate_parent, _try_benzonitrile_parent,
        _try_benzamide_parent,
    ):
        if (p := fn(info)) is not None:
            return p
    return None
