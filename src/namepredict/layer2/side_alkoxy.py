"""Ring outer linear n-alkoxy C1–C4 and PEG tails -(OCH2CH2)k-OR (L2 topology)."""
from __future__ import annotations

from rdkit.Chem import Mol


def _heavies(atom) -> list:
    return [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]


def _open_c_ok(mol: Mol, idx: int) -> bool:
    a = mol.GetAtomWithIdx(idx)
    return a.GetAtomicNum() == 6 and not a.IsInRing() and not a.GetIsAromatic()


def _fwd_heavies(mol: Mol, cur: int, prev: int) -> list:
    return [x for x in _heavies(mol.GetAtomWithIdx(cur)) if x.GetIdx() != prev]


def _alkyl_c_ok(mol: Mol, idx: int, prev: int) -> bool:
    """Open C: non-H neighbors are only C or prev (ether O allowed as prev)."""
    if not _open_c_ok(mol, idx):
        return False
    for n in mol.GetAtomWithIdx(idx).GetNeighbors():
        if n.GetAtomicNum() not in (1, 6) and n.GetIdx() != prev:
            return False
    return True


def _one_fwd_c(mol: Mol, cur: int, prev: int) -> int | None:
    fwd = [x for x in _fwd_heavies(mol, cur, prev) if x.GetAtomicNum() == 6]
    return fwd[0].GetIdx() if len(fwd) == 1 else None


def _no_extra_heavy(mol: Mol, cur: int, prev: int) -> bool:
    return not [x for x in _fwd_heavies(mol, cur, prev) if x.GetAtomicNum() != 1]


def _try_append_c(mol: Mol, path: list[int], back: int) -> list[int] | None:
    nxt = _one_fwd_c(mol, path[-1], back)
    if nxt is None:
        return path if _no_extra_heavy(mol, path[-1], back) else None
    if not _alkyl_c_ok(mol, nxt, path[-1]) or len(path) >= 4:
        return None
    return path + [nxt]


def _linear_alkyl_atoms(mol: Mol, start: int, prev: int) -> list[int] | None:
    """Linear C1–C4 chain from start; prev may be ether O."""
    if not _alkyl_c_ok(mol, start, prev):
        return None
    path, back = [start], prev
    for _ in range(4):
        out = _try_append_c(mol, path, back)
        if out is None or len(out) == len(path):
            return out
        back, path = path[-1], out
    return None


def _peg_o_after(mol: Mol, c2: int, c1: int) -> int | None:
    if not _open_c_ok(mol, c2):
        return None
    fwd = [x for x in _fwd_heavies(mol, c2, c1) if x.GetAtomicNum() == 8]
    if len(fwd) != 1 or mol.GetAtomWithIdx(fwd[0].GetIdx()).IsInRing():
        return None
    return fwd[0].GetIdx()


def _peg_unit(mol: Mol, c1: int, prev: int) -> tuple[int, list[int]] | None:
    """One -CH2-CH2-O-: return (O_idx, [c1, c2, o])."""
    if not _open_c_ok(mol, c1):
        return None
    c2 = _one_fwd_c(mol, c1, prev)
    if c2 is None:
        return None
    o2 = _peg_o_after(mol, c2, c1)
    return (o2, [c1, c2, o2]) if o2 is not None else None


def _term_after_o(mol: Mol, o_idx: int, prev_c: int) -> tuple[int, list[int]] | None:
    """Me/Et after ether O (PEG terminal)."""
    nbs = [x for x in _fwd_heavies(mol, o_idx, prev_c) if x.GetAtomicNum() == 6]
    if len(nbs) != 1:
        return None
    path = _linear_alkyl_atoms(mol, nbs[0].GetIdx(), o_idx)
    return (len(path), path) if path is not None and len(path) in (1, 2) else None


def _fix_peg_code(peg_k: int, n_term: int) -> int:
    if peg_k == 0:
        return n_term
    if n_term not in (1, 2) or not (1 <= peg_k <= 3):
        return 0
    return int("1" * (peg_k - 1) + str(n_term) + "2")


def _peg_result(k: int, atoms: list[int], term: tuple[int, list[int]]) -> dict | None:
    n_term, t_atoms = term
    code = _fix_peg_code(k, n_term)
    return None if not code else {
        "n_terminal": n_term, "peg_k": k, "atoms": atoms + t_atoms, "code": code,
    }


def _next_peg_c(mol: Mol, o_idx: int, prev_c: int) -> int | None:
    nbs = [x for x in _fwd_heavies(mol, o_idx, prev_c) if x.GetAtomicNum() == 6]
    return nbs[0].GetIdx() if len(nbs) == 1 else None


def _peg_after_unit(
    mol: Mol, o_next: int, unit: list[int], atoms: list[int], k: int,
) -> tuple[dict | None, tuple | None]:
    """After one PEG unit: (done, cont_state)."""
    atoms = atoms + unit
    term = _term_after_o(mol, o_next, unit[1])
    if term is not None:
        return _peg_result(k, atoms, term), None
    nxt = _next_peg_c(mol, o_next, unit[1])
    return (None, None) if nxt is None else (None, (nxt, o_next, atoms))


def _apply_peg_unit(
    mol: Mol, cur_c: int, back: int, atoms: list[int], k: int,
) -> tuple[dict | None, tuple | None]:
    step = _peg_unit(mol, cur_c, back)
    if step is None:
        return None, None
    return _peg_after_unit(mol, step[0], step[1], atoms, k)


def _walk_peg_from(mol: Mol, c1: int, prev: int) -> dict | None:
    """O_ring-(CH2CH2O)k-R, k=1..3, R=Me/Et."""
    atoms: list[int] = []
    cur_c, back = c1, prev
    for k in range(1, 4):
        done, cont = _apply_peg_unit(mol, cur_c, back, atoms, k)
        if done is not None or cont is None:
            return done
        cur_c, back, atoms = cont
    return None


def _parse_simple_alkoxy(mol: Mol, start: int, o_idx: int) -> dict | None:
    path = _linear_alkyl_atoms(mol, start, o_idx)
    if path is None:
        return None
    n = len(path)
    return {"n_terminal": n, "peg_k": 0, "atoms": path, "code": n}


def _parse_outer_alkoxy(mol: Mol, start: int, o_idx: int) -> dict | None:
    """Linear n-alkoxy C1–C4 or PEG -(OCH2CH2)k-OR (k≤3, R=Me/Et)."""
    if not _open_c_ok(mol, start):
        return None
    return _parse_simple_alkoxy(mol, start, o_idx) or _walk_peg_from(mol, start, o_idx)


def _outer_alkoxy_n(mol: Mol, start: int, o_idx: int) -> int:
    """Compat: alkoxy code or 0 if unsupported."""
    p = _parse_outer_alkoxy(mol, start, o_idx)
    return 0 if p is None else int(p["code"])


def _outer_atoms(mol: Mol, start: int, o_idx: int, n: int) -> list[int]:
    p = _parse_outer_alkoxy(mol, start, o_idx)
    if p is not None and int(p["code"]) == n:
        return list(p["atoms"])
    return [start]
