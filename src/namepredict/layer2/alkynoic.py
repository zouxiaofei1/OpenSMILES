"""Open-chain mono alkynoic acid / alkynoate parents (IUPAC P-65.1.1 / P-31.1).

Mono FG + mono C≡C (no C=C); FG carbon and triple-bond ends acyclic.
Parent keeps kind=acid/ester and carries triple_bond for L4/L5.
"""
from __future__ import annotations


def _open_chain_ynsat_atoms(mol, fg_c: int, tb: dict) -> bool:
    """True iff principal FG attach carbon and both C≡C ends are acyclic."""
    atoms = (int(fg_c), int(tb["c1"]), int(tb["c2"]))
    return all(not mol.GetAtomWithIdx(a).IsInRing() for a in atoms)


def _ok_ynsat_fg(info: dict, flag: str, ekey: str, bad: tuple) -> bool:
    """Mono open-chain yne FG: FG carbon + C≡C ends not in ring."""
    from namepredict.layer2.parent_core import _is_mono_alkyne, _is_mono_fg, _no_fgs
    if not _is_mono_fg(info, flag, ekey) or not _is_mono_alkyne(info):
        return False
    c, tb, mol = info[ekey][0]["c_idx"], info["triple_bonds"][0], info["mol"]
    return _open_chain_ynsat_atoms(mol, c, tb) and _no_fgs(info, bad)


def _ynsat_cover_atoms(info, c_idx: int, tb: dict) -> list[int]:
    return [c_idx, tb["c1"], tb["c2"]]


def _try_ynsat_fg(info, flag, ekey, bad, kind, ckey, **extra) -> dict | None:
    from namepredict.layer2.parent_core import _best_cover_pair, _parent_dict
    if not _ok_ynsat_fg(info, flag, ekey, bad):
        return None
    c_idx, tb = info[ekey][0]["c_idx"], info["triple_bonds"][0]
    chain = _best_cover_pair(info["mol"], _ynsat_cover_atoms(info, c_idx, tb))
    if not chain or c_idx not in chain:
        return None
    kw = {ckey: c_idx, "triple_bond": (tb["c1"], tb["c2"]), "mol": info["mol"], **extra}
    return _parent_dict(chain, kind, **kw)


def ynsat_or_unsat_or_sat(info, flag, ekey, bad, ukind, skind, ckey, **extra):
    """Try mono-yne FG parent, then multi-ene, then mono-ene, then saturated."""
    from namepredict.layer2.parent_core import (
        _fg_chain, _try_polyunsat_fg, _try_unsat_fg,
    )
    for fn in (_try_ynsat_fg, _try_polyunsat_fg, _try_unsat_fg):
        p = fn(info, flag, ekey, bad, ukind, ckey, **extra)
        if p is not None:
            return p
    return _fg_chain(info, ekey, skind, ckey, **extra)
