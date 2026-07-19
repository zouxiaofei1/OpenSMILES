"""Terminal parent ownership: immutable owned_atoms finalization."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.ring_parent import _dbl_o_idx


def _chain_atoms(parent: dict) -> set[int]:
    return set(parent.get("chain") or [])


def _add_opt(out: set[int], idx: int | None) -> set[int]:
    if idx is not None:
        out.add(idx)
    return out


def _amide_n_from_c(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for n in carbon.GetNeighbors():
        if n.GetAtomicNum() == 7:
            return n.GetIdx()
    return None


def _amide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    c_idx = parent.get("amide_c_idx")
    if c_idx is None:
        return set()
    out = {c_idx}
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _amide_n_from_c(mol, c_idx))


def _single_o_idx(mol: Mol, c_idx: int) -> int | None:
    carbon = mol.GetAtomWithIdx(c_idx)
    for bond in carbon.GetBonds():
        if bond.GetBondType().name != "SINGLE":
            continue
        other = bond.GetOtherAtom(carbon)
        if other.GetAtomicNum() == 8:
            return other.GetIdx()
    return None


def _acid_o_atoms(mol: Mol, c_idx: int) -> set[int]:
    out: set[int] = set()
    _add_opt(out, _dbl_o_idx(mol, c_idx))
    return _add_opt(out, _single_o_idx(mol, c_idx))


def _cooh_c_idxs(parent: dict) -> list[int]:
    multi = parent.get("cooh_c_idxs")
    if multi:
        return [int(x) for x in multi]
    one = parent.get("cooh_c_idx")
    return [int(one)] if one is not None else []


def _acid_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    out: set[int] = set()
    for c in _cooh_c_idxs(parent):
        out.add(c)
        out |= _acid_o_atoms(mol, c)
    return out


def _aldehyde_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    c_idx = parent.get("aldehyde_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    return _add_opt(out, _dbl_o_idx(mol, int(c_idx)))


def _ether_arm_atoms(mol: Mol, o_idx: int, c_idx: int) -> set[int]:
    from namepredict.layer2.chain_walk import _longest_from
    return {c_idx, *_longest_from(mol, c_idx, {o_idx})}


def _ether_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Ether O + both carbon arms (short arm not always in chain)."""
    if parent.get("kind") != "ether" or parent.get("o_idx") is None:
        return set()
    o_idx = int(parent["o_idx"])
    out = {o_idx}
    for n in mol.GetAtomWithIdx(o_idx).GetNeighbors():
        if n.GetAtomicNum() == 6:
            out |= _ether_arm_atoms(mol, o_idx, n.GetIdx())
    return out


def _sulfide_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Sulfide S + both carbon arms (short arm not always in chain)."""
    if parent.get("kind") != "sulfide" or parent.get("s_idx") is None:
        return set()
    s_idx = int(parent["s_idx"])
    out = {s_idx}
    for n in mol.GetAtomWithIdx(s_idx).GetNeighbors():
        if n.GetAtomicNum() == 6:
            out |= _ether_arm_atoms(mol, s_idx, n.GetIdx())
    return out


def _sulfonic_acid_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Sulfonic acid S + all O neighbours."""
    s_idx = parent.get("s_idx")
    if s_idx is None:
        return set()
    out = {int(s_idx)}
    for n in mol.GetAtomWithIdx(int(s_idx)).GetNeighbors():
        if n.GetAtomicNum() == 8:
            out.add(n.GetIdx())
    return out


def _hydroxy_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Alcohol/diol/triol: attachment carbon(s) + OH oxygen(s)."""
    c_idxs = parent.get("oh_c_idxs") or ([parent.get("oh_c_idx")] if parent.get("oh_c_idx") is not None else [])
    if not c_idxs:
        return set()
    out: set[int] = set()
    for c in c_idxs:
        out.add(int(c))
        _add_opt(out, _single_o_idx(mol, int(c)))
    return out


def _ketone_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Ketone carbonyl C + double-bonded O (+ acetyl methyl when present)."""
    c_idx = parent.get("ketone_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    _add_opt(out, _dbl_o_idx(mol, int(c_idx)))
    return _add_opt(out, parent.get("acetyl_methyl_idx"))


def _amine_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Aniline/amine: ring/chain attachment C + amine N (primary)."""
    c_idx = parent.get("amine_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 7:
            out.add(n.GetIdx())
    return out


def _nitrile_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Nitrile/benzonitrile: CN carbon + triple-bonded N."""
    c_idx = parent.get("nitrile_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 7:
            out.add(n.GetIdx())
    return out


def _thiol_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Thiol: attachment carbon + SH sulfur."""
    c_idx = parent.get("sh_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    for n in mol.GetAtomWithIdx(int(c_idx)).GetNeighbors():
        if n.GetAtomicNum() == 16:
            out.add(n.GetIdx())
    return out


def _add_ester_alkoxy_arm(mol: Mol, c_idx: int, dbl_o: int | None, out: set[int]) -> None:
    """Find ester -O- neighbor of c_idx and walk its alkoxy arm."""
    for n in mol.GetAtomWithIdx(c_idx).GetNeighbors():
        if n.GetAtomicNum() != 8 or n.GetIdx() == dbl_o:
            continue
        out.add(n.GetIdx())
        for nn in n.GetNeighbors():
            if nn.GetIdx() != c_idx and nn.GetAtomicNum() == 6:
                out |= _ether_arm_atoms(mol, n.GetIdx(), nn.GetIdx())


def _one_ester_fg(mol: Mol, c_idx: int, out: set[int]) -> None:
    """Add one ester carbonyl C + =O + -O- + alkoxy arm to out."""
    out.add(int(c_idx))
    dbl_o = _dbl_o_idx(mol, int(c_idx))
    _add_opt(out, dbl_o)
    _add_ester_alkoxy_arm(mol, int(c_idx), dbl_o, out)


def _single_ester_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Mono-ester/benzoate: carbonyl C + =O + -O- + alkoxy arm."""
    c_idx = parent.get("ester_c_idx")
    if c_idx is None:
        return set()
    out: set[int] = set()
    _one_ester_fg(mol, c_idx, out)
    return out


def _diester_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Diester: both ester carbonyls + O atoms + alkoxy arms."""
    c_idxs = parent.get("ester_c_idxs")
    if not c_idxs:
        return set()
    out: set[int] = set()
    for c in c_idxs:
        _one_ester_fg(mol, c, out)
    return out


def _ester_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Ester/benzoate/diester: dispatch to mono or diester handler."""
    return _single_ester_fg_atoms(mol, parent) or _diester_fg_atoms(mol, parent)


def _anhydride_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Anhydride: both acyl carbons + bridge O + both carbonyl oxygens."""
    out: set[int] = set()
    for key in ("acyl_c_idx", "other_acyl_c_idx"):
        c = parent.get(key)
        if c is None:
            continue
        out.add(int(c))
        _add_opt(out, _dbl_o_idx(mol, int(c)))
    return _add_opt(out, parent.get("o_idx"))


def _acyl_chloride_fg_atoms(mol: Mol, parent: dict) -> set[int]:
    """Acyl chloride: acyl C + carbonyl O + Cl."""
    c_idx = parent.get("acyl_c_idx")
    if c_idx is None:
        return set()
    out = {int(c_idx)}
    _add_opt(out, _dbl_o_idx(mol, int(c_idx)))
    return _add_opt(out, parent.get("cl_idx") if parent.get("cl_idx") is not None else parent.get("hal_idx"))


def _kind_fg_atoms(parent: dict, mol: Mol) -> set[int]:
    """FG ownership is field-driven: own every heavy atom of the principal FG."""
    parts = (
        _amide_fg_atoms(mol, parent) if parent.get("amide_c_idx") is not None else set(),
        _aldehyde_fg_atoms(mol, parent) if parent.get("aldehyde_c_idx") is not None else set(),
        _acid_fg_atoms(mol, parent),
        _ether_fg_atoms(mol, parent),
        _sulfide_fg_atoms(mol, parent),
        _sulfonic_acid_fg_atoms(mol, parent),
        _hydroxy_fg_atoms(mol, parent),
        _ketone_fg_atoms(mol, parent),
        _amine_fg_atoms(mol, parent),
        _nitrile_fg_atoms(mol, parent),
        _thiol_fg_atoms(mol, parent),
        _ester_fg_atoms(mol, parent),
        _anhydride_fg_atoms(mol, parent),
        _acyl_chloride_fg_atoms(mol, parent),
    )
    out: set[int] = set()
    for part in parts:
        out |= part
    return out


def compute_owned_atoms(parent: dict, mol: Mol) -> frozenset[int]:
    """Union chain + kind-specific FG atoms (terminal ownership set)."""
    return frozenset(_chain_atoms(parent) | _kind_fg_atoms(parent, mol))


def finalize_parent_ownership(parent: dict, mol: Mol) -> dict:
    """Copy candidate once with immutable owned_atoms frozenset."""
    if isinstance(parent.get("owned_atoms"), frozenset):
        return parent
    return {**parent, "owned_atoms": compute_owned_atoms(parent, mol)}
