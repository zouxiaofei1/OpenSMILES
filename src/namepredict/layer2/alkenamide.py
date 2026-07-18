"""Open-chain monoamide / alkenamide parents (IUPAC P-66.1.1 / P-31.1)."""
from __future__ import annotations

from rdkit.Chem import Mol

from namepredict.layer2.chain_walk import _longest_from
from namepredict.layer2.ring_parent import _unsub_phenyl_at

# Higher-priority FGs + peers that block open-chain alkenamide (amide allowed).
_ALKENAMIDE_BAD = (
    "has_acid", "has_ester", "has_ketone", "has_amine",
    "has_acyl_chloride", "has_anhydride", "has_thiol", "has_alcohol",
    "has_aldehyde", "has_nitrile",
)


def _amide_n_alkyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    from namepredict.layer2.parent_core import _arm_ok
    arms = [_longest_from(mol, c, set()) for c in cs]
    if not all(_arm_ok(mol, a, am["n_idx"]) for a in arms):
        return {}
    if len(arms) == 1:
        return {"n_alkyl_n": len(arms[0])}
    return {"n_alkyl_ns": [len(a) for a in arms]} if arms else {}


def _amide_n_phenyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    if len(cs) != 1:
        return {}
    return {"n_phenyl": True} if _unsub_phenyl_at(mol, cs[0], am["n_idx"]) else {}


def _amide_n_benzyl(mol: Mol, am: dict, cs: list[int]) -> dict:
    """N–CH2–Ph (Ph may carry simple leaves) → n_benzyl flag + ch2 idx."""
    from namepredict.layer2.aryl_sub import _ch2_ph_at
    if len(cs) != 1:
        return {}
    ch2 = cs[0]
    ph = _ch2_ph_at(mol, ch2, am["n_idx"])
    return {"n_benzyl": True, "n_benzyl_ch2": ch2} if ph is not None else {}


def _amide_core_atoms(mol: Mol, am: dict) -> frozenset[int]:
    """Amide C + N + dbl O (parent core for N-block cut before parent exists)."""
    from namepredict.layer2.ring_parent import _dbl_o_idx
    o = _dbl_o_idx(mol, am["c_idx"])
    return frozenset({am["c_idx"], am["n_idx"]} | ({o} if o is not None else set()))


def _amide_n_block(mol: Mol, am: dict, cs: list[int]) -> dict:
    """Single N–C with non-empty cut_block → topology claim for complex N."""
    from namepredict.layer2.block_cut import cut_block
    if len(cs) != 1:
        return {}
    root = cs[0]
    block = cut_block(mol, root, _amide_core_atoms(mol, am))
    return {"n_block": True, "n_block_root": root} if block else {}


def _amide_n_meta(info: dict) -> dict:
    ams = info.get("amides") or []
    if len(ams) != 1:
        return {}
    mol, am, cs = info["mol"], ams[0], ams[0].get("n_c_idxs") or []
    return (
        _amide_n_benzyl(mol, am, cs)
        or _amide_n_phenyl(mol, am, cs)
        or _amide_n_alkyl(mol, am, cs)
        or _amide_n_block(mol, am, cs)
    )


def _amide_parent(info: dict) -> dict:
    from namepredict.layer2.parent_core import _unsat_or_sat
    return _unsat_or_sat(
        info, "has_amide", "amides", _ALKENAMIDE_BAD, "amide", "amide",
        "amide_c_idx", **_amide_n_meta(info),
    )
