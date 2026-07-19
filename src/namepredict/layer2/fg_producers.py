"""Bootstrap FG parent producers into kind_registry (single registration site).

candidates._fg_candidates iterates kind_registry.fg_try_fns(); add new
open-chain / principal FG parents here instead of editing a hand tuple in
candidates. Keep try order stable — scoring (_FG_RANK) arbitrates, not or.
"""
from __future__ import annotations

from namepredict.layer2.arene_carbonyl import _try_arene_other_fg
from namepredict.layer2.benzenediamine import (
    _benzenediamine_parent,
    _is_simple_benzenediamine,
)
from namepredict.layer2.boronic import _try_boronic
from namepredict.layer2.carbamate import _carbamate_parent
from namepredict.layer2.carbonate import _carbonate_parent
from namepredict.layer2.diester import _diester_parent
from namepredict.layer2.guanidine import _guanidine_parent
from namepredict.layer2.hydrazine import _hydrazine_parent
from namepredict.layer2.isocyanate import _try_isocyanate, _try_isothiocyanate
from namepredict.layer2.kind_registry import register_fg_try
from namepredict.layer2.parent_core import _is_mono_fg
from namepredict.layer2.parent_selector import (
    _acid_parent,
    _acyl_chloride_parent,
    _alcohol_parent,
    _aldehyde_parent,
    _amide_parent,
    _amine_parent,
    _anhydride_parent,
    _ester_parent,
    _ether_parent,
    _is_sym_anhydride,
    _ketone_parent,
    _nitrile_parent,
    _sulfide_parent,
    _thiol_parent,
    _try_izcn,
    _try_pycn,
)
from namepredict.layer2.phosphate import _phosphate_parent, _phosphonic_parent
from namepredict.layer2.sulfonamide import _sulfonamide_parent
from namepredict.layer2.sulfonate import _sulfonate_parent
from namepredict.layer2.sulfone import _sulfone_parent
from namepredict.layer2.sulfonic_acid import _sulfonic_acid_parent
from namepredict.layer2.sulfonyl_chloride import _sulfonyl_chloride_parent
from namepredict.layer2.sulfoxide import _sulfoxide_parent
from namepredict.layer2.urea import _urea_parent
from namepredict.layer2.cyclo_polycarboxylic import _try_cycloalkane_polycarboxylic
from namepredict.layer2.cyclo_carboxylic import _try_cycloalkane_exocyclic_fg
from namepredict.layer2.sat_hetero_one import _try_sat_hetero_one_parent
from namepredict.layer2.chromenone import _try_chromenone_parent


def _try_acid(info: dict) -> dict | None:
    if info.get("has_acid") and info.get("carboxyls"):
        return _acid_parent(info)
    return None


def _try_anhydride(info: dict) -> dict | None:
    return _anhydride_parent(info) if _is_sym_anhydride(info) else None


def _try_cyclo_exo_or(info: dict, open_fn) -> dict | None:
    """Prefer cycloalkane exocyclic FG parent over open-chain formic collapse."""
    if (c := _try_cycloalkane_exocyclic_fg(info)) is not None:
        return c
    return open_fn()


def _try_acyl_chloride(info: dict) -> dict | None:
    """Acyl chloride / bromide (kind set inside parent by halogen)."""
    if not _is_mono_fg(info, "has_acyl_chloride", "acyl_chlorides"):
        return None
    return _try_cyclo_exo_or(info, lambda: _acyl_chloride_parent(info))


def _try_diester(info: dict) -> dict | None:
    return _diester_parent(info)


def _try_lactone_or(info: dict, open_fn) -> dict | None:
    """Prefer fused/sat lactone/lactam over open-chain ester/amide collapse."""
    if (c := _try_chromenone_parent(info)) is not None:
        return c
    if (c := _try_sat_hetero_one_parent(info)) is not None:
        return c
    return _try_cyclo_exo_or(info, open_fn)


def _try_ester(info: dict) -> dict | None:
    if not _is_mono_fg(info, "has_ester", "esters"):
        return None
    return _try_lactone_or(info, lambda: _ester_parent(info))


def _try_amide(info: dict) -> dict | None:
    if not _is_mono_fg(info, "has_amide", "amides"):
        return None
    return _try_lactone_or(info, lambda: _amide_parent(info))


def _is_aryl_nitrile_c(mol, c_idx: int) -> bool:
    """True if nitrile carbon's only carbon neighbor is aromatic (Ar–CN leaf)."""
    atom = mol.GetAtomWithIdx(c_idx)
    cs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(cs) == 1 and cs[0].GetIsAromatic()


def _try_heteroaryl_carbonitrile(info: dict) -> dict | None:
    """Retained heteroaryl–CN parents (indazole / pyridine carbonitrile)."""
    return _try_izcn(info) or _try_pycn(info)


def _try_nitrile(info: dict) -> dict | None:
    if not _is_mono_fg(info, "has_nitrile", "nitriles"):
        return None
    if (h := _try_heteroaryl_carbonitrile(info)) is not None:
        return h
    nit = info["nitriles"][0]
    if _is_aryl_nitrile_c(info["mol"], nit["c_idx"]):
        return None  # Ar–CN leaf; benzonitrile via arene path
    return _try_cyclo_exo_or(info, lambda: _nitrile_parent(info))


def _try_aldehyde(info: dict) -> dict | None:
    if not (info.get("has_aldehyde") and info.get("aldehydes")):
        return None
    return _try_cyclo_exo_or(info, lambda: _aldehyde_parent(info))


def _try_ketone(info: dict) -> dict | None:
    if info.get("has_ketone") and info.get("ketones"):
        return _ketone_parent(info)
    return None


def _try_alcohol(info: dict) -> dict | None:
    if info.get("has_alcohol") and info.get("hydroxyls"):
        return _alcohol_parent(info)
    return None


def _try_thiol(info: dict) -> dict | None:
    if info.get("has_thiol") and info.get("thiols"):
        return _thiol_parent(info)
    return None


def _try_amine(info: dict) -> dict | None:
    if info.get("has_amine") and info.get("amines"):
        return _amine_parent(info)
    return None


def _try_phosphate(info: dict) -> dict | None:
    return _phosphate_parent(info)


def _try_phosphonic(info: dict) -> dict | None:
    return _phosphonic_parent(info)


def _try_carbamate(info: dict) -> dict | None:
    return _carbamate_parent(info)


def _try_carbonate(info: dict) -> dict | None:
    return _carbonate_parent(info)


def _try_urea(info: dict) -> dict | None:
    return _urea_parent(info)


def _try_hydrazine(info: dict) -> dict | None:
    return _hydrazine_parent(info)


def _try_guanidine(info: dict) -> dict | None:
    return _guanidine_parent(info)


def _try_benzenediamine(info: dict) -> dict | None:
    if _is_simple_benzenediamine(info):
        return _benzenediamine_parent(info)
    return None


def _try_sulfone(info: dict) -> dict | None:
    if info.get("has_sulfone") and info.get("sulfones"):
        return _sulfone_parent(info)
    return None


# Order = historical candidates._FG_TRY; do not reorder without dual check.
_FG_PRODUCERS = (
    _try_cycloalkane_polycarboxylic,
    _try_acid,
    _try_anhydride,
    _try_arene_other_fg,
    _try_acyl_chloride,
    _try_diester,
    _try_ester,
    _try_carbamate,
    _try_carbonate,
    _try_urea,
    _try_guanidine,
    _try_amide,
    _try_nitrile,
    _try_aldehyde,
    _try_ketone,
    _try_alcohol,
    _try_thiol,
    _try_benzenediamine,
    _try_amine,
    _try_hydrazine,
    _try_phosphate,
    _try_phosphonic,
    _ether_parent,
    _sulfide_parent,
    _sulfoxide_parent,
    _try_sulfone,
    _try_isocyanate,
    _try_isothiocyanate,
    _sulfonamide_parent,
    _sulfonate_parent,
    _sulfonyl_chloride_parent,
    _sulfonic_acid_parent,
    _try_boronic,
)


def _bootstrap() -> None:
    for fn in _FG_PRODUCERS:
        register_fg_try(fn)


_bootstrap()
