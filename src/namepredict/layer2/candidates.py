"""Layer2 parent-candidate collection (scored by layer2.scoring, P-44).

FG and ring classes each try independently (no short-circuit `or`).
Scoring then picks the best among all viable parents so `_FG_RANK`
actually arbitrates acid vs alcohol vs amine, etc.

Ring parent producers are registered in `layer2.ring_producers` →
`kind_registry.ring_try_fns()` (not a hand-written tuple here).
"""
from __future__ import annotations

from namepredict.layer2 import kind_registry as _kr
from namepredict.layer2 import ring_producers as _ring_producers  # noqa: F401  bootstrap
from namepredict.layer2.arene_carbonyl import _try_arene_other_fg
from namepredict.layer2.benzenediamine import (
    _benzenediamine_parent,
    _is_simple_benzenediamine,
)
from namepredict.layer2.phosphate import _phosphate_parent, _phosphonic_parent
from namepredict.layer2.parent_selector import (
    _acid_parent,
    _acyl_chloride_parent,
    _alcohol_parent,
    _aldehyde_parent,
    _alkene_parent,
    _alkyne_parent,
    _amide_parent,
    _amine_parent,
    _anhydride_parent,
    _benzene_parent,
    _ester_parent,
    _ether_parent,
    _is_mono_alkene,
    _is_mono_alkyne,
    _is_mono_fg,
    _is_polyene,
    _is_sym_anhydride,
    _ketone_parent,
    _longest_chain,
    _nitrile_parent,
    _parent_dict,
    _polyene_parent,
    _sulfide_parent,
    _thiol_parent,
)
from namepredict.layer2.ring_parent import (
    _is_benzene_core,
    _is_simple_benzene,
)
from namepredict.layer2.scoring import _pick_best

# silence unused: import side-effect registers producers
assert _ring_producers is not None


def _benzene_candidate(info: dict) -> dict | None:
    if not _is_benzene_core(info):
        return None
    cand = _benzene_parent(info)
    if not _is_simple_benzene(info):
        cand["n_unhandled"] = 1
    return cand


def _alkane_fallback(info: dict) -> dict:
    return _parent_dict(_longest_chain(info["mol"]), "alkane")


def _try_acid(info: dict) -> dict | None:
    if info.get("has_acid") and info.get("carboxyls"):
        return _acid_parent(info)
    return None


def _try_anhydride(info: dict) -> dict | None:
    return _anhydride_parent(info) if _is_sym_anhydride(info) else None


def _try_acyl_chloride(info: dict) -> dict | None:
    if _is_mono_fg(info, "has_acyl_chloride", "acyl_chlorides"):
        return _acyl_chloride_parent(info)
    return None


def _try_ester(info: dict) -> dict | None:
    if _is_mono_fg(info, "has_ester", "esters"):
        return _ester_parent(info)
    return None


def _try_amide(info: dict) -> dict | None:
    if _is_mono_fg(info, "has_amide", "amides"):
        return _amide_parent(info)
    return None


def _is_aryl_nitrile_c(mol, c_idx: int) -> bool:
    """True if nitrile carbon's only carbon neighbor is aromatic (Ar–CN leaf)."""
    atom = mol.GetAtomWithIdx(c_idx)
    cs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
    return len(cs) == 1 and cs[0].GetIsAromatic()


def _try_nitrile(info: dict) -> dict | None:
    if not _is_mono_fg(info, "has_nitrile", "nitriles"):
        return None
    nit = info["nitriles"][0]
    if _is_aryl_nitrile_c(info["mol"], nit["c_idx"]):
        return None  # Ar–CN is a leaf; benzonitrile via arene path
    return _nitrile_parent(info)


def _try_aldehyde(info: dict) -> dict | None:
    if info.get("has_aldehyde") and info.get("aldehydes"):
        return _aldehyde_parent(info)
    return None


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


def _try_benzenediamine(info: dict) -> dict | None:
    if _is_simple_benzenediamine(info):
        return _benzenediamine_parent(info)
    return None


_FG_TRY = (
    _try_acid, _try_anhydride, _try_arene_other_fg, _try_acyl_chloride,
    _try_ester, _try_amide, _try_nitrile, _try_aldehyde, _try_ketone,
    _try_alcohol, _try_thiol, _try_benzenediamine, _try_amine,
    _try_phosphate, _try_phosphonic, _ether_parent, _sulfide_parent,
)


def _fg_candidates(info: dict) -> list[dict]:
    return [c for fn in _FG_TRY if (c := fn(info)) is not None]


def _ring_candidates(info: dict) -> list[dict]:
    return [c for fn in _kr.ring_try_fns() if (c := fn(info)) is not None]


def _unsat_candidates(info: dict) -> list[dict]:
    out: list[dict] = []
    if _is_mono_alkyne(info):
        out.append(_alkyne_parent(info))
    if _is_polyene(info):
        out.append(_polyene_parent(info))
    if _is_mono_alkene(info):
        out.append(_alkene_parent(info))
    return out


def _dedupe_parents(cands: list[dict]) -> list[dict]:
    seen: set[tuple] = set()
    out: list[dict] = []
    for c in cands:
        key = (c.get("kind"), tuple(c.get("chain") or []))
        if key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


def _collect_candidates(info: dict) -> list[dict]:
    raw = (
        _fg_candidates(info) + _ring_candidates(info)
        + [_benzene_candidate(info)] + _unsat_candidates(info)
        + [_alkane_fallback(info)]
    )
    return _dedupe_parents([c for c in raw if c is not None])


def _fg_parent(info: dict) -> dict | None:
    """Compat: best FG among independent class tries."""
    return _pick_best(info, _fg_candidates(info))


def _ring_parent(info: dict) -> dict | None:
    """Compat: best ring among independent class tries."""
    return _pick_best(info, _ring_candidates(info))
