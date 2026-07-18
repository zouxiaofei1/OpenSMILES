"""Bootstrap ring parent producers into kind_registry (single registration site).

candidates._ring_candidates iterates kind_registry.ring_try_fns(); add new
fused/hetero ring parents here instead of editing a hand tuple in candidates.
"""
from __future__ import annotations

from namepredict.layer2.anthracene import _try_anthracene_parent
from namepredict.layer2.anthraquinone import _try_anthraquinone_parent
from namepredict.layer2.azole13 import _try_azole13_parent
from namepredict.layer2.benzoquinone import _try_benzoquinone_parent
from namepredict.layer2.ortho_benzoquinone import _try_ortho_benzoquinone_parent
from namepredict.layer2.benzimidazole import (
    _try_benzimidazolamine_parent,
    _try_benzimidazole_parent,
)
from namepredict.layer2.benzodiazine import (
    _try_quinazoline_parent,
    _try_quinoxaline_parent,
)
from namepredict.layer2.benzofuran import _try_benzofuran_parent
from namepredict.layer2.benzothiazole import _try_benzothiazole_parent
from namepredict.layer2.benzothiophene import _try_benzothiophene_parent
from namepredict.layer2.benzoxazole import _try_benzoxazole_parent
from namepredict.layer2.chromenone import _try_chromenone_parent
from namepredict.layer2.heteroarene5 import (
    _try_diazine_parent,
    _try_hetero5_parent,
    _try_imidazole_parent,
    _try_pyrazole_parent,
)
from namepredict.layer2.indazole import _try_indazole_parent
from namepredict.layer2.indole import _try_indole_parent
from namepredict.layer2.kind_registry import register_ring_try
from namepredict.layer2.naphthalene import _try_naphthalene_parent
from namepredict.layer2.parent_core import _parent_dict
from namepredict.layer2.parent_selector import (
    _benzene_parent,
    _cycloalkane_parent,
)
from namepredict.layer2.pyridine import _try_pyridine_parent
from namepredict.layer2.quinoline import (
    _try_isoquinoline_parent,
    _try_quinoline_parent,
)
from namepredict.layer2.ring_parent import (
    _endocyclic_doubles,
    _is_simple_benzene,
    _is_simple_cycloalkane,
    _is_simple_cyclopolyene,
)
from namepredict.layer2.sat_hetero import _try_sat_hetero_parent
from namepredict.layer2.scaffold.builders.sat_hetero_repl import try_sat_hetero_repl


def _try_simple_benzene(info: dict) -> dict | None:
    return _benzene_parent(info) if _is_simple_benzene(info) else None


def _try_simple_cycloalkane(info: dict) -> dict | None:
    return _cycloalkane_parent(info) if _is_simple_cycloalkane(info) else None


def _cyclopolyene_parent(info: dict) -> dict:
    chain = list(info["rings"][0]["atom_ids"])
    bonds = _endocyclic_doubles(info, set(chain))
    return _parent_dict(chain, "cyclopolyene", double_bonds=bonds)


def _try_simple_cyclopolyene(info: dict) -> dict | None:
    return _cyclopolyene_parent(info) if _is_simple_cyclopolyene(info) else None


_RING_PRODUCERS = (
    _try_anthraquinone_parent,
    _try_benzoquinone_parent,
    _try_ortho_benzoquinone_parent,
    _try_anthracene_parent,
    _try_chromenone_parent,
    _try_quinazoline_parent,
    _try_quinoxaline_parent,
    _try_naphthalene_parent,
    _try_indole_parent,
    _try_indazole_parent,
    _try_benzofuran_parent,
    _try_benzothiophene_parent,
    _try_benzothiazole_parent,
    _try_benzoxazole_parent,
    _try_benzimidazolamine_parent,
    _try_benzimidazole_parent,
    _try_quinoline_parent,
    _try_isoquinoline_parent,
    _try_pyridine_parent,
    _try_diazine_parent,
    _try_imidazole_parent,
    _try_pyrazole_parent,
    _try_azole13_parent,
    _try_hetero5_parent,
    _try_sat_hetero_parent,
    try_sat_hetero_repl,
    _try_simple_benzene,
    _try_simple_cycloalkane,
    _try_simple_cyclopolyene,
)


def _bootstrap() -> None:
    for fn in _RING_PRODUCERS:
        register_ring_try(fn)


_bootstrap()
