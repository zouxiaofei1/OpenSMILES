"""Composable parent-candidate scoring (simplified IUPAC P-44 seniority).

Each candidate parent dict is scored into a tuple (bigger wins):
(has_principal_fg, fg_class_rank, sides_ok, is_hetero_ring, is_carbo_ring,
 n_rings, ring_size, retained_bonus, n_unsat, n_carbons, -n_unhandled_side)

`sides_ok` sits above the ring bits on purpose: a ring candidate whose
side chains cannot be expressed as substituents (parent["n_unhandled"]>0)
must lose to a chain fallback instead of emitting a bare ring name.
"""
from __future__ import annotations

_FG_RANK = {
    "acid": 13, "diacid": 13, "alkenoic_acid": 13, "alkenedioic": 13,
    "benzoic": 13, "pyridinecarboxylic": 13, "cycloalkanecarboxylic": 13,
    "anhydride": 12,
    "ester": 11, "alkenoate": 11, "benzoate": 11,
    "acyl_chloride": 10, "benzoyl_chloride": 10,
    "amide": 9,
    "nitrile": 8, "alkenenitrile": 8, "benzonitrile": 8,
    "aldehyde": 7, "alkenal": 7, "benzaldehyde": 7,
    "ketone": 6, "dione": 6, "cycloketone": 6, "acetophenone": 6,
    "alcohol": 5, "alkenol": 5, "diol": 5, "triol": 5,
    "cycloalcohol": 5, "phenol": 5, "benzenediol": 5, "pyridinol": 5,
    "thiol": 4,
    "amine": 3, "diamine": 3, "sec_amine": 3, "tert_amine": 3,
    "cycloamine": 3, "aniline": 3, "pyridinamine": 3, "pyrimidinamine": 3,
    "ether": 2, "sulfide": 2,
}
_HETERO_RING = frozenset(
    {"pyridine", "furan", "thiophene", "pyrrole", "imidazole", "pyrazole",
     "pyrimidine", "pyrazine", "pyridazine", "indole",
     "aziridine", "oxirane", "oxolane", "oxane", "pyrrolidine", "piperidine",
     "morpholine", "piperazine", "dioxolane", "dioxane", "thiolane"}
)
_CARBO_RING = frozenset({"benzene", "cycloalkane", "cycloalkene", "naphthalene"})
_RETAINED = frozenset(
    {"phenol", "aniline", "benzoic", "benzaldehyde", "acetophenone",
     "benzoate", "benzonitrile", "benzoyl_chloride", "benzene", "naphthalene"}
    | set(_HETERO_RING)
)


def _has_principal_fg(kind: str) -> int:
    return 1 if kind in _FG_RANK else 0


def _fg_class_rank(kind: str) -> int:
    return _FG_RANK.get(kind, 0)


def _n_unhandled(parent: dict) -> int:
    return int(parent.get("n_unhandled") or 0)


def _sides_ok(parent: dict) -> int:
    return 0 if _n_unhandled(parent) else 1


def _is_hetero_ring(kind: str) -> int:
    return 1 if kind in _HETERO_RING else 0


def _is_carbo_ring(kind: str) -> int:
    return 1 if kind in _CARBO_RING else 0


def _n_rings(kind: str) -> int:
    if kind in ("naphthalene", "indole"):
        return 2
    return 1 if kind in _HETERO_RING or kind in _CARBO_RING else 0


def _ring_size(parent: dict, kind: str) -> int:
    if not _n_rings(kind):
        return 0
    return len(parent.get("chain") or [])


def _n_unsat(parent: dict) -> int:
    kind = parent.get("kind")
    if kind == "polyene":
        return len(parent.get("double_bonds") or [])
    return 1 if kind in ("alkene", "alkyne", "cycloalkene") else 0


def _n_carbons(parent: dict) -> int:
    return int(parent.get("n_carbons") or len(parent.get("chain") or []))


def _retained_bonus(kind: str) -> int:
    return 1 if kind in _RETAINED else 0


def _score_parent(info: dict, parent: dict) -> tuple:
    kind = parent.get("kind") or ""
    return (
        _has_principal_fg(kind), _fg_class_rank(kind), _sides_ok(parent),
        _is_hetero_ring(kind), _is_carbo_ring(kind),
        _n_rings(kind), _ring_size(parent, kind),
        _retained_bonus(kind), _n_unsat(parent), _n_carbons(parent),
        -_n_unhandled(parent),
    )


def _better_parent(info: dict, a: dict, b: dict) -> bool:
    return _score_parent(info, a) > _score_parent(info, b)


def _pick_best(info: dict, candidates: list[dict]) -> dict | None:
    best: dict | None = None
    for cand in candidates:
        if cand is None:
            continue
        if best is None or _better_parent(info, cand, best):
            best = cand
    return best
