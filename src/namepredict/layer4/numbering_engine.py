"""Candidate-based parent numbering engine (P-14.4).

Replaces kind-dispatched orienters with a uniform pipeline:

  1. enumerate numbering candidates
       chain: forward / reversed
       ring : every atom as locant 1, both directions (2n candidates)
  2. narrow candidates by P-14.4 rules in order:
       (c) principal characteristic group  → lowest locant set
       (e) multiple bonds                  → lowest locant set (double first)
       (f) substituents                    → lowest locant set
  3. early-stop when a single candidate survives

Each candidate is an ``{atom: locant}`` dict; the survivor is emitted back as
an atom-order chain (locant = index + 1), the same shape `_orient_chain`
returned, so downstream locant packing is untouched.
"""
from __future__ import annotations

from namepredict.layer4.locants.adapt import plan_from_chain


# ── candidate generation ──────────────────────────────────────────────────

def _numbered(chain: list[int]) -> dict[int, int]:
    return {a: i + 1 for i, a in enumerate(chain)}


def _to_chain(cand: dict[int, int]) -> list[int]:
    return sorted(cand, key=cand.get)


def _chain_cands(chain: list[int]) -> list[dict[int, int]]:
    return [_numbered(chain), _numbered(list(reversed(chain)))]


def _ring_cands(chain: list[int]) -> list[dict[int, int]]:
    n = len(chain)
    out: list[dict[int, int]] = []
    for i in range(n):
        fwd = chain[i:] + chain[:i]
        out.append(_numbered(fwd))
        out.append(_numbered(list(reversed(fwd))))
    return out


# ── locant-set keys ───────────────────────────────────────────────────────

def _locant_set(cand: dict[int, int], atoms: list[int]) -> tuple[int, ...] | None:
    locs = sorted(cand[a] for a in atoms if a in cand)
    return tuple(locs) if locs else None


def _bond_locants(cand: dict[int, int], bonds) -> tuple[int, ...] | None:
    if not bonds:
        return None
    mins = []
    for b in bonds:
        if b[0] in cand and b[1] in cand:
            mins.append(min(cand[b[0]], cand[b[1]]))
    return tuple(sorted(mins)) if mins else None


def _narrow(cands: list[dict], key_fn) -> list[dict]:
    """Keep candidates with the lowest locant-set key; early-stop at one."""
    if len(cands) <= 1:
        return cands
    keys = [key_fn(c) for c in cands]
    if any(k is None for k in keys):
        return cands  # feature absent everywhere → rule does not apply
    best = min(keys)
    return [c for c, k in zip(cands, keys) if k == best]


# ── P-14.4 feature extraction from the parent dict ────────────────────────

def _principal_atoms(parent: dict) -> list[int]:
    """P-14.4(c): principal characteristic group attachment atoms."""
    atoms = parent.get("principal_attachment_atoms")
    if atoms:
        return list(atoms)
    facts = parent.get("principal_expression_facts")
    if facts is not None:
        attach = getattr(facts, "attachment_atoms", None)
        if attach:
            return sorted(attach)
    return []


def _unsat_bonds(parent: dict) -> tuple[list, list]:
    """(all multiple bonds, double bonds) endpoint pairs."""
    all_bonds, doubles = [], []
    for key, target in (("double_bond", doubles), ("triple_bond", all_bonds)):
        v = parent.get(key)
        if v:
            pair = (v[0], v[1])
            (target if key == "double_bond" else all_bonds).append(pair)
    for key in ("double_bonds", "triple_bonds"):
        for b in parent.get(key) or []:
            pair = (b[0], b[1])
            all_bonds.append(pair)
            if key == "double_bonds":
                doubles.append(pair)
    return all_bonds, doubles


def _is_ring(parent: dict) -> bool:
    return bool(parent.get("scaffold_id"))


# P-14.4(a): parent dict fields naming an atom that must be locant 1
# (heteroatom ring starts, exocyclic carbonyl attach, radical centre).
_FIXED_START_KEYS = (
    "ring_attach_idx", "n_idx", "nh_idx", "hetero_idx", "radical_c_idx",
)


def _ring_hetero_start(parent: dict, chain: list[int]) -> int | None:
    """Monohetero ring: the single heteroatom is locant 1 (P-14.4, e.g. pyridine)."""
    mol = parent.get("mol")
    if mol is None:
        return None
    heteros = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    return heteros[0] if len(heteros) == 1 else None


def _fixed_start(parent: dict) -> int | None:
    for key in _FIXED_START_KEYS:
        v = parent.get(key)
        if v is not None:
            return v
    return _ring_hetero_start(parent, parent.get("chain") or [])


def _fixed_numbering(parent: dict, chain: list[int]) -> list[int] | None:
    """P-14.4(a): retained-scaffold fixed numbering via L2 numbering_scaffold."""
    plan = plan_from_chain(
        chain, parent.get("scaffold_id") or parent.get("kind"),
        parent.get("numbering_scaffold"),
        required=parent.get("numbering_scaffold_required", False),
    )
    return list(plan.atom_order) if plan is not None else None


# ── entry ─────────────────────────────────────────────────────────────────

def orient_numbering(parent: dict, substituents: list) -> list[int] | None:
    """Return the P-14.4-oriented atom order, or None if not applicable."""
    chain = parent.get("chain") or []
    if not chain:
        return None
    fixed = _fixed_numbering(parent, chain)
    if fixed is not None:
        return fixed
    if _is_ring(parent):
        cands = _ring_cands(chain)
    else:
        cands = _chain_cands(chain)
    start = _fixed_start(parent)
    if start is not None:
        cands = [c for c in cands if c.get(start) == 1]
    principal = _principal_atoms(parent)
    if principal:
        cands = _narrow(cands, lambda c: _locant_set(c, principal))
    bonds, doubles = _unsat_bonds(parent)
    if bonds:
        cands = _narrow(cands, lambda c: (_bond_locants(c, bonds), _bond_locants(c, doubles)))
    subs = [s["attach_idx"] for s in substituents if s.get("attach_idx") in chain]
    if subs:
        cands = _narrow(cands, lambda c: _locant_set(c, subs))
    if len(cands) > 1 and substituents:
        # P-14.4(f) tie: lowest locant set already equal — assign lowest locant
        # to the alphabetically-first substituent (stem-alpha pairs, P-14.5).
        from namepredict.layer4._chain_orient import _stem_loc_pairs
        cands = _narrow(cands, lambda c: _stem_loc_pairs(_to_chain(c), substituents))
    return _to_chain(cands[0])
