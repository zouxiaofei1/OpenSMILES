"""Typed L3 substituent namer: ordered retained → rooted_tree → recursive."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence

from rdkit.Chem import BondType

from namepredict.layer2.claimable_block import ClaimedBlock
from namepredict.layer3.as_substituent import name_as_substituent


@dataclass(frozen=True)
class SubstituentName:
    claim: ClaimedBlock
    en: str
    zh: str
    requires_parentheses: bool
    backend: str  # "retained", "rooted_tree", or "recursive"


class SubstituentBackend(Protocol):
    name: str

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None: ...


def _retained_hit(claim: ClaimedBlock, en: str, zh: str, paren: bool) -> SubstituentName:
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="retained",
    )


def _try_cycloalkyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.side_cycloalkyl import _cycloalkyl_names, _is_monocycloalkyl

    parent = {claim.attach_parent}
    atoms = _is_monocycloalkyl(mol, claim.root, parent)
    if atoms is None or set(atoms) != set(claim.atoms):
        return None
    names = _cycloalkyl_names(len(atoms))
    return None if names is None else _retained_hit(claim, names[0], names[1], False)


def _try_sub_phenyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.aryl_sub import _phenyl_at
    from namepredict.layer3.ring_namer import recursive_ph_name

    ph = _phenyl_at(mol, claim.root, claim.attach_parent)
    if ph is None:
        return None
    en, zh, paren, atoms = recursive_ph_name(mol, ph, claim.root, claim.attach_parent)
    return None if set(atoms) != set(claim.atoms) or not en else _retained_hit(claim, en, zh, paren)


def _methoxy_o(mol, claim: ClaimedBlock) -> int | None:
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() == 8:
        return claim.root
    return next((i for i in claim.atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 8), None)


def _is_plain_methyl_on_o(mol, o_idx: int, c_idx: int) -> bool:
    if mol.GetAtomWithIdx(c_idx).GetAtomicNum() != 6:
        return False
    heavies = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == o_idx


def _try_methoxy(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """O–CH3 side attached at claim.attach_parent (ether O or chain C via O root)."""
    if len(claim.atoms) != 2:
        return None
    o_idx = _methoxy_o(mol, claim)
    if o_idx is None:
        return None
    c_idxs = [i for i in claim.atoms if i != o_idx]
    if len(c_idxs) != 1 or not _is_plain_methyl_on_o(mol, o_idx, c_idxs[0]):
        return None
    return _retained_hit(claim, "methoxy", "甲氧基", False)


def _methylsulfanyl_s(mol, claim: ClaimedBlock) -> int | None:
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() == 16:
        return claim.root
    return next((i for i in claim.atoms if mol.GetAtomWithIdx(i).GetAtomicNum() == 16), None)


def _is_plain_methyl_on_s(mol, s_idx: int, c_idx: int) -> bool:
    if mol.GetAtomWithIdx(c_idx).GetAtomicNum() != 6:
        return False
    heavies = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() if n.GetAtomicNum() != 1]
    return len(heavies) == 1 and heavies[0].GetIdx() == s_idx


def _try_methylsulfanyl(mol, claim: ClaimedBlock) -> SubstituentName | None:
    """CH3–S side attached at claim.attach_parent (thioether S or chain C via S root)."""
    if len(claim.atoms) != 2:
        return None
    s_idx = _methylsulfanyl_s(mol, claim)
    if s_idx is None:
        return None
    c_idxs = [i for i in claim.atoms if i != s_idx]
    if len(c_idxs) != 1 or not _is_plain_methyl_on_s(mol, s_idx, c_idxs[0]):
        return None
    # Reject if S has any =O neighbour (sulfinyl/sulfonyl territory)
    if _s_dbl_o(mol, s_idx) != 0:
        return None
    return _retained_hit(claim, "methylsulfanyl", "甲硫基", False)


def _s_dbl_o(mol, s_idx: int) -> int:
    """Count doubly-bonded oxygen neighbours on sulfur."""
    return sum(1 for n in mol.GetAtomWithIdx(s_idx).GetNeighbors()
               if n.GetAtomicNum() == 8
               and mol.GetBondBetweenAtoms(s_idx, n.GetIdx()).GetBondType() == BondType.DOUBLE)


def _try_methylsulfinyl(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """CH3–S(=O)– : S with methyl and exactly 1 oxo neighbour."""
    s_idx = _methylsulfanyl_s(mol, claim)
    if s_idx is None:
        return None
    c_idx = next((i for i in claim.atoms if i != s_idx
                  and mol.GetAtomWithIdx(i).GetAtomicNum() == 6), None)
    if c_idx is None or not _is_plain_methyl_on_s(mol, s_idx, c_idx):
        return None
    if _s_dbl_o(mol, s_idx) != 1:
        return None
    from namepredict.layer3.retained_substituents import resolve_name
    en, zh = resolve_name("methylsulfinyl", name_mode=name_mode)
    return _retained_hit(claim, en, zh, False)


def _try_methylsulfonyl(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """CH3–S(=O)(=O)– PIN name 'methylsulfonyl'."""
    s_idx = _methylsulfanyl_s(mol, claim)
    if s_idx is None:
        return None
    c_idx = next((i for i in claim.atoms if i != s_idx
                  and mol.GetAtomWithIdx(i).GetAtomicNum() == 6), None)
    if c_idx is None or not _is_plain_methyl_on_s(mol, s_idx, c_idx):
        return None
    if _s_dbl_o(mol, s_idx) != 2:
        return None
    from namepredict.layer3.retained_substituents import resolve_name
    en, zh = resolve_name("methylsulfonyl", name_mode=name_mode)
    return _retained_hit(claim, en, zh, False)


def _try_tosyl(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """4-methylphenyl–S(=O)(=O)– retained name 'tosyl'."""
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() != 16:
        return None
    s_idx = claim.root
    if _s_dbl_o(mol, s_idx) != 2:
        return None
    # Aromatic C must be in-claim (not parent-owned)
    arom_c = [n for n in mol.GetAtomWithIdx(s_idx).GetNeighbors()
              if n.GetAtomicNum() == 6 and n.GetIsAromatic()
              and n.GetIdx() in claim.atoms]
    if len(arom_c) != 1:
        return None
    # Require at least 6 aromatic ring carbons in claim
    ring_c = [i for i in claim.atoms
              if mol.GetAtomWithIdx(i).GetIsAromatic() and mol.GetAtomWithIdx(i).GetAtomicNum() == 6]
    if len(ring_c) < 6:
        return None
    from namepredict.layer3.retained_substituents import resolve_name
    en, zh = resolve_name("tosyl", name_mode=name_mode)
    return _retained_hit(claim, en, zh, False)


def _try_triflyl(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """CF3–S(=O)(=O)– retained name 'triflyl'."""
    if mol.GetAtomWithIdx(claim.root).GetAtomicNum() != 16:
        return None
    # Triflyl needs at least S+2O+C+3F = 7 atoms
    if len(claim.atoms) < 5:
        return None
    s_idx = claim.root
    if _s_dbl_o(mol, s_idx) != 2:
        return None
    c_nb = [n for n in mol.GetAtomWithIdx(s_idx).GetNeighbors() if n.GetAtomicNum() == 6
            and sum(1 for nn in n.GetNeighbors() if nn.GetAtomicNum() == 9) == 3]
    if len(c_nb) != 1:
        return None
    from namepredict.layer3.retained_substituents import resolve_name
    en, zh = resolve_name("triflyl", name_mode=name_mode)
    return _retained_hit(claim, en, zh, False)


def _try_registry_leaf(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    """Match claim against registry entries that have leaf_atoms topology fingerprint."""
    from namepredict.layer3.retained_substituents import _REGISTRY, resolve_name

    claim_z = sorted(mol.GetAtomWithIdx(i).GetAtomicNum() for i in claim.atoms)
    root_z = mol.GetAtomWithIdx(claim.root).GetAtomicNum()
    _BO = {BondType.SINGLE: 1.0, BondType.DOUBLE: 2.0, BondType.TRIPLE: 3.0, BondType.AROMATIC: 1.5}

    for key, entry in _REGISTRY.items():
        if entry.leaf_atoms is None:
            continue
        if tuple(claim_z) != entry.leaf_atoms:
            continue
        if entry.leaf_root_z is not None and root_z != entry.leaf_root_z:
            continue
        if entry.leaf_bond_order is not None:
            if len(claim.atoms) != 2:
                continue
            bond = mol.GetBondBetweenAtoms(*claim.atoms)
            if bond is None or _BO.get(bond.GetBondType()) != entry.leaf_bond_order:
                continue
        if entry.validate is not None and not entry.validate(mol, claim):
            continue
        en, zh = resolve_name(key, name_mode=name_mode)
        return _retained_hit(claim, en, zh, False)
    return None


_ALKENYL_CHECKS = (
    ("vinyl",),
    ("allyl",),
    ("isopropenyl",),
)


def _try_alkenyl_retained(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    from namepredict.layer2 import side_facts
    from namepredict.layer3.retained_substituents import resolve_name

    SHAPES = {
        "vinyl": side_facts.AlkylShape.C2_VINYL,
        "allyl": side_facts.AlkylShape.C3_ALLYL,
        "isopropenyl": side_facts.AlkylShape.C3_ISOPROPENYL,
    }
    root, parent = claim.root, {claim.attach_parent}
    for (key,) in _ALKENYL_CHECKS:
        fact = side_facts.alkyl_shape(mol, root, parent, SHAPES[key])
        if fact and set(fact.atoms) == set(claim.atoms):
            en, zh = resolve_name(key, name_mode=name_mode)
            return _retained_hit(claim, en, zh, False)
    return None


def _retained_name(mol, claim: ClaimedBlock, *, name_mode: str = "general") -> SubstituentName | None:
    return (
        _try_registry_leaf(mol, claim, name_mode=name_mode)
        or _try_cycloalkyl(mol, claim)
        or _try_sub_phenyl(mol, claim)
        or _try_alkenyl_retained(mol, claim, name_mode=name_mode)
        or _try_methoxy(mol, claim)
        or _try_methylsulfanyl(mol, claim)
        or _try_methylsulfinyl(mol, claim, name_mode=name_mode)
        or _try_methylsulfonyl(mol, claim, name_mode=name_mode)
        or _try_tosyl(mol, claim, name_mode=name_mode)
        or _try_triflyl(mol, claim, name_mode=name_mode)
    )


class RetainedBackend:
    """Retained simple leaves: cycloalkyl, substituted phenyl, methoxy."""

    name = "retained"

    def __init__(self, *, name_mode: str = "general") -> None:
        self._name_mode = name_mode

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return _retained_name(mol, claim, name_mode=self._name_mode)


def _rooted_tree_name(mol, claim: ClaimedBlock) -> SubstituentName | None:
    from namepredict.layer2.side_alkyl_sys import build_rooted_alkyl_tree
    from namepredict.layer3.alkyl_sys_names import name_rooted_alkyl

    tree = build_rooted_alkyl_tree(mol, root=claim.root, atoms=claim.atoms)
    hit = None if tree is None else name_rooted_alkyl(tree)
    if hit is None or not hit[0] or not hit[1]:
        return None
    en, zh, paren = hit
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="rooted_tree",
    )


class RootedTreeBackend:
    """Pure saturated-carbon rooted tree (max_atoms=12, max_depth=3)."""

    name = "rooted_tree"

    def __init__(self, *, name_mode: str = "general") -> None:
        self._name_mode = name_mode

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        return _rooted_tree_name(mol, claim)


class RecursiveBackend:
    """Bounded recursive cut → free-name → yl_form."""

    name = "recursive"

    def __init__(self, *, name_mode: str = "general") -> None:
        self._name_mode = name_mode

    def try_name(self, mol, claim: ClaimedBlock, *, depth: int) -> SubstituentName | None:
        hit = name_as_substituent(mol, claim.root, claim.atoms, depth=depth, name_mode=self._name_mode)
        return None if hit is None else _from_yl(claim, hit)


def _from_yl(claim: ClaimedBlock, hit: tuple[str, str, bool]) -> SubstituentName | None:
    en, zh, paren = hit
    if not en or not zh:
        return None
    return SubstituentName(
        claim=claim, en=en, zh=zh, requires_parentheses=paren, backend="recursive",
    )


def _default_backends(name_mode: str = "general") -> list[SubstituentBackend]:
    return [RetainedBackend(name_mode=name_mode), RootedTreeBackend(name_mode=name_mode), RecursiveBackend(name_mode=name_mode)]


class SubstituentNamer:
    def __init__(self, backends: Sequence[SubstituentBackend] | None = None, *, name_mode: str = "general") -> None:
        self._backends = list(backends) if backends is not None else _default_backends(name_mode)

    def name(self, mol, claim: ClaimedBlock, *, depth: int = 0) -> SubstituentName | None:
        for backend in self._backends:
            hit = backend.try_name(mol, claim, depth=depth)
            if hit is not None:
                return hit
        return None
