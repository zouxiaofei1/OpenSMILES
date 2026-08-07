"""Layer 3 的侧链拓扑事实：Layer 2 主链选择与 Layer 3 取代基命名共用。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


from rdkit.Chem import Atom, Mol

from namepredict.layer3 import aryl_sub, side_alkoxy, side_alkyl
from namepredict.layer3.leaves.protocol import ArylLeafKind, LeafTopology
from namepredict.layer3.leaves.registry import match_leaf_topology


class AlkylShape(Enum):
    C2_VINYL = auto()
    C3_ALLYL = auto()
    C3_ISOPROPENYL = auto()
    C3_BRANCH_AT_ROOT = auto()
    C4_BRANCH_AT_SECOND = auto()
    C4_TRIPLE_BRANCH_AT_ROOT = auto()
    C4_BRANCH_AFTER_ROOT = auto()
    C5_ASYMMETRIC_ROOT_BRANCH = auto()
    C5_BRANCH_NEAR_LEAF = auto()
    C5_DOUBLE_BRANCH_AFTER_ROOT = auto()
    C5_PRENYL = auto()
    C1_THREE_HALOGEN_LEAVES = auto()


class ArylArmKind(Enum):
    DIRECT_C = auto()
    METHYLENE_C = auto()
    DIRECT_O = auto()
    O_METHYLENE_C = auto()


class HeteroarylKind(Enum):
    SIX_MEMBER_ONE_N = auto()
    FUSED_TEN_MEMBER_C = auto()


class HeteroarylLeaf(Enum):
    HALOGEN = auto()
    METHYL = auto()
    ALKOXY = auto()


@dataclass(frozen=True)
class ArylLeafFact:
    kind: ArylLeafKind
    site: int
    atoms: frozenset[int]
    value: int
    child_ring: frozenset[int]
    child_attach: int
    child_parent: int
    extra_atom: int


@dataclass(frozen=True)
class CarboxyalkylArm:
    attachment: int
    atoms: tuple[int, ...]
    carboxyl: int


# Compatibility alias for Round B1 callers.
CarboxymethylArm = CarboxyalkylArm


@dataclass(frozen=True)
class SidePath:
    root: int
    atoms: tuple[int, ...]


@dataclass(frozen=True)
class ArylArmFact:
    kind: ArylArmKind
    attachment: int
    outer: int
    bridge: int | None
    ring: frozenset[int]
    atoms: tuple[int, ...]


@dataclass(frozen=True)
class AlkoxyFact:
    root: int
    oxygen: int
    code: int
    atoms: tuple[int, ...]


_SHAPE_MATCHERS = {
    AlkylShape.C2_VINYL: side_alkyl._is_vinyl,
    AlkylShape.C3_ALLYL: side_alkyl._is_allyl,
    AlkylShape.C3_ISOPROPENYL: side_alkyl._is_isopropenyl,
    AlkylShape.C3_BRANCH_AT_ROOT: side_alkyl._is_isopropyl,
    AlkylShape.C4_BRANCH_AT_SECOND: side_alkyl._is_sec_butyl,
    AlkylShape.C4_TRIPLE_BRANCH_AT_ROOT: side_alkyl._is_tert_butyl,
    AlkylShape.C4_BRANCH_AFTER_ROOT: side_alkyl._is_isobutyl,
    AlkylShape.C5_ASYMMETRIC_ROOT_BRANCH: side_alkyl._is_2_methylbutan_2_yl,
    AlkylShape.C5_PRENYL: side_alkyl._is_prenyl,
    AlkylShape.C5_BRANCH_NEAR_LEAF: side_alkyl._is_isopentyl,
    AlkylShape.C5_DOUBLE_BRANCH_AFTER_ROOT: side_alkyl._is_neopentyl,
    AlkylShape.C1_THREE_HALOGEN_LEAVES: side_alkyl._is_trifluoromethyl,
}


def _path(root: int, atoms: list[int] | None) -> SidePath | None:
    return SidePath(root, tuple(atoms)) if atoms is not None else None


def _leaf_fact(match: LeafTopology) -> ArylLeafFact:
    return ArylLeafFact(match.kind, match.site, match.atoms, match.value,
                        match.child_ring, match.child_attach,
                        match.child_parent, match.extra_atom)


def ring_leaf(mol: Mol, atom: Atom, ring_atom: int, depth: int) -> ArylLeafFact | None:
    got = match_leaf_topology(mol, atom, ring_atom, depth)
    return _leaf_fact(got) if got is not None else None


def carbon_neighbors(mol: Mol, atom: int) -> list[int]:
    return side_alkyl._c_neighbors(mol, atom)


def alkyl_shape(mol: Mol, root: int, parent: set[int], shape: AlkylShape) -> SidePath | None:
    return _path(root, _SHAPE_MATCHERS[shape](mol, root, parent))


def outer_alkoxy(mol: Mol, root: int, oxygen: int) -> AlkoxyFact | None:
    code = side_alkoxy._outer_alkoxy_n(mol, root, oxygen)
    atoms = side_alkoxy._outer_atoms(mol, root, oxygen, code) if code else None
    return AlkoxyFact(root, oxygen, code, tuple(atoms)) if atoms else None


def phenyl_ring(mol: Mol, root: int, parent: int) -> frozenset[int] | None:
    ring = aryl_sub._phenyl_at(mol, root, parent)
    return frozenset(ring) if ring else None


def benzyl_ring(mol: Mol, root: int, parent: int) -> frozenset[int] | None:
    ring = aryl_sub._ch2_ph_at(mol, root, parent)
    return frozenset(ring) if ring else None


def aryl_leaves(mol: Mol, ring: frozenset[int]) -> frozenset[int]:
    return frozenset(aryl_sub._halo_atoms_on(mol, set(ring)))


def _arm_path(mol: Mol, acid: int, parent: set[int]) -> tuple[int, tuple[int, ...]] | None:
    path, previous, current = [], acid, acid
    while True:
        nexts = [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != previous]
        if len(nexts) != 1: return None
        current = nexts[0]
        if current in parent: return (current, tuple(reversed(path))) if path else None
        path.append(current)
        previous = path[-2] if len(path) > 1 else acid


def _carboxyalkyl_arm(mol: Mol, acid: int, parent: set[int]) -> CarboxyalkylArm | None:
    path = _arm_path(mol, acid, parent)
    return CarboxyalkylArm(path[0], path[1], acid) if path else None


def carboxyalkyl_arms(mol: Mol, chain: list[int], acids: list[int]) -> tuple[CarboxyalkylArm, ...]:
    parent = set(chain)
    facts = (_carboxyalkyl_arm(mol, acid, parent) for acid in acids)
    return tuple(fact for fact in facts if fact is not None)


carboxymethyl_arms = carboxyalkyl_arms
