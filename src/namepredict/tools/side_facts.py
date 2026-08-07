"""Layer 2 向 Layer 3 暴露的侧链拓扑事实。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


from rdkit.Chem import Atom, Mol

from namepredict.tools import aryl_sub, heteroaryl_sub, side_alkyl
from namepredict.tools import side_alkoxy, side_cycloalkyl, side_sat_hetero
from namepredict.tools.leaves.protocol import ArylLeafKind, LeafTopology
from namepredict.tools.leaves.registry import match_leaf_topology
from namepredict.tools.leaves.topo import nb_out


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


@dataclass(frozen=True)
class SaturatedHeterocycleFact:
    root: int
    ring: frozenset[int]
    signature: tuple


@dataclass(frozen=True)
class HeteroarylFact:
    kind: HeteroarylKind
    attachment: int
    outer: int
    ring: frozenset[int]
    atoms: tuple[int, ...]
    locant: int
    leaves: tuple[tuple[int, HeteroarylLeaf, int], ...]


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


def _phenyl_arm(item: dict) -> ArylArmFact:
    return ArylArmFact(ArylArmKind.DIRECT_C, item.get("attach"), item.get("outer_c"), None,
                       frozenset(item.get("ph")), tuple(item.get("atoms")))


def _benzyl_arm(item: dict) -> ArylArmFact:
    return ArylArmFact(ArylArmKind.METHYLENE_C, item.get("attach"), item.get("outer_c"), item.get("ch2"),
                       frozenset(item.get("ph")), tuple(item.get("atoms")))


def _phenoxy_arm(item: dict) -> ArylArmFact:
    return ArylArmFact(ArylArmKind.DIRECT_O, item.get("ring_c"), item.get("outer_c"), item.get("o_idx"),
                       frozenset(item.get("ph")), tuple(item.get("atoms")))


def _benzyloxy_arm(item: dict) -> ArylArmFact:
    return ArylArmFact(ArylArmKind.O_METHYLENE_C, item.get("ring_c"), item.get("outer_c"), item.get("ch2"),
                       frozenset(item.get("ph")), tuple(item.get("atoms")))


def ring_outside(mol: Mol, atom: int, ring: set[int]) -> list[Atom]:
    return nb_out(mol, atom, ring)


def _leaf_fact(match: LeafTopology) -> ArylLeafFact:
    return ArylLeafFact(match.kind, match.site, match.atoms, match.value,
                        match.child_ring, match.child_attach,
                        match.child_parent, match.extra_atom)


def ring_leaf(mol: Mol, atom: Atom, ring_atom: int, depth: int) -> ArylLeafFact | None:
    got = match_leaf_topology(mol, atom, ring_atom, depth)
    return _leaf_fact(got) if got is not None else None


def carbon_neighbors(mol: Mol, atom: int) -> list[int]:
    return side_alkyl._c_neighbors(mol, atom)


def linear_alkyl(mol: Mol, root: int, parent: set[int], limit: int = 12) -> SidePath | None:
    return _path(root, side_alkyl._walk_linear_n(mol, root, parent, limit))


def omega_halo_alkyl(mol: Mol, root: int, parent: set[int]) -> SidePath | None:
    return _path(root, side_alkyl._walk_omega_halo(mol, root, parent))


def terminal_halogen(mol: Mol, atom: int) -> int | None:
    return side_alkyl._terminal_halo_z(mol, atom)


def alkyl_shape(mol: Mol, root: int, parent: set[int], shape: AlkylShape) -> SidePath | None:
    return _path(root, _SHAPE_MATCHERS[shape](mol, root, parent))


def outer_alkoxy(mol: Mol, root: int, oxygen: int) -> AlkoxyFact | None:
    code = side_alkoxy._outer_alkoxy_n(mol, root, oxygen)
    atoms = side_alkoxy._outer_atoms(mol, root, oxygen, code) if code else None
    return AlkoxyFact(root, oxygen, code, tuple(atoms)) if atoms else None


def cycloalkyl_side(mol: Mol, root: int, parent: set[int]) -> SidePath | None:
    return _path(root, side_cycloalkyl._is_monocycloalkyl(mol, root, parent))


def cycloalkylethyl_side(mol: Mol, root: int, parent: set[int]) -> SidePath | None:
    return _path(root, side_cycloalkyl._is_1_cycloalkylethyl(mol, root, parent))


def saturated_heterocycle_side(
    mol: Mol, root: int, parent: set[int],
) -> SaturatedHeterocycleFact | None:
    for ring, signature in side_sat_hetero._sat_hetero_rings(mol):
        if root in ring and not ring & parent and side_sat_hetero._exo_only_parent(mol, ring, root, parent):
            return SaturatedHeterocycleFact(root, frozenset(ring), signature)
    return None


def aryl_arms(info: dict, parent: set[int]) -> list[ArylArmFact]:
    mol = info.get("mol")
    phenoxy = map(_phenoxy_arm, aryl_sub._ring_phenoxys(info, parent))
    phenyl = map(_phenyl_arm, aryl_sub._ring_phenyls(mol, parent))
    benzyloxy = map(_benzyloxy_arm, aryl_sub._ring_benzyloxys(info, parent))
    benzyl = map(_benzyl_arm, aryl_sub._ring_benzyls(mol, parent))
    return [*phenoxy, *phenyl, *benzyloxy, *benzyl]


def phenyl_ring(mol: Mol, root: int, parent: int) -> frozenset[int] | None:
    ring = aryl_sub._phenyl_at(mol, root, parent)
    return frozenset(ring) if ring else None


def benzyl_ring(mol: Mol, root: int, parent: int) -> frozenset[int] | None:
    ring = aryl_sub._ch2_ph_at(mol, root, parent)
    return frozenset(ring) if ring else None


def aryl_leaves(mol: Mol, ring: frozenset[int]) -> frozenset[int]:
    return frozenset(aryl_sub._halo_atoms_on(mol, set(ring)))


def heteroaryl_outers(mol: Mol, parent: set[int]) -> set[int]:
    return heteroaryl_sub.heteroaryl_outers(mol, parent)


def _leaf_topology(mol: Mol, ring: set[int], order: list[int], locant: int,
                   outer: int, attach: int) -> tuple[int, HeteroarylLeaf, int]:
    site = order[locant - 1]
    leaves = heteroaryl_sub._nb_out(mol, site, ring)
    leaves = [atom for atom in leaves if not (site == outer and atom.GetIdx() == attach)]
    atom = leaves[0]
    if atom.GetAtomicNum() in (9, 17, 35, 53):
        return locant, HeteroarylLeaf.HALOGEN, atom.GetAtomicNum()
    if atom.GetAtomicNum() == 8:
        return locant, HeteroarylLeaf.ALKOXY, heteroaryl_sub._alkoxy_n(mol, atom, site)
    return locant, HeteroarylLeaf.METHYL, 6


def _pyridinyl_fact(mol: Mol, item: dict) -> HeteroarylFact:
    ring, outer, attach = set(item.get("ring")), item.get("outer_c"), item.get("attach")
    order = heteroaryl_sub._ring_order(mol, ring, outer)
    raw = heteroaryl_sub._leaf_items(mol, ring, outer, attach, order)
    leaves = tuple(_leaf_topology(mol, ring, order, loc, outer, attach) for loc, _, _ in raw)
    return HeteroarylFact(HeteroarylKind.SIX_MEMBER_ONE_N, attach, outer, frozenset(ring),
                          tuple(item.get("atoms")), order.index(outer) + 1, leaves)


def _naphthyl_fact(mol: Mol, item: dict) -> HeteroarylFact:
    ready = heteroaryl_sub._naph_ready(mol, item.get("outer_c"), item.get("attach"))
    locant = heteroaryl_sub._naph_locant(ready[1], item.get("outer_c"))
    ring = frozenset(item.get("atoms"))
    return HeteroarylFact(HeteroarylKind.FUSED_TEN_MEMBER_C, item.get("attach"), item.get("outer_c"), ring,
                          tuple(item.get("atoms")), locant, ())


def naphthyl_facts(mol: Mol, parent: set[int]) -> list[HeteroarylFact]:
    return [_naphthyl_fact(mol, item)
            for item in heteroaryl_sub.ring_naphthyls(mol, parent)]


def pyridinyl_facts(mol: Mol, parent: set[int]) -> list[HeteroarylFact]:
    return [_pyridinyl_fact(mol, item)
            for item in heteroaryl_sub.ring_pyridinyls(mol, parent)]


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
