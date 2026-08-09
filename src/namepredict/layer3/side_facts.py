"""Layer 3 的侧链拓扑事实：Layer 2 主链选择与 Layer 3 取代基命名共用。"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto


from rdkit.Chem import Atom, Mol

from namepredict.layer3 import side_alkyl
from namepredict.layer3.leaves.protocol import ArylLeafKind


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


def ring_leaf(mol: Mol, atom: Atom, ring_atom: int, depth: int) -> ArylLeafFact | None:
    """Leaf matching retired in 446ce69; kept as typed contract stub."""
    return None


from namepredict.constants import C, F, H, HALO_Z as _HALO_Z
def carbon_neighbors(mol: Mol, atom: int) -> list[int]:
    #return side_alkyl._c_neighbors(mol, atom)
    atom = mol.GetAtomWithIdx(atom)
    return [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == C]



