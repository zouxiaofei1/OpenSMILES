"""对所有权和已命名 claim 的最终重原子覆盖门控。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol

from opensmiles.layer3.substituent_namer import SubstituentName


@dataclass(frozen=True)
class CoverageLedger:
    """所有权与命名 claim 的重原子覆盖台账：缺口与重叠原子集。"""
    owned_atoms: frozenset[int]
    named_claims: tuple[SubstituentName, ...]
    gap: frozenset[int]
    overlap: frozenset[int]

    @property
    def complete(self) -> bool:
        """是否没有缺口也没有重叠。"""
        return not self.gap and not self.overlap


def build_coverage_ledger(
    mol: Mol,
    *,
    owned_atoms: frozenset[int],
    names: list[SubstituentName],
) -> CoverageLedger:
    """仅针对重原子构建 gap/overlap 台账；排除 H。"""
    heavy = frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1)
    counts: Counter[int] = Counter(owned_atoms)  # 归属次数：所有权 +1、每个命名 claim +1
    for name in names:
        counts.update(name.claim.atoms)
    return CoverageLedger(
        owned_atoms=owned_atoms,
        named_claims=tuple(names),
        gap=frozenset(heavy - owned_atoms),
        overlap=frozenset(i for i, n in counts.items() if n > 1),
    )
