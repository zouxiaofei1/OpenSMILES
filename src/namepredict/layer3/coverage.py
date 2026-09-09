"""对所有权和已命名 claim 的最终重原子覆盖门控。"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol

from namepredict.layer3.substituent_namer import SubstituentName


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


def _heavy_atoms(mol: Mol) -> frozenset[int]:
    """收集所有非氢原子的索引集合。"""
    return frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1)


def _membership_counts(
    owned: frozenset[int], names: list[SubstituentName]
) -> Counter[int]:
    """统计所有权与命名 claim 的原子归属次数。"""
    counts: Counter[int] = Counter()
    for idx in owned:
        counts[idx] += 1
    for name in names:
        for idx in name.claim.atoms:
            counts[idx] += 1
    return counts


def _overlap_atoms(counts: Counter[int]) -> frozenset[int]:
    """取出归属次数大于 1 的重叠原子。"""
    return frozenset(i for i, n in counts.items() if n > 1)


def _covered_atoms(owned: frozenset[int], names: list[SubstituentName]) -> frozenset[int]:
    """合并所有权与命名 claim 的覆盖原子集合。"""
    covered = set(owned)
    for name in names:
        covered |= set(name.claim.atoms)
    return frozenset(covered)


def build_coverage_ledger(
    mol: Mol,
    *,
    owned_atoms: frozenset[int],
    names: list[SubstituentName],
) -> CoverageLedger:
    """仅针对重原子构建 gap/overlap 台账；排除 H。"""
    heavy = _heavy_atoms(mol)
    covered = _covered_atoms(owned_atoms, names)
    counts = _membership_counts(owned_atoms, names)
    return CoverageLedger(
        owned_atoms=owned_atoms,
        named_claims=tuple(names),
        gap=frozenset(heavy - covered),
        overlap=_overlap_atoms(counts),
    )
