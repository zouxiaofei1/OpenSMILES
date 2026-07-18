"""Final heavy-atom coverage gate over ownership and named claims."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from rdkit.Chem import Mol

from namepredict.layer3.substituent_namer import SubstituentName


@dataclass(frozen=True)
class CoverageLedger:
    owned_atoms: frozenset[int]
    named_claims: tuple[SubstituentName, ...]
    gap: frozenset[int]
    overlap: frozenset[int]

    @property
    def complete(self) -> bool:
        return not self.gap and not self.overlap


def _heavy_atoms(mol: Mol) -> frozenset[int]:
    return frozenset(a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1)


def _membership_counts(
    owned: frozenset[int], names: list[SubstituentName]
) -> Counter[int]:
    counts: Counter[int] = Counter()
    for idx in owned:
        counts[idx] += 1
    for name in names:
        for idx in name.claim.atoms:
            counts[idx] += 1
    return counts


def _overlap_atoms(counts: Counter[int]) -> frozenset[int]:
    return frozenset(i for i, n in counts.items() if n > 1)


def _covered_atoms(owned: frozenset[int], names: list[SubstituentName]) -> frozenset[int]:
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
    """Build gap/overlap ledger for heavy atoms only; H excluded."""
    heavy = _heavy_atoms(mol)
    covered = _covered_atoms(owned_atoms, names)
    counts = _membership_counts(owned_atoms, names)
    return CoverageLedger(
        owned_atoms=owned_atoms,
        named_claims=tuple(names),
        gap=frozenset(heavy - covered),
        overlap=_overlap_atoms(counts),
    )
