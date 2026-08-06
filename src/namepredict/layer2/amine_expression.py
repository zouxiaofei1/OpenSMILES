"""Typed parent-arm facts for N-substituted neutral amines."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class NArm:
    atom_ids: tuple[int, ...]
    attachment_atom: int


@dataclass(frozen=True)
class NeutralAmineArmFacts:
    nitrogen_atom: int
    degree: int
    parent_arm: NArm
    substituent_arms: tuple[NArm, ...]


def build_amine_arm_facts(nitrogen: int, degree: int, parent: list[int], rest: list[list[int]]) -> NeutralAmineArmFacts:
    make = lambda arm: NArm(tuple(arm), arm[0])
    return NeutralAmineArmFacts(nitrogen, degree, make(parent), tuple(make(arm) for arm in rest))
