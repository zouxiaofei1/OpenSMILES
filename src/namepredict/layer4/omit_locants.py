"""FG locant omit rules (L4; P-14.3.4 / cyclo mono FG)."""
from __future__ import annotations


def _omit_cyclo_fg(kind: str | None, target: str, n_subs: int) -> bool:
    """Unsubstituted cyclo FG omits locant; any ring sub keeps FG@1."""
    return kind == target and n_subs == 0


def _keep_cyclo_ene_fg(kind: str | None, target: str, parent, has_ene) -> bool:
    return kind == target and bool(parent and has_ene and has_ene(parent))


def omit_oh(
    oh_pos: int | None, n_carbons: int, kind: str | None = None,
    parent: dict | None = None, n_subs: int = 0, *,
    has_ene=None, has_yne=None,
) -> bool:
    if _keep_cyclo_ene_fg(kind, "cycloalcohol", parent, has_ene):
        return False
    if kind == "cycloalcohol":
        return n_subs == 0
    if kind == "alcohol" and parent and (
        (has_ene and has_ene(parent)) or (has_yne and has_yne(parent))
    ):
        return False
    return oh_pos == 1 and n_carbons <= 2


def omit_sh(sh_pos: int | None, n_carbons: int) -> bool:
    return sh_pos == 1 and n_carbons <= 2


def omit_amine(
    am_pos: int | None, n_carbons: int, kind: str | None = None, n_subs: int = 0,
) -> bool:
    if _omit_cyclo_fg(kind, "cycloamine", n_subs):
        return True
    if kind == "cycloamine":
        return False
    return am_pos == 1 and n_carbons <= 2


def omit_ketone(
    kind: str | None, n_subs: int, parent: dict | None = None, *, has_ene=None,
) -> bool:
    if _keep_cyclo_ene_fg(kind, "cycloketone", parent, has_ene):
        return False
    return _omit_cyclo_fg(kind, "cycloketone", n_subs)


def omit_unsat(
    n_carbons: int, kind: str | None = None, parent: dict | None = None, *,
    has_ene=None, has_yne=None,
) -> bool:
    if kind == "cycloalkene":
        return True
    if kind in ("cycloalcohol", "cycloketone") and parent and has_ene and has_ene(parent):
        return False
    if kind == "alcohol" and parent and has_yne and has_yne(parent):
        return False
    if parent and has_ene and has_ene(parent) and kind != "alkene":
        return False
    return n_carbons <= 3
