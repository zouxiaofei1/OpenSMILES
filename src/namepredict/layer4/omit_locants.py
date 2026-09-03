"""FG 位次省略规则（L4；P-14.3.4 / 环单 FG）。"""
from __future__ import annotations


def _is_cyclo(parent: dict | None) -> bool:
    """环状单环 carbocycle（scaffold_id 承载环系）；有 fused_tree 的未注册稠环不作单环环烷烃（由 L5 按稠合名处理）。"""
    return (parent or {}).get("scaffold_id") == "carbocycle" and not (parent or {}).get("fused_tree")


def omit_oh(
    oh_pos: int | None, n_carbons: int, kind: str | None = None,
    parent: dict | None = None, n_subs: int = 0, *,
    has_ene=None, has_yne=None,
) -> bool:
    """判定醇羟基位次是否省略（P-14.3.4）。"""
    # 环单醇：有烯保留位次；无取代省略位次（P-14.3.4）。
    if _is_cyclo(parent) and oh_pos is not None:
        if has_ene and has_ene(parent):
            return False
        return n_subs == 0
    if kind == "alcohol" and parent and (
        (has_ene and has_ene(parent)) or (has_yne and has_yne(parent))
    ):
        return False
    return oh_pos == 1 and n_carbons <= 2


def omit_sh(sh_pos: int | None, n_carbons: int) -> bool:
    """判定巯基位次是否省略（短链首位）。"""
    return sh_pos == 1 and n_carbons <= 2


def omit_amine(
    am_pos: int | None, n_carbons: int, kind: str | None = None, n_subs: int = 0,
    parent: dict | None = None,
) -> bool:
    """判定氨基位次是否省略（环单胺或短链）。"""
    # 环单胺：无取代省略位次；有取代保留（cycloamine 规则，按 scaffold 判断）。
    if _is_cyclo(parent) and am_pos is not None:
        return n_subs == 0
    return am_pos == 1 and n_carbons <= 2


def omit_ketone(
    kind: str | None, n_subs: int, parent: dict | None = None, *,
    has_ene=None, single: bool = True,
) -> bool:
    """判定酮位次是否省略（环单酮，按 scaffold 判断）。"""
    # 环单酮：有烯保留位次；无取代省略（cycloketone 规则，按 scaffold 判断）。
    if _is_cyclo(parent) and single:
        if has_ene and has_ene(parent):
            return False
        return n_subs == 0
    return False


def omit_unsat(
    n_carbons: int, kind: str | None = None, parent: dict | None = None, *,
    has_ene=None, has_yne=None,
) -> bool:
    """判定烯/炔位次是否省略（环单烯或短链）。"""
    # 纯烃环单烯（cycloalkene，kind 收敛为 alkane + carbocycle scaffold）位次隐含省略；环多烯（double_bonds）保留位次。
    if kind == "alkane" and (parent or {}).get("scaffold_id") == "carbocycle":
        if not (parent or {}).get("double_bonds"):
            return True
    if kind == "alcohol" and parent and has_yne and has_yne(parent):
        return False
    if parent and has_ene and has_ene(parent) and kind != "alkene":
        return False
    return n_carbons <= 3
