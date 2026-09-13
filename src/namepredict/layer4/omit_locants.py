"""FG 位次省略规则（L4；P-14.3.4 / 环单 FG）。"""
from __future__ import annotations


def _is_cyclo(parent: dict | None) -> bool:
    """环状单环 carbocycle（scaffold_id 承载环系）；有 fused_tree 的未注册稠环不作单环环烷烃（由 L5 按稠合名处理）。"""
    return (parent or {}).get("scaffold_id") == "carbocycle" and not (parent or {}).get("fused_tree")


def omit_fg_locant(
    pos: int | None, n_carbons: int, parent: dict | None = None, n_subs: int = 0, *,
    single: bool = True, has_ene=None, has_yne=None,
) -> bool:
    """判定主官能团位次是否省略（P-14.3.4 / 环单 FG）：母体另有烯/炔必留位次；环状无取代省位次；开链仅 C1–C2 首位省。"""
    if parent and ((has_ene and has_ene(parent)) or (has_yne and has_yne(parent))):
        return False
    if _is_cyclo(parent) and pos is not None and single:  # 环单 FG：无取代省位次，有取代保留；多官能团（single=False）一律保留
        return n_subs == 0
    return pos == 1 and n_carbons <= 2


def omit_unsat(
    n_carbons: int, kind: str | None = None, parent: dict | None = None, *,
    has_ene=None, has_yne=None, triple: bool = False,
) -> bool:
    """判定烯/炔位次是否省略（环单烯或短链）；triple 选择炔规则。"""
    if kind == "alkane" and (parent or {}).get("scaffold_id") == "carbocycle":  # 纯烃环单烯（cycloalkene，kind 收敛为 alkane + carbocycle scaffold）位次隐含省略；环多烯（double_bonds）保留位次
        if not (parent or {}).get("double_bonds"):
            return True
    if kind == "alcohol" and parent and has_yne and has_yne(parent):
        return False
    if parent and has_ene and has_ene(parent) and kind not in ("alkene", "alkane"):  # 开链烃骨架（kind 为 alkane）不因自身双键而保留位次，交回短链规则
        return False
    # 未取代二核烯（ethene，P-14.3.4.2(d)）与二/三核炔（ethyne/propyne）位次 '1' 隐含省略；C3 烯（丙-1-烯）仍保留。
    return n_carbons <= (3 if triple else 2)
