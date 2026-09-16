"""submol 切割 → free-name → -yl（P-29）。
无宿主接线，纯 cut→pipeline→yl。
"""
from __future__ import annotations

import copy

from namepredict.constants import (
    AMIDO_RETAINED_EN, DIATOMIC_BRIDGE_YL, SIMPLE_ALKOXY_NO_PAREN, SIMPLE_BRIDGE_YL_NO_PAREN,
)
from namepredict.tools.common_names import CommonNameCache
from namepredict.layer3.submol_build import build_anchor_submol


def _fix_rs_with_real(root_mol, block_root_order: list[int], anchored, hit):
    """用原始根分子 CIP 校正取代基 R/S（异头位被 * 顶替后会翻转）。"""
    if not (hit.success and hit.en):
        return hit
    from namepredict.layer5.stereo import _cip_on_chain, _with_rs
    parent = hit.meta or {}
    if parent.get("parent_kind") != "radical":
        return hit
    chain = parent.get("parent_chain") or []
    if len(chain) < 2 or any(i >= len(block_root_order) for i in chain):  # [0]-collapsed 或含 * 的链：本层未产出 R/S，跳过。
        return hit
    real = [block_root_order[i] for i in chain]
    rs_anch = _cip_on_chain(anchored, chain)      # * 锚定算出的（可能错误）R/S
    if not rs_anch:
        return hit                      # 本层未贡献 R/S（由内层 -yl 携带），不重复处理
    try:
        rs_real = _cip_on_chain(root_mol, real)   # 原始根分子 CIP
    except Exception:
        return hit
    if rs_real == rs_anch:
        return hit
    labels = parent.get("parent_labels") or []
    rs_lab = [(labels[pos - 1] if 0 < pos <= len(labels) else pos, code)  # 位次改用整体编号标签（稠环桥头 4aS/8aS），与 _chain_locant 同约定
              for pos, code in rs_real]  # 链序号 pos（1 起）→ 整体标签；无标签表或越界时退回 pos

    out = copy.copy(hit)
    out.en = _with_rs(hit.en, rs_lab)
    out.zh = _with_rs(hit.zh, rs_lab)
    return out


def _obridge_front_simple(mol, atoms, attach_old, *, cache, root_ctx):
    """O/S 桥前端是否为简单取代基；无法判定时返回 None。"""
    a = mol.GetAtomWithIdx(attach_old)
    if a.GetAtomicNum() not in (8, 16) or a.GetDegree() != 2:
        return None
    ins = [n.GetIdx() for n in a.GetNeighbors()
           if n.GetAtomicNum() != 1 and n.GetIdx() in atoms]
    if len(ins) != 1:
        return None
    front = frozenset(atoms) - {attach_old}
    f = ins[0]
    from namepredict.tools.anchored_table import anchored_lookup
    ret = anchored_lookup(mol, front, f)
    if ret is not None:
        return not ret[2]
    fr = name_as_substituent(mol, f, front, cache=cache, root_ctx=root_ctx)
    return None if fr is None else (not fr[2])


def _radical_yl_from_sub(
    mol, atoms: frozenset, attach_old: int, *,
    cache: CommonNameCache | None, root_ctx: tuple | None = None,
) -> tuple[str, str, bool] | None:
    """碳连接点：锚定 * 走 radical 管线，L4 位次 + P-22.2.4。"""
    from namepredict.namer import _cache_put, _name_mol
    from rdkit import Chem
    anchored = build_anchor_submol(mol, atoms, attach_old)
    if anchored is None:
        return None
    root_mol, to_root = root_ctx if root_ctx is not None else (mol, None)  # 根分子上下文：块原子→原始根分子索引，供 R/S 在完整分子上重算（糖苷异头碳 CIP 随配基翻转，切断碎片会算反）。
    order = sorted(atoms)
    block_root_order = order if to_root is None else [to_root[o] for o in order]
    anchored_to_root = block_root_order + [-1]          # 锚定子结构按 order 复制 + 末尾 dummy
    smiles = Chem.MolToSmiles(anchored)
    hit = cache.get(smiles) if cache is not None else None
    if hit is None:
        hit = _name_mol(anchored, cache=cache,
                        root_ctx=(root_mol, anchored_to_root))
        if cache is not None and hit.success and hit.en:
            _cache_put(cache, smiles, copy.copy(hit))  # 只缓存片段自身自由基名（保留锚定链、不做宿主校正）；立体随宿主根变化，须按当前根重算，不能跨根共享。
    hit = _fix_rs_with_real(root_mol, block_root_order, anchored, hit)  # R/S 取决于宿主根分子：fresh 与 cache 命中都按当前根分子校正一次。
    if not hit.success or not hit.en:
        return None
    composite = int((hit.meta or {}).get("parent_substituent_count") or 0) > 0
    # P-16.5.1.1 复合前缀必括；amido(P-66.1.1.4.3)免括
    need_paren = composite and hit.en not in (
        "phenyl", *SIMPLE_ALKOXY_NO_PAREN, *AMIDO_RETAINED_EN, *SIMPLE_BRIDGE_YL_NO_PAREN)
    if (hit.meta or {}).get("bridge_self_enclosed"):  # S 桥复合前端名已自含围栏（(4-甲氧基苯基)磺酰基），L5 不得再整体加括号
        need_paren = False
    # P-63.2.1/.2.2 简单前端 + O/S 桥不加围栏；双原子桥（…disulfanyl）围栏改由 L5 按前端定形
    if need_paren and hit.en.endswith(("oxy", "sulfanyl")) and not hit.en.endswith(DIATOMIC_BRIDGE_YL):
        if  _obridge_front_simple(mol, atoms, attach_old, cache=cache, root_ctx=root_ctx):
            need_paren = False
    return hit.en, hit.zh, need_paren


def name_as_substituent(
    mol, attach_old: int, atoms, *,
    cache: CommonNameCache | None = None, root_ctx: tuple | None = None,
) -> tuple[str, str, bool] | None:
    """在 attach_old 处切割，free-name 后输出 -yl 双语名。"""
    atoms = frozenset(atoms)
    if mol.GetAtomWithIdx(attach_old).GetAtomicNum() <= 1:
        return None  # 连接点非重原子（dummy/氢）：无 -yl 名；dummy 锚点会自复制致无限递归
    return _radical_yl_from_sub(mol, atoms, attach_old,  # print(_radical_yl_from_sub(...)) 调试用
                                cache=cache, root_ctx=root_ctx)
