"""submol 切割 → free-name → -yl（P-29）。
无宿主接线，纯 cut→pipeline→yl。
"""
from __future__ import annotations

import copy
import re

from opensmiles.constants import (
    AMIDO_RETAINED_EN, DIATOMIC_BRIDGE_YL, SIMPLE_ALKOXY_NO_PAREN, SIMPLE_BRIDGE_YL_NO_PAREN,
    retained_dehydro_yl,
)
from opensmiles.tools.common_names import CommonNameCache
from opensmiles.layer3.submol_build import build_anchor_submol


def _fix_rs_with_real(root_mol, block_root_order: list[int], anchored, hit):
    """用原始根分子 CIP 校正取代基 R/S（异头位被 * 顶替后会翻转）。"""
    if not (hit.success and hit.en):
        return hit
    from opensmiles.layer5.stereo import _cip_on_chain, _with_rs
    parent = hit.meta or {}
    if parent.get("parent_kind") != "radical":
        return hit
    chain = parent.get("parent_chain") or []
    if not chain or any(i >= len(block_root_order) for i in chain):  # 空链或含 * 的链：无法映射回根分子，跳过。
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
    single = len(chain) == 1  # 单原子链即唯一位次，省略位次号（同 stereo._rs_parts）
    rs_lab = [(None if single else (labels[pos - 1] if 0 < pos <= len(labels) else pos), code)  # 位次改用整体编号标签（稠环桥头 4aS/8aS）
              for pos, code in rs_real]  # 链序号 pos（1 起）→ 整体标签；无标签表或越界时退回 pos

    out = copy.copy(hit)
    out.en = _with_rs(hit.en, rs_lab)
    out.zh = _with_rs(hit.zh, rs_lab)
    return out


_P_ACYL_STEM_TAIL = ("phosphoryl", "phosphanyl", "phosphinothioyl")  # P 酰基词干名尾（P-67.1.4.1.1）

_AMIDINE_N_AMINO_EN = re.compile(r"^diaminomethylidene\((.+)\)amino$")
_AMIDINE_N_AMINO_ZH = re.compile(r"^二氨基亚甲基\((.+)\)氨基$")


def _carbamimidoyl_prefix_fixup(en: str, zh: str) -> tuple[str, str]:
    """N-取代脒自由名改写为 carbamimidoyl(R)amino。"""
    m = _AMIDINE_N_AMINO_EN.match(en or "")
    if m is None:
        return en, zh
    mz = _AMIDINE_N_AMINO_ZH.match(zh or "")
    return (f"carbamimidoyl({m.group(1)})amino",
            f"氨基甲亚氨酰基({mz.group(1)})氨基" if mz else zh)


def _obridge_front_simple(mol, atoms, attach_old, *, cache, root_ctx):
    """O/S 桥前端是否为简单取代基；无法判定时返回 None。"""
    a = mol.GetAtomWithIdx(attach_old)
    if a.GetAtomicNum() not in (8, 16) or a.GetDegree() != 2:
        return None
    ins = [n.GetIdx() for n in a.GetNeighbors()
           if n.GetAtomicNum() != 1 and n.GetIdx() in atoms]
    if len(ins) != 1:
        return None
    from opensmiles.layer3.claimable_block import ClaimedBlock, SideSlot
    from opensmiles.layer3.substituent_namer import SubstituentNamer
    claim = ClaimedBlock(slot=SideSlot.OTHER, attach_parent=attach_old, root=ins[0],
                         atoms=frozenset(atoms) - {attach_old})
    named = SubstituentNamer(cache=cache, root_ctx=root_ctx).name(mol, claim)
    if named is None:
        return None
    en = named.en or ""
    if en.endswith(_P_ACYL_STEM_TAIL) and not re.search(r"[-\-()\[\]]", en):  # noqa: RUF001
        return True
    return not named.requires_parentheses


def _radical_yl_from_sub(
    mol, atoms: frozenset, attach_old: int, *,
    cache: CommonNameCache | None, root_ctx: tuple | None = None,
) -> tuple[str, str, bool] | None:
    """碳连接点：锚定 * 走 radical 管线，L4 位次 + P-22.2.4。"""
    from opensmiles.namer import _cache_put, _name_mol
    from rdkit import Chem
    anchored = build_anchor_submol(mol, atoms, attach_old)
    if anchored is None:
        return None
    root_mol, to_root = root_ctx if root_ctx is not None else (mol, None)  # 根分子上下文：块原子→根分子索引，供 R/S 在完整分子上重算
    order = sorted(atoms)
    block_root_order = order if to_root is None else [to_root[o] for o in order]
    anchored_to_root = block_root_order + [-1]          # 锚定子结构按 order 复制 + 末尾 dummy
    smiles = Chem.MolToSmiles(anchored)
    hit = cache.get(smiles) if cache is not None else None
    if hit is None:
        hit = _name_mol(anchored, cache=cache,
                        root_ctx=(root_mol, anchored_to_root))
        if cache is not None and hit.success and hit.en:
            _cache_put(cache, smiles, copy.copy(hit))  # 只缓存片段自身自由基名；立体须按当前根重算，不能跨根共享
    hit = _fix_rs_with_real(root_mol, block_root_order, anchored, hit)  # R/S 取决于宿主根分子：fresh 与 cache 命中都校正一次
    if not hit.success or not hit.en:
        return None
    composite = int((hit.meta or {}).get("parent_substituent_count") or 0) > 0
    need_paren = composite and hit.en not in (  # P-16.5.1.1 复合前缀必括；amido 类免括
        "phenyl", *SIMPLE_ALKOXY_NO_PAREN, *AMIDO_RETAINED_EN, *SIMPLE_BRIDGE_YL_NO_PAREN)
    if (hit.meta or {}).get("bridge_self_enclosed"):  # S 桥复合前端名已自含围栏，L5 不得再整体加括号
        need_paren = False
    if need_paren and hit.en.endswith(("oxy", "sulfanyl")) and not hit.en.endswith(DIATOMIC_BRIDGE_YL):
        if  _obridge_front_simple(mol, atoms, attach_old, cache=cache, root_ctx=root_ctx):
            need_paren = False
    en, zh = retained_dehydro_yl(hit.en, hit.zh)  # 保留名取代基改去氢前缀（adamantan-2-yl）
    en, zh = _carbamimidoyl_prefix_fixup(en, zh)  # N-取代脒片段改 carbamimidoyl 形式
    return en, zh, need_paren


def name_as_substituent(
    mol, attach_old: int, atoms, *,
    cache: CommonNameCache | None = None, root_ctx: tuple | None = None,
) -> tuple[str, str, bool] | None:
    """在 attach_old 处切割，free-name 后输出 -yl 双语名。"""
    atoms = frozenset(atoms)
    if mol.GetAtomWithIdx(attach_old).GetAtomicNum() <= 1:
        return None  # 非重原子连接点（dummy/氢）：无 -yl 名，dummy 会无限递归
    return _radical_yl_from_sub(mol, atoms, attach_old,
                                cache=cache, root_ctx=root_ctx)
