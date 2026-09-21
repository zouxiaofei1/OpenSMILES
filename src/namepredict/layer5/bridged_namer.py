"""P-23 von Baeyer 桥环母体名：环数词头 + 描述符 + 烷词干 + 骨架置换 'a' 前缀。

同时给出两种词干形态：FG 分支要带 "ane/烷" 的完整名（后缀去 e 用），
链式词干引擎要裸词干（自行拼 "ane/烷" 或 "a…-triene"）。L5 不得 import L2，
因此 node 只按鸭子类型读 descriptor / locant_pairs。
"""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer5.skeleton_replacement import prefix_from_chain
from namepredict.layer5.stems import alkane_en, alkane_zh, stem_forms


def ring_count_prefix(n_rings: int) -> tuple[str, str] | None:
    """环数词头：1 无前缀、2 用 bi/双，3 起用数量词（不是 dicyclo）。"""
    if n_rings == 1:
        return "", ""
    if n_rings == 2:
        return "bi", "双"
    en, zh = MULT_EN.get(n_rings), MULT_ZH.get(n_rings)
    return (en, zh) if en and zh else None


def descriptor_str(descriptor: tuple[int, ...], locant_pairs) -> str:
    """描述符串：次级桥位次对直接续在长度数字后（判分口径无 ^{}）。"""
    head = len(descriptor) - len(locant_pairs)
    parts = [str(x) for x in descriptor[:head]]
    parts += [f"{d}{a},{b}" for d, (a, b) in zip(descriptor[head:], locant_pairs)]
    return ".".join(parts)


def bridged_body_names(node):
    """von Baeyer 主体名（不含 'a' 前缀）→ ((完整英, 完整中), (裸词干英, 裸词干中))。

    P-24.5.2 / P-24.3.4：'a' 前缀须挂在 spiro/spirobi 之前而非组分名内，
    故组分式螺环要取不带 'a' 的主体名，前缀由调用方另拼。
    """
    prefix = ring_count_prefix(len(node.descriptor) - 1)
    if prefix is None:
        return None
    stem_en, stem_zh = alkane_en(sum(node.descriptor) + 2), alkane_zh(sum(node.descriptor) + 2)
    if not stem_en or not stem_zh:
        return None
    desc = descriptor_str(node.descriptor, node.locant_pairs)
    body_en, body_zh = f"{prefix[0]}cyclo[{desc}]", f"{prefix[1]}环[{desc}]"
    return stem_forms(body_en, body_zh, stem_en, stem_zh)


def bridged_parent_names(mol, node, chain: list[int]):
    """桥环名 → ((完整英, 完整中), (裸词干英, 裸词干中))；不可组装返回 None。"""
    body = bridged_body_names(node)
    if body is None:
        return None
    a_en, a_zh = prefix_from_chain(mol, chain)
    if a_en is None:
        return None
    (full_en, full_zh), (bare_en, bare_zh) = body
    return ((f"{a_en}{full_en}", f"{a_zh}{full_zh}"),
            (f"{a_en}{bare_en}", f"{a_zh}{bare_zh}"))
