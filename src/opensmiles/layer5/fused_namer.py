"""P-25.3.2 稠合名组装: 由 fused_info 树生成 base 名。
(未注册稠环; 组分词干由 L2 打包, 编号走 L4)。"""
from __future__ import annotations

import re
from collections import Counter
from opensmiles.constants import RETAINED_FUSION_ALIASES
from opensmiles.layer1.ring_systems import build_ring_systems, sssr_rings
from opensmiles.layer4.numbering_engine import fused_component_numbering

_LEAD_LOCANT_RUN_RE = re.compile(r"^(\d+(?:,\d+)*[a-z]?)-(?=[^\d])")


def _bracket_lead_locants(stem: str) -> str:
    """组分名前导位次串加方括号（P-25.3.1.3）。"""
    m = _LEAD_LOCANT_RUN_RE.match(stem or "")
    return f"[{m.group(1)}]{stem[m.end():]}" if m else stem


def _component_prefix(node) -> tuple[str, str] | None:
    """附加组分前缀: 保留前缀, 或「去尾 e 加 o」（P-25.3.2.2.2）。"""
    if node.fused_prefix:
        return node.fused_prefix
    stem_en, stem_zh = node.fused_stem or (None, None)
    en = (stem_en or "").rstrip("e") + "o"
    zh = (stem_zh or "") + "并"
    return (en, zh) if en else None


def _component_numbering(mol, node, rings, fusion_edges, shared=None, *, side_letter=False):
    """组分自身编号：委托 L4 按 shared 原子取向（P-25.4）。"""
    rset = sorted(node.ring_indices)
    if not rset:
        return None, None
    sub_rings = [rings[i] for i in rset]
    if shared is None:
        shared = node.attached[0].fusion_shared[0] if node.attached and node.attached[0].fusion_shared else None
    idx_map = {i: k for k, i in enumerate(rset)}
    sub_edges = [(idx_map[i], idx_map[j], sh) for i, j, sh in fusion_edges
                 if i in idx_map and j in idx_map]
    return fused_component_numbering(mol, node.scaffold_id, sub_rings, shared, sub_edges,
                                     side_letter=side_letter)


def _inner_atoms(node, rings) -> set[int]:
    """节点内出现在 ≥3 环的原子(perifused 中心, 不在外周边界)。"""
    counts: Counter = Counter()
    for i in node.ring_indices:
        for a in rings[i]:
            counts[a] += 1
    return {a for a, c in counts.items() if c >= 3}


def _fusion_letter(parent_chain, shared) -> str | None:
    """共享边在母体外周位次序中的侧字母。"""
    shared = [a for a in shared if a in parent_chain]
    if len(shared) != 2:
        return None
    ia, ib = parent_chain.index(shared[0]), parent_chain.index(shared[1])
    n = len(parent_chain)
    i, j = min(ia, ib), max(ia, ib)
    if (i + 1) % n == j:
        return chr(97 + i)
    if i == 0 and j == n - 1:  # P-25.3.1.3：收尾侧 (n,1) 记末位字母，勿判为畸变
        return chr(97 + n - 1)
    return None  # 共享边须为母体外周边


def _fusion_numbers(child_chain, child_labels, parent_chain, shared) -> tuple | None:
    """附加组分共享原子位次(顺序沿母体低位次端→高位次端)。"""
    labels = dict(zip(child_chain, child_labels))
    shared = [a for a in shared if a in labels]
    if len(shared) != 2:
        return None  # 位次缺失即判失败，交由上层回退，勿吐出 [-字母] 畸形名
    ordered = sorted(shared, key=lambda a: parent_chain.index(a))
    return tuple(labels[a] for a in ordered)


def _fused_one(mol, parent_node, child_node, rings, fusion_edges) -> tuple[str, str] | None:
    """单级: 附加组分前缀 + 融合描述符(数字-字母)。"""
    shared = child_node.fusion_shared[0] if child_node.fusion_shared else None  # 双方编号同以稠合原子集作取代基（P-25.3.1.3）
    parent_chain, _ = _component_numbering(mol, parent_node, rings, fusion_edges, shared,
                                           side_letter=True)
    if parent_chain:  # 母体外周边界只需保留外侧原子：剔除出现在 ≥3 环的 perifused 中心
        inner = _inner_atoms(parent_node, rings)
        parent_chain = [a for a in parent_chain if a not in inner]
    if not parent_chain:
        return None
    prefix = _component_prefix(child_node)
    if prefix is None:
        return None
    if child_node.fused_omit_numbers:  # 一级单环烃附加组分省略数字位次（P-25.3.8.1）
        letter = next((ltr for sh in child_node.fusion_shared
                       if (ltr := _fusion_letter(parent_chain, sh))), None)
        if not letter:
            return None
        return prefix[0] + f"[{letter}]", prefix[1] + f"[{letter}]"
    child_chain, child_labels = _component_numbering(mol, child_node, rings, fusion_edges, shared)
    if not child_labels:
        return None
    letter = numbers = None
    for sh in child_node.fusion_shared:
        letter = _fusion_letter(parent_chain, sh)
        numbers = _fusion_numbers(child_chain, child_labels, parent_chain, sh)
        if letter and numbers:
            break
    if not letter or not numbers:
        return None
    desc = f"[{','.join(map(str, numbers))}-{letter}]"
    return prefix[0] + desc, prefix[1] + desc


def _collect_attached(mol, parent_node, rings, fusion_edges) -> tuple[str, str] | None:
    """递归收集 parent_node 的全部附加组分前缀(二级在附着的一级前)。"""
    en_parts, zh_parts = [], []
    for child in sorted(parent_node.attached, key=lambda n: n.scaffold_id):
        sub = _collect_attached(mol, child, rings, fusion_edges)
        pre = _fused_one(mol, parent_node, child, rings, fusion_edges)
        if pre is None:
            return None
        en_parts.append((sub[0] if sub else "") + pre[0])
        zh_parts.append((sub[1] if sub else "") + pre[1])
    return "".join(en_parts), "".join(zh_parts)


def _supported_fusion_tree(node) -> bool:
    """只实现双组分稠合名（P-25.3.2）：一个母体 + 一个一级附加组分。"""
    return len(node.attached) == 1 and not node.attached[0].attached


def fused_parent_names(mol, node) -> tuple[str, str] | None:
    """稠合名组装入口: node 为 FusedNode 根。"""
    if not node.attached or not _supported_fusion_tree(node):
        return None  # 单节点保留名走 _parent_stem_names；多组分超出双组分式能力
    root_en, root_zh = node.fused_stem or (None, None)
    if not root_en:
        return None
    rings = list(sssr_rings(mol))
    fusion_edges = ()
    root_rings = set(node.ring_indices)
    for s in build_ring_systems(mol):
        if root_rings <= set(s.get("sssr_indices") or ()):
            fusion_edges = s.get("fusion_edges") or ()
            break
    parts = _collect_attached(mol, node, rings, fusion_edges)
    if parts is None:
        return None
    alias = RETAINED_FUSION_ALIASES.get(parts[0] + root_en)  # 保留名按未括形态匹配
    if alias:
        return alias
    return (parts[0] + _bracket_lead_locants(root_en),
            parts[1] + _bracket_lead_locants(root_zh))
