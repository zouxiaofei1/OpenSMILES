"""P-25.3.2 稠合名称组装: fused_info 拆解树生成 benzo[a]/naphtho[...] 稠合 base 名
(未注册稠环; 组分词干/保留前缀由 L2 打包进 FusedNode, 编号走 L4 fused_numbering)。"""
from __future__ import annotations

from collections import Counter
from namepredict.layer1.ring_systems import sssr_rings


_RETAINED_FUSION_ALIASES: dict[str, tuple[str, str]] = {  # 稠合组装名 → 保留名（P-25.1.1）；只收位次形态与稠合名完全相同的条目，命中即整名替换
    "benzo[c]furan":       ("2-benzofuran",  "2-苯并呋喃"),    # 异苯并呋喃：O 占 2 位，1/3 位为环碳
    "benzo[c]pyrrole":     ("isoindole",     "异吲哚"),        # N 占 2 位
    "benzo[b]benzofuran":  ("dibenzofuran",  "二苯并呋喃"),
    "benzo[b]quinoxaline": ("phenazine",     "菲嗪"),
    "benzo[a]indene":      ("fluorene",      "芴"),
    "benzo[d]1,2-oxazole": ("1,2-benzoxazole", "1,2-苯并噁唑"),
}


def _stem_of(node) -> tuple[str | None, str | None]:
    """节点的组分词干(stem_en/stem_zh); L2 未标注的组分(非稠合零件)返回 (None, None)。"""
    return node.fused_stem or (None, None)


def _prefix_of(node, stem_en: str | None, stem_zh: str | None) -> tuple[str, str] | None:
    """附加组分前缀: L2 标注的保留前缀, 否则通用「去尾 e 加 o / 中文加并」(P-25.3.2.2.2)。"""
    if node.fused_prefix:
        return node.fused_prefix
    en = (stem_en or "").rstrip("e") + "o"
    zh = (stem_zh or "") + "并"
    return (en, zh) if en else None


def _component_numbering(mol, node, rings, fusion_edges, shared=None):
    """组分自身编号：委托 L4 fused_component_numbering(P-25.4/P-25.3.3)，把稠合掉的 shared 原子当取代基最小化位次，使取向与规范稠合描述符一致。"""
    rset = sorted(node.ring_indices)
    if not rset:
        return None, None
    sub_rings = [rings[i] for i in rset]
    if shared is None:  # 兜底: 沿用旧逻辑(节点自身的附加组分之一, 无附加时为空)。螺环附加组分只共享 1 个原子、fusion_shared 为空，须与 :79 同样守卫，否则 [0] 越界。
        shared = node.attached[0].fusion_shared[0] if node.attached and node.attached[0].fusion_shared else None
    idx_map = {i: k for k, i in enumerate(rset)}
    sub_edges = [(idx_map[i], idx_map[j], sh) for i, j, sh in fusion_edges
                 if i in idx_map and j in idx_map]
    from namepredict.layer4.numbering_engine import fused_component_numbering
    return fused_component_numbering(mol, node.scaffold_id, sub_rings, shared, sub_edges)


def _inner_atoms(node, rings) -> set[int]:
    """节点内出现在 ≥3 环的原子(perifused 中心, 不在外周边界)。"""
    counts: Counter = Counter()
    for i in node.ring_indices:
        for a in rings[i]:
            counts[a] += 1
    return {a for a, c in counts.items() if c >= 3}


def _outer_chain_labels(node, rings, chain, labels):
    """过滤内原子(≥3 环)后的外周边界 chain/labels。"""
    inner = _inner_atoms(node, rings)
    kept = [(a, lbl) for a, lbl in zip(chain, labels) if a not in inner]
    if not kept:
        return None, None
    return [a for a, _ in kept], [lbl for _, lbl in kept]


def _fusion_letter(parent_chain, shared) -> str | None:
    """共享边在母体外周位次序中的侧字母: 侧(chain[i], chain[i+1]) → chr(97+i)。"""
    shared = [a for a in shared if a in parent_chain]
    if len(shared) != 2:
        return None
    ia, ib = parent_chain.index(shared[0]), parent_chain.index(shared[1])
    i = min(ia, ib)
    if (i + 1) % len(parent_chain) != max(ia, ib):
        return None  # 共享边须为母体外周边
    return chr(97 + i)


def _fusion_numbers(child_chain, child_labels, parent_chain, shared) -> tuple:
    """附加组分共享原子位次(顺序沿母体低位次端→高位次端)。"""
    labels = dict(zip(child_chain, child_labels))
    shared = [a for a in shared if a in labels]
    if len(shared) != 2:
        return ()
    ordered = sorted(shared, key=lambda a: parent_chain.index(a))
    return tuple(labels[a] for a in ordered)


def _fused_one(mol, parent_node, child_node, rings, fusion_edges) -> tuple[str, str] | None:
    """单级: 附加组分前缀 + 融合描述符(数字-字母)。"""
    shared = child_node.fusion_shared[0] if child_node.fusion_shared else None  # 母体与附加组分各自编号都以同一稠合原子集作取代基(P-25.3.1.3: 位次尽可能低)。
    parent_chain, _ = _component_numbering(mol, parent_node, rings, fusion_edges, shared)
    parent_chain, _ = _outer_chain_labels(parent_node, rings, parent_chain, [""] * len(parent_chain)) \
        if parent_chain else (None, None)
    if not parent_chain:
        return None
    prefix = _prefix_of(child_node, *_stem_of(child_node))
    if prefix is None:
        return None
    if child_node.fused_omit_numbers:  # 一级单环烃附加组分(benzo 及 P-25.3.2.2.1 的 cyclopenta 等)省略数字位次(P-25.3.8.1), 故无需附加组分自身编号。
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
        if letter:
            break
    if not letter:
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


def fused_parent_names(mol, node) -> tuple[str, str] | None:
    """稠合名组装入口: node 为 FusedNode 根; 返回 (en, zh) 或 None(无法组装)。"""
    if not node.attached:
        return None  # 单节点保留名走 _parent_stem_names
    root_en, root_zh = _stem_of(node)
    if not root_en:
        return None
    rings = list(sssr_rings(mol))
    from namepredict.layer1.ring_systems import build_ring_systems
    fusion_edges = ()
    root_rings = set(node.ring_indices)
    for s in build_ring_systems(mol):
        if root_rings <= set(s.get("sssr_indices") or ()):
            fusion_edges = s.get("fusion_edges") or ()
            break
    parts = _collect_attached(mol, node, rings, fusion_edges)
    if parts is None:
        return None
    core_en = parts[0] + root_en
    alias = _RETAINED_FUSION_ALIASES.get(core_en)
    return alias if alias else (core_en, parts[1] + root_zh)
