"""P-25.3.2 稠合名称组装: 用 fused_info 拆解树生成 benzo[a].../naphtho[...]... 类稠合名。

供未注册稠环(母体/附加均为已注册保留组分)的 base 名; 单边融合支持, 多边/跨位放后续。
L5 分层纯净: 组件词干用本地表(与 L2 _TEMPLATES 需同步), 编号走 L4 fused_numbering。
"""
from __future__ import annotations

from collections import Counter

# 附加组分保留前缀(组分 id → 前缀; P-25.3.2.2.3 保留前缀)。
_FUSED_PREFIX = {
    "benzene": ("benzo", "苯并"),
    "naphthalene": ("naphtho", "萘并"),
    "anthracene": ("anthra", "蒽并"),
    "phenanthrene": ("phenanthro", "菲并"),
    "furan": ("furo", "呋喃并"),
    "thiophene": ("thieno", "噻吩并"),
    "pyridine": ("pyrido", "吡啶并"),
    "pyrimidine": ("pyrimido", "嘧啶并"),
    "imidazole": ("imidazo", "咪唑并"),
}

# 组分词干表(L5 本地, 与 layer2.ring_scaffold._TEMPLATES 词干需同步)。
_COMPONENT_STEM = {
    "benzene": ("benzene", "苯"),
    "naphthalene": ("naphthalene", "萘"),
    "anthracene": ("anthracene", "蒽"),
    "phenanthrene": ("phenanthrene", "菲"),
    "pyrene": ("pyrene", "芘"),
    "furan": ("furan", "呋喃"),
    "thiophene": ("thiophene", "噻吩"),
    "pyrrole": ("pyrrole", "吡咯"),
    "pyridine": ("pyridine", "吡啶"),
    "pyrimidine": ("pyrimidine", "嘧啶"),
    "pyrazine": ("pyrazine", "吡嗪"),
    "pyridazine": ("pyridazine", "哒嗪"),
    "imidazole": ("imidazole", "咪唑"),
    "pyrazole": ("pyrazole", "吡唑"),
    "oxazole": ("oxazole", "噁唑"),
    "thiazole": ("thiazole", "噻唑"),
    "quinoline": ("quinoline", "喹啉"),
    "isoquinoline": ("isoquinoline", "异喹啉"),
    "indole": ("indole", "吲哚"),
    "indazole": ("indazole", "吲唑"),
    "benzimidazole": ("benzimidazole", "苯并咪唑"),
    "benzofuran": ("benzofuran", "苯并呋喃"),
    "benzothiophene": ("benzothiophene", "苯并噻吩"),
    "benzothiazole": ("benzothiazole", "苯并噻唑"),
    "benzoxazole": ("benzoxazole", "苯并噁唑"),
    "quinazoline": ("quinazoline", "喹唑啉"),
    "quinoxaline": ("quinoxaline", "喹喔啉"),
    "benzodioxole": ("benzodioxole", "苯并二氧杂环戊烯"),
    "pyrrolidine": ("pyrrolidine", "吡咯烷"),
    "piperidine": ("piperidine", "哌啶"),
    "morpholine": ("morpholine", "吗啉"),
    "piperazine": ("piperazine", "哌嗪"),
    "oxolane": ("oxolane", "四氢呋喃"),
    "oxane": ("oxane", "四氢吡喃"),
    "carbazole": ("carbazole", "咔唑"),
    "acridine": ("acridine", "吖啶"),
    "phenothiazine": ("phenothiazine", "吩噻嗪"),
}


def _stem_of(sid: str) -> tuple[str | None, str | None]:
    """组分的保留词干(stem_en/stem_zh)。"""
    return _COMPONENT_STEM.get(sid, (None, None))


def _prefix_of(sid: str, stem_en: str | None, stem_zh: str | None) -> tuple[str, str] | None:
    """附加组分前缀: 保留前缀表, 否则通用「去尾 e 加 o / 中文加并」(P-25.3.2.2.2)。"""
    p = _FUSED_PREFIX.get(sid)
    if p:
        return p
    en = (stem_en or "").rstrip("e") + "o"
    zh = (stem_zh or "") + "并"
    return (en, zh) if en else None


def _component_numbering(mol, node, rings, fusion_edges):
    """组分自身编号: 委托 L4 fused_component_numbering(P-25.4/P-25.3.3), 稠合点作取代基。"""
    rset = sorted(node.ring_indices)
    if not rset:
        return None, None
    sub_rings = [rings[i] for i in rset]
    shared = node.attached[0].fusion_shared[0] if node.attached else None
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
    parent_chain, _ = _component_numbering(mol, parent_node, rings, fusion_edges)
    parent_chain, _ = _outer_chain_labels(parent_node, rings, parent_chain, [""] * len(parent_chain)) \
        if parent_chain else (None, None)
    child_chain, child_labels = _component_numbering(mol, child_node, rings, fusion_edges)
    # print(child_chain, child_labels )
    if not parent_chain or not child_labels:
        return None
    prefix = _prefix_of(child_node.scaffold_id, *_stem_of(child_node.scaffold_id))
    if prefix is None:
        return None
    letter = numbers = None
    for sh in child_node.fusion_shared:
        letter = _fusion_letter(parent_chain, sh)
        numbers = _fusion_numbers(child_chain, child_labels, parent_chain, sh)
        if letter:
            break
    if not letter:
        return None
    # 单环烃附加(benzo 等)省略数字位次(P-25.3.8.1)。
    desc = f"[{letter}]" if child_node.scaffold_id == "benzene" else f"[{','.join(map(str, numbers))}-{letter}]"
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
    root_en, root_zh = _stem_of(node.scaffold_id)
    if not root_en:
        return None
    rings = list(mol.GetRingInfo().AtomRings())
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
    return parts[0] + root_en, parts[1] + root_zh
