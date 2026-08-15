"""ParentKind 元数据注册表 —— 评分 + L5 词干权威（P-44）。"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class KindMeta:
    kind: str
    en: str | None = None
    zh: str | None = None
    ring: str = "none"  # 可选值："none" | "hetero" | "carbo"
    n_rings: int = 0
    retained: bool = False


_REG: dict[str, KindMeta] = {}

def get(kind: str) -> KindMeta | None:
    """按 kind 查 KindMeta（未注册返回 None）。"""
    return _REG.get(kind)


def is_hetero_ring(kind: str) -> int:
    """kind 是否为杂环（1/0）。"""
    m = get(kind)
    return 1 if m and m.ring == "hetero" else 0


def is_carbo_ring(kind: str) -> int:
    """kind 是否为碳环（1/0）。"""
    m = get(kind)
    return 1 if m and m.ring == "carbo" else 0


def n_rings_of(kind: str) -> int:
    """返回 kind 的环数。"""
    m = get(kind)
    return m.n_rings if m else 0


def retained_bonus(kind: str) -> int:
    """kind 是否为保留名（加分 1/0）。"""
    m = get(kind)
    return 1 if m and m.retained else 0


def parent_names(kind: str) -> tuple[str, str] | None:
    """返回 kind 的 (en, zh) 母体名，缺失返回 None。"""
    m = get(kind)
    if m is None or m.en is None or m.zh is None:
        return None
    return m.en, m.zh


def _attach_numbering_scaffold(packed: dict) -> dict:
    """附加编号 scaffold facts（解析失败则原样返回）。"""
    from namepredict.layer2.ring_scaffold import numbering_scaffold_facts
    facts = numbering_scaffold_facts(
        packed.get("scaffold_id") or packed.get("kind"), len(packed.get("chain") or ()),
    )
    return packed if facts is None else {
        **packed, "numbering_scaffold": facts, "numbering_scaffold_required": True,
    }


def pack_parent_stem(parent: dict, mol=None) -> dict:
    """补齐母体词干与编号 scaffold 字段。

    优先按 scaffold_id 查词干（环+FG 收敛为 FG kind 时，词干属 scaffold 而非 FG 类别）。
    """
    packed = parent if parent.get("mol") is not None else {**parent, "mol": mol}
    names = (parent_names(packed.get("scaffold_id") or "")
             or parent_names(packed.get("kind") or ""))
    if names is not None and not (packed.get("stem_en") or packed.get("stem_zh")):
        packed = {**packed, "stem_en": names[0], "stem_zh": names[1]}
    return _attach_numbering_scaffold(packed)

def _load_from_scaffold_specs() -> None:
    """唯一词干权威：ScaffoldSpec → KindMeta（若已存在则覆盖）。"""
    from namepredict.layer2.ring_scaffold import all_specs

    for sp in all_specs():
        if sp.stem_en is None or sp.stem_zh is None:
            continue
        meta = KindMeta( sp.id, sp.stem_en, sp.stem_zh, sp.ring, sp.n_rings, sp.retained,)
        _REG[meta.kind] = meta


def _bootstrap() -> None:
    """启动时加载 scaffold 词干注册表。"""
    _load_from_scaffold_specs()  # 最后：Spec 是词干权威


_bootstrap()
