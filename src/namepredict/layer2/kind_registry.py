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


def _ring_keeps_nh_prefix(mol, chain) -> bool:
    """环内是否含未取代的芳香 NH（决定是否保留 1H- 前缀，P-61.2.4 指示氢）。

    N 全被取代（如 1-甲基咪唑 / 1-苄基吲哚）时 1H- 应省略；mol/chain 缺失时
    保守保留前缀（与既有 indole 硬编码行为一致，不引入新破坏）。
    """
    if mol is None or not chain:
        return True
    for i in chain:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
            return True
    return False


def _strip_locant_prefix(name: str, prefix: str) -> str:
    """去掉词干已带的前导 locant 前缀（'1H-indole' → 'indole'）。"""
    return name[len(prefix):] if prefix and name.startswith(prefix) else name


def pack_parent_stem(parent: dict, mol=None) -> dict:
    """补齐母体词干与编号 scaffold 字段。

    优先按 scaffold_id 查词干（环+FG 收敛为 FG kind 时，词干属 scaffold 而非 FG 类别）。
    五元杂环 locant 前缀（1H- / 1,3-）在此统一成终态：1,3- 二唑无条件注入；1H- 吡咯型
    仅当环含未取代 NH 注入（N 全取代则省略，避免 N-取代环误标 1H-）。词干若已带前缀
    （注册表 indole = "1H-indole"）先剥离再按条件加回，保证 N-取代 indole 输出 "indol-…"。
    """
    from namepredict.layer2.ring_scaffold import locant_prefix

    packed = parent if parent.get("mol") is not None else {**parent, "mol": mol}
    names = (parent_names(packed.get("scaffold_id") or "")
             or parent_names(packed.get("kind") or ""))
    if names is not None and not (packed.get("stem_en") or packed.get("stem_zh")):
        stem_en, stem_zh = names
        pref_en, pref_zh, nh_cond = locant_prefix(packed.get("scaffold_id") or "")
        if pref_en:
            keep = not nh_cond or _ring_keeps_nh_prefix(packed.get("mol"), packed.get("chain"))
            bare_en, bare_zh = _strip_locant_prefix(stem_en, pref_en), _strip_locant_prefix(stem_zh, pref_zh)
            stem_en = pref_en + bare_en if keep else bare_en
            stem_zh = pref_zh + bare_zh if keep else bare_zh
        packed = {**packed, "stem_en": stem_en, "stem_zh": stem_zh}
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
