"""ParentKind 元数据注册表 —— 评分 + L5 词干权威（P-44）。"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class KindMeta:
    """母体 kind 元数据：中英词干、环类型、环数与是否保留名。"""
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
    return packed if facts is None else {**packed, "numbering_scaffold": facts}


def _ring_keeps_nh_prefix(mol, chain) -> bool:
    """环内是否有未取代芳香 NH（1H- 前缀判据，P-61.2.4）。"""
    if mol is None or not chain:
        return True
    for i in chain:
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 7 and atom.GetIsAromatic() and atom.GetTotalNumHs() > 0:
            return True
    return False


def _strip_locant_prefix(name: str, prefix: str) -> str:
    """去掉词干已带的前导 locant 前缀。"""
    return name[len(prefix):] if prefix and name.startswith(prefix) else name


def pack_parent_stem(parent: dict, mol=None) -> dict:
    """补齐母体词干与编号 scaffold 字段，并定五元杂环 locant 前缀。"""
    from namepredict.layer2.ring_scaffold import locant_prefix
    from namepredict.layer2.hantzsch_widman import parent_names as hw_parent_names

    packed = parent if parent.get("mol") is not None else {**parent, "mol": mol}
    names = (parent_names(packed.get("scaffold_id") or "")
             or hw_parent_names(packed.get("scaffold_id"), packed.get("mol"), packed.get("chain"))
             or parent_names(packed.get("kind") or ""))
    if names is not None and not (packed.get("stem_en") or packed.get("stem_zh")):
        stem_en, stem_zh = names
        pref_en, pref_zh, nh_cond = locant_prefix(packed.get("scaffold_id") or "")
        if pref_en and f"-{pref_en}" not in stem_en:  # 词干自带加氢前缀，locant 前缀前移会重复位次
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


_load_from_scaffold_specs()  # 启动时加载 scaffold 词干注册表（唯一事实来源）
