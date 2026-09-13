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
    return packed if facts is None else {
        **packed, "numbering_scaffold": facts, "numbering_scaffold_required": True,
    }


def _ring_keeps_nh_prefix(mol, chain) -> bool:
    """环内是否有未取代的芳香 NH（决定是否保留 1H- 前缀，P-61.2.4）；mol/chain 缺失时保守返回 True。"""
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


def _embeds_locant_prefix(name: str, prefix: str) -> bool:
    """词干是否把 locant 前缀嵌在词中而非词首：前缀须留在组分名前，不再前移或剥离。"""
    return f"-{prefix}" in name


def pack_parent_stem(parent: dict, mol=None) -> dict:
    """补齐母体词干与编号 scaffold 字段（词干属 scaffold 而非 FG 类别）；五元杂环 locant 前缀在此统一成终态：1,3- 二唑无条件注入、1H- 吡咯型仅含未取代 NH 时注入。"""
    from namepredict.layer2.ring_scaffold import locant_prefix

    packed = parent if parent.get("mol") is not None else {**parent, "mol": mol}
    names = (parent_names(packed.get("scaffold_id") or "")
             or parent_names(packed.get("kind") or ""))
    if names is not None and not (packed.get("stem_en") or packed.get("stem_zh")):
        stem_en, stem_zh = names
        pref_en, pref_zh, nh_cond = locant_prefix(packed.get("scaffold_id") or "")
        if pref_en and not _embeds_locant_prefix(stem_en, pref_en):  # 局部不饱和保留名的词干自带加氢前缀（2,5-dihydro-1H-pyrrole），locant 前缀已在其末段，前移会重复位次
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
    """启动时加载 scaffold 词干注册表（唯一事实来源）。"""
    _load_from_scaffold_specs()  # Spec 是词干权威


_bootstrap()
