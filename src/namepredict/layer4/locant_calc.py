"""L4 位次计算：将官能团/取代基附着原子映射为链上位次。"""
from __future__ import annotations
from namepredict.layer1.fg_registry import FG_SPECS
import re
from namepredict.layer4._chain_orient import _bond_min_locs
from namepredict.layer4.omit_locants import omit_fg_locant as _omit_fg

def locant_key(x) -> tuple[int, str]:
    """locant → 排序键：数字按数值，字母尾作次级键。"""
    m = re.match(r"(\d+)([a-z]*)", str(x))
    return (int(m.group(1)), m.group(2)) if m else (0, "")


def locant_str_sort(locs) -> list:
    """按 locant_key 排序 locant 集合(兼容 int 与字母位)。"""
    return sorted(locs, key=locant_key)


def _typed_group_atoms(parent: dict, group: str) -> list[int]:
    """返回 parent 中指定基团类型的附着原子。"""
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts and facts.group_class.value == group else []

def _atom_locant(chain: list[int], atom: int | None, facts=None) -> int | str | None:
    """按固定标签定位次(数字 int / 字母位 "4a" 原样)，否则用链编号。"""
    if atom is None or atom not in chain:
        return None
    labels = (facts or {}).get("labels")
    if labels and len(labels) == len(chain):
        lbl = labels[chain.index(atom)]
        return int(lbl) if str(lbl).isdigit() else str(lbl)
    return chain.index(atom) + 1

def _atom_locants(oriented: dict, atoms) -> list[int]:
    """把一组骨架原子映射为位次列表（保留 fused 字母位，跳过不在链表中的原子）。"""
    chain, facts = oriented.get("chain") or [], oriented.get("numbering_scaffold")
    return [loc for atom in atoms
            if (loc := _atom_locant(chain, atom, facts)) is not None]


def _typed_atom_locants(oriented: dict, group: str) -> list[int]:
    """返回指定基团全部附着原子在链上的位次列表。"""
    return locant_str_sort(_atom_locants(oriented, _typed_group_atoms(oriented, group)))


def _single_locant(oriented: dict, group: str) -> int | None:
    """主官能团唯一位次；多原子或缺失返回 None。"""
    locs = _typed_atom_locants(oriented, group)
    return locs[0] if len(locs) == 1 else None



def _omit_unsat(n_carbons, kind=None, parent=None, triple=False):
    """委托 omit_locants.omit_unsat 判定不饱和位次省略。"""
    from namepredict.layer4.omit_locants import omit_unsat as _core
    return _core(n_carbons, kind, parent, has_ene=None, has_yne=None,
                 triple=triple)

def _sub_locant(chain: list[int], attach: int, facts=None) -> int:
    """返回取代基附着原子的位次；无位次时取 0。"""
    loc = _atom_locant(chain, attach, facts)
    return 0 if loc is None else loc

def _with_locants(chain: list[int], substituents: list, facts=None) -> list:
    """为每个取代基附加其 locant 后返回新列表。"""
    return [{**s, "locant": _sub_locant(chain, s["attach_idx"], facts)} for s in substituents]

_UNSAT_BOND_KEY = {"ene": "double", "yne": "triple"}  # 不饱和键类别 → 父字典字段前缀。

def _unsat_bonds(oriented: dict, b: str) -> list | None:
    """取某类不饱和键的边列表，单键标量与多键列表统一成列表。"""
    key = _UNSAT_BOND_KEY[b]
    scalar = oriented.get(f"{key}_bond")
    return [scalar] if scalar else oriented.get(f"{key}_bonds") or None


def _bond_locants(oriented: dict, b: str) -> list[int] | None:
    """返回某类不饱和键较小端点位次的排序列表；缺失返回 None。"""
    bonds = _unsat_bonds(oriented, b)
    if not bonds:
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None


def ene_locants(oriented: dict) -> list[int] | None:
    """返回全部双键端点较小位次的排序列表（单烯亦为单元素列表）。"""
    return _bond_locants(oriented, "ene")


def yne_locants(oriented: dict) -> list[int] | None:
    """返回全部三键端点较小位次的排序列表（单炔亦为单元素列表）。"""
    return _bond_locants(oriented, "yne")

def _unsat_locants(oriented: dict, n: int) -> dict:
    """打包烯/炔位次及其省略标志。"""
    kind = oriented.get("kind")
    return {
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": _omit_unsat(n, kind, oriented),
        "yne_locants": yne_locants(oriented),
        "omit_yne_locant": _omit_unsat(n, kind, oriented, triple=True),
    }

_FG_GROUP = {"oh": "alcohol", "amine": "amine", "ketone": "ketone", "sh": "thiol"}  # 记录 kind → principal_expression_facts 类别


def _omit_for(kind: str, oriented: dict, n: int, n_subs: int) -> bool:
    """FG 记录 omit 标志：环状判断交由 omit_locants 完成。"""
    group = _FG_GROUP.get(kind)
    if group is None:
        return False
    single = kind != "ketone" or len(_typed_group_atoms(oriented, "ketone")) == 1
    return _omit_fg(_single_locant(oriented, group), n, oriented, n_subs, single=single,
                    has_ene=None, has_yne=None)


def _anchor_field_locants(oriented: dict, key: str) -> list[int] | None:
    """由固定 locant 1 锚点字段（radical_c_idx）算位次列表。"""
    value = oriented.get(key)
    if value is None:
        return None
    return locant_str_sort(_atom_locants(oriented, [value])) or None


def _exocyclic_only(oriented: dict) -> bool:
    """principal 基团是否以环外方式表达（位次落在环附着原子）。"""
    facts = oriented.get("principal_expression_facts")
    return bool(facts) and facts.relation.value == "exocyclic"


def _locants_for(oriented: dict, spec) -> list[int]:
    """取 spec 对应主官能团的位次列表；该 FG 非主官能团时为空。"""
    if spec.locant_source == "anchor_field":
        return _anchor_field_locants(oriented, spec.parent_anchor_fields[0]) or []
    if spec.locant_source == "attachment_exocyclic" and not _exocyclic_only(oriented):
        return []
    return _typed_atom_locants(oriented, spec.fg)

_FG_LOCANTS = tuple((sp.locant_kind, sp) for sp in FG_SPECS if sp.locant_kind is not None)  # (记录 kind, spec)：由 fg_registry 承载跨层一致性。

def _fg_locants(oriented: dict, n_subs: int = 0) -> list[dict]:
    """FG 位次记录 [{kind, locants, omit}]，仅产存在的项。"""
    n = oriented.get("n_carbons", 0)
    records = []
    for kind, spec in _FG_LOCANTS:
        locs = _locants_for(oriented, spec)
        if not locs:
            continue
        records.append({
            "kind": kind, "locants": locant_str_sort(locs), "omit": _omit_for(kind, oriented, n, n_subs),
        })
    return records

def _pack(oriented: dict, substituents: list) -> dict:
    """组装 parent/取代基/FG 位次/不饱和位次。"""
    n_subs = len(substituents or [])
    result =  {
        "parent": oriented, "substituents": substituents,
        "fg_locants": _fg_locants(oriented, n_subs),
        **_unsat_locants(oriented, oriented.get("n_carbons", 0)),
    }
    return result
