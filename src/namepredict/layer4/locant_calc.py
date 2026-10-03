"""L4 位次计算：附着原子 → 链上位次，含 FG/烯炔位次省略与候选比较键。"""
from __future__ import annotations
from namepredict.layer1.fg_registry import FG_SPECS
import re

def locant_key(x) -> tuple[int, str]:
    """locant → 排序键：数字按数值，字母尾与撇号作次级键。"""
    m = re.match(r"(\d+)([a-z]*)(\'*)", str(x))
    return (int(m.group(1)), m.group(2) + m.group(3)) if m else (0, "")


def locant_str_sort(locs) -> list:
    """按 locant_key 排序 locant 集合(兼容 int 与字母位)。"""
    return sorted(locs, key=locant_key)


def _as_bond_pairs(value) -> tuple:
    """把键取值统一成键对序列：单个键对或键对序列均可。"""
    if not value or isinstance(value, (int, str)):
        return ()
    return tuple(value) if isinstance(value[0], (tuple, list)) else (tuple(value),)


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

atom_locant = _atom_locant  # 公开别名：各层统一用此名，_atom_locant 保留供既有 import

def _atom_locants(oriented: dict, atoms) -> list[int]:
    """把一组骨架原子映射为位次列表（保留 fused 字母位，跳过不在链表中的原子）。"""
    chain, facts = oriented.get("chain") or [], oriented.get("numbering_scaffold")
    return [loc for atom in atoms
            if (loc := _atom_locant(chain, atom, facts)) is not None]


def _typed_atom_locants(oriented: dict, group: str) -> list[int]:
    """返回指定基团全部附着原子在链上的位次列表。"""
    return locant_str_sort(_atom_locants(oriented, _typed_group_atoms(oriented, group)))

def _occ_attachment(oriented: dict, anchors: set, atoms: set) -> int | None:
    """occurrence 锚点 → 骨架内附着原子（骨架外取骨架内邻居）。"""
    inside = anchors & atoms
    if inside:
        return min(inside)
    mol = oriented.get("mol")
    if mol is None:
        return min(anchors) if anchors else None
    neighbours = [n.GetIdx() for i in anchors
                  for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in atoms]
    return min(neighbours) if neighbours else None


def _occurrence_locants(oriented: dict, spec) -> list[int]:
    """逐 occurrence 求附着位次（同一原子承载多个同类 FG 时保留重数）。"""
    facts = oriented.get("principal_expression_facts")
    if facts is None or facts.group_class.value != spec.fg:
        return []
    ids = set(oriented.get("covered_principal_ids") or ())
    chain, atoms = oriented.get("chain") or [], set(oriented.get("chain") or ())
    scaffold = oriented.get("numbering_scaffold")
    out = []
    for occ in oriented.get("principal_occurrences") or ():
        if occ.id not in ids:
            continue
        anchors = {i for key in spec.anchors for i in (occ.payload.get(key) or ())}
        atom = _occ_attachment(oriented, anchors, atoms)
        loc = _atom_locant(chain, atom, scaffold)
        if loc is not None:
            out.append(loc)
    return out


def _expand_shared_locants(oriented: dict, spec, locs: list) -> list:
    """同位次重复 FG 补回重复位次（偕二醇 propane-2,2-diol）。"""
    facts = oriented.get("principal_expression_facts")
    mult = facts.multiplicity if facts is not None and facts.group_class.value == spec.fg else None
    if not mult or mult <= len(locs):
        return locs
    per_occ = locant_str_sort(_occurrence_locants(oriented, spec))
    return per_occ if len(per_occ) == mult else locs


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
    pairs = _as_bond_pairs(oriented.get(f"{key}_bond")) or _as_bond_pairs(oriented.get(f"{key}_bonds"))
    return list(pairs) or None


def _bond_locant(chain: list[int], bond) -> int | str:
    """单键位次：相邻取较小者，否则写复合位次 x(y)。"""
    lo, hi = sorted((chain.index(bond[0]) + 1, chain.index(bond[1]) + 1))
    return lo if hi - lo == 1 else f"{lo}({hi})"


def _bond_min_locs(chain: list[int], bonds) -> tuple | None:
    """返回全部键位次的排序元组（跨位次键为复合位次）；有键无位次则 None。"""
    if not bonds:
        return None
    locs = []
    for b in bonds:
        if not b or b[0] not in chain or b[1] not in chain:
            return None
        locs.append(_bond_locant(chain, b))
    return tuple(sorted(locs, key=locant_key))


def _bond_locants(oriented: dict, b: str) -> list | None:
    """返回某类不饱和键位次的排序列表；缺失返回 None。"""
    bonds = _unsat_bonds(oriented, b)
    if not bonds:
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None


def ene_locants(oriented: dict) -> list | None:
    """返回全部双键位次的排序列表（单烯亦为单元素列表；跨位次键写作 x(y)）。"""
    return _bond_locants(oriented, "ene")


def yne_locants(oriented: dict) -> list | None:
    """返回全部三键位次的排序列表（单炔亦为单元素列表；跨位次键写作 x(y)）。"""
    return _bond_locants(oriented, "yne")

def _unsat_locants(oriented: dict, n: int) -> dict:
    """打包烯/炔位次及其省略标志。"""
    kind = oriented.get("kind")
    return {
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": omit_unsat(n, kind, oriented),
        "yne_locants": yne_locants(oriented),
        "omit_yne_locant": omit_unsat(n, kind, oriented, triple=True),
    }

_FG_GROUP = {"alcohol": "alcohol", "amine": "amine", "ketone": "ketone", "thiol": "thiol"}  # 记录 kind 到 FG 类别映射（FG_SPECS）


def _omit_for(kind: str, oriented: dict, n: int, n_subs: int) -> bool:
    """FG 记录 omit 标志；环状判断交 omit_fg_locant。"""
    group = _FG_GROUP.get(kind)
    if group is None:
        return False
    single = kind != "ketone" or len(_typed_group_atoms(oriented, "ketone")) == 1
    return omit_fg_locant(True, n, oriented, n_subs, single=single)

_FG_LOCANTS = tuple((sp.fg, sp) for sp in FG_SPECS if sp.fg is not None)  # (记录 kind, spec)：由 fg_registry 承载跨层一致性。

def _fg_locants(oriented: dict, n_subs: int = 0) -> list[dict]:
    """FG 位次记录 [{kind, locants, omit}]，仅产存在的项。"""
    n = oriented.get("n_carbons", 0)
    records = []
    for kind, spec in _FG_LOCANTS:
        locs = _typed_atom_locants(oriented, spec.fg)
        if not locs:
            continue
        if spec.locant_source == "attachment":
            locs = _expand_shared_locants(oriented, spec, locs)
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


def omit_fg_locant(
    pos: int | None, n_carbons: int, parent: dict | None = None, n_subs: int = 0, *,
    single: bool = True,
) -> bool:
    """判定主官能团位次是否省略（P-14.3.4 / 环单 FG）：环状无取代省。"""
    if ((parent or {}).get("scaffold_id") == "carbocycle" and not (parent or {}).get("fused_tree")
            and pos is not None and single):
        return n_subs == 0
    return pos == 1 and n_carbons <= 2


def omit_unsat(
    n_carbons: int, kind: str | None = None, parent: dict | None = None, *,
    triple: bool = False,
) -> bool:
    """判定烯/炔位次是否省略（环单烯或短链）；triple 选择炔规则。"""
    if kind == "heterane":  # 杂原子链：二核与三核的单一不饱和键省略位次（P-14.3.4.2(d)）
        return n_carbons <= 3
    if kind == "alkane" and (parent or {}).get("scaffold_id") == "carbocycle":  # 纯烃环单烯位次隐含省略；环多烯保留位次。
        if not (parent or {}).get("double_bonds"):
            return True
    return n_carbons <= (3 if triple else 2)  # 烯 ≤C2、炔 ≤C3 时位次 '1' 省略（P-14.3.4.2(d)）。


def suffix_locant_set(numbered: dict) -> tuple:
    """P-44.1.1：principal 特征基团位次集合（P-14.3.5）。"""
    from namepredict.layer4.numbering_engine import _principal_atoms  # 函数内导入：避开与 numbering_engine 的环

    parent = numbered.get("parent") or {}
    locs = _atom_locants(parent, _principal_atoms(parent))
    return tuple(sorted(locant_key(x) for x in locs))


def prefix_locant_set(numbered: dict) -> tuple:
    """P-45.2.2：前缀取代基位次集合（P-14.3.5）；无位次前缀不入键。"""
    locs = [s["locant"] for s in (numbered.get("substituents") or []) if s.get("locant") is not None]
    return tuple(sorted(locant_key(x) for x in locs))
