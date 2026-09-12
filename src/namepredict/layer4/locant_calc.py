"""L4 位次计算：将官能团/取代基附着原子映射为链上位次。"""
from __future__ import annotations
from namepredict.constants import AMINE_KINDS, OH_KINDS
from namepredict.layer1.fg_registry import FG_SPECS
from namepredict.layer4._chain_orient import _edge_min_locant, _pair_locants
from namepredict.layer4.locant_key import locant_str_sort
from namepredict.layer4.omit_locants import (
    omit_amine as _omit_amine, omit_ketone as _omit_ketone, omit_sh as _omit_sh,
)


def _typed_group_atoms(parent: dict, group: str) -> list[int]:
    """返回 principal_expression_facts 中属于指定基团类型的附着原子。"""
    facts = parent.get("principal_expression_facts")
    return sorted(facts.attachment_atoms) if facts and facts.group_class.value == group else []

def _atom_locant(chain: list[int], atom: int | None, kind: str | None, facts=None, required=False) -> int | str | None:
    """按固定标签定位次(数字 int / 字母位 "4a" 原样 str)，否则保留普通链编号——fused 桥头 4a/8a 字母位由此保留。"""
    if atom is None or atom not in chain:
        return None
    labels = (facts or {}).get("labels")
    if labels and len(labels) == len(chain):
        lbl = labels[chain.index(atom)]
        return int(lbl) if str(lbl).isdigit() else str(lbl)
    return chain.index(atom) + 1

def _fg_locant(oriented: dict, kinds: tuple, key: str) -> int | None:
    """按 key 取单个官能团附着原子并计算其链上位次。"""
    if oriented.get("kind") not in kinds:
        return None
    return _atom_locant(oriented.get("chain") or [], oriented.get(key), oriented.get("kind"), oriented.get("numbering_scaffold"), oriented.get("numbering_scaffold_required", False))

def _typed_atom_locants(oriented: dict, group: str) -> list[int]:
    """返回指定基团全部附着原子在链上的位次列表。"""
    chain = oriented.get("chain") or []
    atoms = _typed_group_atoms(oriented, group)
    kind, facts = oriented.get("kind"), oriented.get("numbering_scaffold")
    required = oriented.get("numbering_scaffold_required", False)
    return locant_str_sort(loc for atom in atoms
                           if (loc := _atom_locant(chain, atom, kind, facts, required)) is not None)


def _oh_locant(oriented: dict) -> int | None:
    """计算醇羟基的唯一 locant；多羟基或缺失时回退到单点字段。"""
    locs = _typed_atom_locants(oriented, "alcohol")
    return locs[0] if len(locs) == 1 else _fg_locant(oriented, OH_KINDS, "oh_c_idx")

def _sh_locant(oriented: dict) -> int | None:
    """计算硫醇巯基的 locant。"""
    return _fg_locant(oriented, ("thiol",), "sh_c_idx")

def _oriented_pair_locants(oriented: dict, kinds, key: str) -> list[int] | None:
    """按 kind 过滤后计算 key 对应的多位次列表。"""
    if oriented.get("kind") not in kinds:
        return None
    locs = _pair_locants(oriented.get("chain") or [], oriented.get(key))
    return list(locs) if locs else None


def _oh_locants(oriented: dict) -> list[int] | None:
    """返回全部醇羟基位次；无 typed 原子时回退到 oh_c_idxs。"""
    locs = _typed_atom_locants(oriented, "alcohol")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("alcohol",), "oh_c_idxs")

def _amine_pair_locants(oriented: dict) -> list[int] | None:
    """返回全部氨基位次；无 typed 原子时回退到 amine_c_idxs。"""
    locs = _typed_atom_locants(oriented, "amine")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("amine",), "amine_c_idxs")

def _amine_locant(oriented: dict) -> int | None:
    """计算氨基的唯一 locant；多氨基或缺失时回退到单点字段。"""
    locs = _typed_atom_locants(oriented, "amine")
    return locs[0] if len(locs) == 1 else _fg_locant(oriented, AMINE_KINDS, "amine_c_idx")

def _ketone_locant(oriented: dict) -> int | None:
    """计算酮羰基的唯一 locant（保留稠环母体按 standard 标签取，如呫吨 9、甾体 3；未登记标签回退链位）。"""
    atoms = _typed_group_atoms(oriented, "ketone")
    if len(atoms) == 1:
        return _atom_locant(oriented.get("chain") or [], atoms[0], oriented.get("kind"),
                            oriented.get("numbering_scaffold"),
                            oriented.get("numbering_scaffold_required", False))
    return _fg_locant(oriented, ("ketone",), "ketone_c_idx")

def _ketone_pair_locants(oriented: dict) -> list[int] | None:
    """返回酮羰基位次列表，多羰基时优先用附着原子（保留母体按 standard 标签）。"""
    locs = _typed_atom_locants(oriented, "ketone")
    if locs:
        return locs
    return _oriented_pair_locants(oriented, ("ketone",), "ketone_c_idxs")

def _has_parent_ene(oriented: dict) -> bool:
    """判断 parent 是否携带双键（单个或列表）。"""
    return bool(oriented.get("double_bond") or oriented.get("double_bonds"))

def _ene_locant(oriented: dict) -> int | None:
    """计算单个双键的 locant；多烯时返回 None。"""
    if not _has_parent_ene(oriented):
        return None
    db = oriented.get("double_bond")
    if db is None:
        return None  # 多烯：位次由 ene_locants 列表承载
    return _edge_min_locant(oriented.get("chain") or [], db)

def _yne_locant(oriented: dict) -> int | None:
    """计算三键的 locant；优先用 yne_locants 列表首位。"""
    locs = yne_locants(oriented)
    if locs:
        return locs[0]
    if oriented.get("kind") == "alkyne" or oriented.get("triple_bond"):
        return _edge_min_locant(oriented.get("chain") or [], oriented.get("triple_bond"))
    return None

def _has_parent_yne(oriented: dict) -> bool:
    """判断 parent 是否携带三键（单个或列表）。"""
    return bool(oriented.get("triple_bond") or oriented.get("triple_bonds"))

def _omit_oh(oh_pos, n_carbons, kind=None, parent=None, n_subs=0):
    """委托 omit_locants.omit_oh 判定羟基位次是否省略。"""
    from namepredict.layer4.omit_locants import omit_oh as _core
    return _core(oh_pos, n_carbons, kind, parent, n_subs,
                 has_ene=_has_parent_ene, has_yne=_has_parent_yne)

def _omit_unsat(n_carbons, kind=None, parent=None, triple=False):
    """委托 omit_locants.omit_unsat 判定不饱和位次是否省略；triple 选择炔规则。"""
    from namepredict.layer4.omit_locants import omit_unsat as _core
    return _core(n_carbons, kind, parent, has_ene=_has_parent_ene, has_yne=_has_parent_yne,
                 triple=triple)

def _sub_locant(chain: list[int], attach: int, kind: str | None, facts=None) -> int:
    """返回取代基附着原子的位次；无位次时取 0。"""
    loc = _atom_locant(chain, attach, kind, facts)
    return 0 if loc is None else loc

def _with_locants(chain: list[int], substituents: list, kind: str | None = None, facts=None) -> list:
    """为每个取代基附加其 locant 后返回新列表。"""
    return [{**s, "locant": _sub_locant(chain, s["attach_idx"], kind, facts)} for s in substituents]

from namepredict.layer4._chain_orient import _bond_min_locs

def ene_locants(oriented: dict) -> list[int] | None:
    """返回全部双键端点较小位次的排序列表。"""
    bonds = oriented.get("double_bonds")
    if not bonds:
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None


def yne_locants(oriented: dict) -> list[int] | None:
    """返回全部三键端点较小位次的排序列表（多炔）。"""
    bonds = oriented.get("triple_bonds")
    if not bonds:
        return None
    locs = _bond_min_locs(oriented.get("chain") or [], bonds)
    return list(locs) if locs else None

def _unsat_locants(oriented: dict, n: int) -> dict:
    """打包烯/炔位次及其省略标志。"""
    kind = oriented.get("kind")
    return {
        "ene_locant": _ene_locant(oriented),
        "ene_locants": ene_locants(oriented),
        "omit_ene_locant": _omit_unsat(n, kind, oriented),
        "yne_locant": _yne_locant(oriented),
        "yne_locants": yne_locants(oriented),
        "omit_yne_locant": _omit_unsat(n, kind, oriented, triple=True),
    }

def _sh_locants_list(oriented: dict) -> list[int] | None:
    """巯基位次列表；多硫醇用全部附着原子，单硫醇回退 sh_c_idx。"""
    locs = _typed_atom_locants(oriented, "thiol")
    if locs:
        return locs
    loc = _sh_locant(oriented)
    return [loc] if loc is not None else None

def _omit_ket_loc(oriented: dict, n_subs: int) -> bool:
    """判定酮位次是否省略：单酮按 scaffold 判断。"""
    single = len(_typed_group_atoms(oriented, "ketone")) == 1
    return _omit_ketone(oriented.get("kind"), n_subs, oriented,
                        has_ene=_has_parent_ene, single=single)

def _omit_for(kind: str, oriented: dict, n: int, n_subs: int) -> bool:
    """FG 记录 omit 标志:环状判断由 omit_locants 基于 scaffold_id 完成（不虚构 cyclo* kind）。"""
    if kind == "oh":
        return _omit_oh(_oh_locant(oriented), n, oriented.get("kind"), oriented, n_subs)
    if kind == "amine":
        return _omit_amine(_amine_locant(oriented), n, oriented.get("kind"), n_subs, oriented)
    if kind == "ketone":
        return _omit_ket_loc(oriented, n_subs)
    if kind == "sh":
        return _omit_sh(_sh_locant(oriented), n)
    return False

def _amine_fg_locants(oriented: dict) -> list[int] | None:
    """返回氨基位次列表；缺省时用单个 locant 包装。"""
    locs = _amine_pair_locants(oriented) or []
    if locs:
        return locs
    loc = _amine_locant(oriented)
    return [loc] if loc is not None else None

def _ketone_fg_locants(oriented: dict) -> list[int] | None:
    """返回酮羰基位次列表；缺省时用单个 locant 包装。"""
    locs = _ketone_pair_locants(oriented) or []
    if locs:
        return locs
    loc = _ketone_locant(oriented)
    return [loc] if loc is not None else None

def _acid_fg_locants(oriented: dict) -> list[int] | None:
    """返回 exocyclic 酸的环上附着原子位次（羧基碳在环外）；多羧酸（multiplicity≥2）取全部附着原子位次列表（与 _aldehyde_fg_locants 一致），单酸回退 ring_attach_idx 单点；无 typed 记录时返回 None。"""
    locs = _typed_atom_locants(oriented, "acid")
    if locs:
        return locs
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


def _amide_fg_locants(oriented: dict) -> list[int] | None:
    """返回 exocyclic 酰胺的环上附着原子位次（羰基碳在环外）。"""
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


def _ester_fg_locants(oriented: dict) -> list[int] | None:
    """返回 exocyclic 酯的环上附着原子位次（酯羰基碳在环外）。"""
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


def _nitrile_fg_locants(oriented: dict) -> list[int] | None:
    """返回 exocyclic 腈的环上附着原子位次（腈基碳在环外，-carbonitrile 用）。"""
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


def _aldehyde_fg_locants(oriented: dict) -> list[int] | None:
    """返回外环醛(-CHO)在环上的附着原子位次（单/多 -carbaldehyde 通用）；开链醛的醛碳在骨架内（in_skeleton），本路径不计数。"""
    facts = oriented.get("principal_expression_facts")
    if not facts or facts.group_class.value != "aldehyde" or facts.relation.value != "exocyclic":
        return None
    return _typed_atom_locants(oriented, "aldehyde")


def _radical_locants(oriented: dict) -> list[int] | None:
    """返回自由基连接点位次列表（radical_c_idx 在定向链上的位次，连接点隐含 1 或按杂环编号）。"""
    idx = oriented.get("radical_c_idx")
    if idx is None:
        return None
    loc = _atom_locant(oriented.get("chain") or [], idx, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None


def _acyl_locants(oriented: dict) -> list[int] | None:
    """返回酰基头位次列表：开链酰基头碳 acyl_c_idx（P-65.1.7.2 酸碳恒 locant 1），exocyclic 环酰基（羰基头在环外）回退环附着原子 ring_attach_idx（furan-2-carbonyl 的 2）。"""
    idx = oriented.get("acyl_c_idx")
    chain = oriented.get("chain") or []
    if idx is not None:
        loc = _atom_locant(chain, idx, oriented.get("kind"),
                           oriented.get("numbering_scaffold"),
                           oriented.get("numbering_scaffold_required", False))
        if loc is not None:
            return [loc]
    attach = oriented.get("ring_attach_idx")
    if attach is None:
        return None
    loc = _atom_locant(chain, attach, oriented.get("kind"),
                       oriented.get("numbering_scaffold"),
                       oriented.get("numbering_scaffold_required", False))
    return [loc] if loc is not None else None

_LOCANT_FNS = {  # 位次记录 kind → 位次函数（kind 为 fg_registry.locant_kind，跨 L4/L5 一致性由 spec 承载）。新增带位次的 FG：在 _LOCANT_FNS 补函数 + fg_registry 设 locant_kind，未登记则 KeyError 显式暴露。
    "oh": _oh_locants,
    "amine": _amine_fg_locants,
    "ketone": _ketone_fg_locants,
    "sh": _sh_locants_list,
    "aldehyde": _aldehyde_fg_locants,
    "acid": _acid_fg_locants,
    "amide": _amide_fg_locants,
    "ester": _ester_fg_locants,
    "nitrile": _nitrile_fg_locants,
    "radical": _radical_locants,
    "acyl": _acyl_locants,
}
_FG_LOCANTS = tuple(
    (sp.locant_kind, _LOCANT_FNS[sp.locant_kind])
    for sp in FG_SPECS if sp.locant_kind is not None
)

def _fg_locants(oriented: dict, n_subs: int = 0) -> list[dict]:
    """FG 位次记录: [{kind, locants, omit}] — 稀疏,只产实际存在的 principal FG."""
    n = oriented.get("n_carbons", 0)
    records = []
    for kind, locs_fn in _FG_LOCANTS:
        locs = locs_fn(oriented)
        if not locs:
            continue
        records.append({
            "kind": kind, "locants": locant_str_sort(locs), "omit": _omit_for(kind, oriented, n, n_subs),
        })
    return records

def _pack(oriented: dict, substituents: list) -> dict:
    """组装 parent/取代基/FG 位次/不饱和位次与相对立体化学。"""
    from namepredict.layer4.cyclo_relative_stereo import relative_stereo_facts
    facts = relative_stereo_facts(oriented)
    merged = {**oriented, **facts}
    n_subs = len(substituents or [])
    return {
        "parent": merged, "substituents": substituents,
        "fg_locants": _fg_locants(merged, n_subs),
        **_unsat_locants(merged, merged.get("n_carbons", 0)), **facts,
    }
