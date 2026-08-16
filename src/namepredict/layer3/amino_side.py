"""非胺母体上的仲氨基取代基（例如乙醇上的 N-苄基）。"""
from __future__ import annotations


def _make_amino(attach: int, n_idx: int, en: str = "amino", zh: str = "氨基",
                atoms: list[int] | None = None, paren: bool = False) -> dict:
    """构造氨基取代基字典（kind=amino）。"""
    return {
        "kind": "amino", "attach_idx": attach,
        "atoms": [n_idx] if atoms is None else list(atoms),
        "en": en, "zh": zh, "paren": paren,
    }


def _one_amino(info: dict, a: dict, chain_set: set[int], owned) -> dict | None:
    """判断单个胺是否位于链上并返回其氨基取代基字典。 """
    if a.get("c_idx") not in chain_set:
        return None
    if any(c not in chain_set for c in (a.get("c_idxs") or [])):
        return None
    return _make_amino(a["c_idx"], a["n_idx"])


def _principal_attachments(parent: dict, group: str) -> frozenset[int]:
    """返回某基团类型的主表达式连接原子集合。"""
    facts = parent.get("principal_expression_facts")
    return facts.attachment_atoms if facts and facts.group_class.value == group else frozenset()


def _extract_aminos(info: dict, parent: dict, parent_nh2_kinds: set) -> list[dict]:
    """提取未被主基团占用的氨基取代基列表。"""
    principal = _principal_attachments(parent, "amine")
    if parent.get("kind") in parent_nh2_kinds and not principal:
        return []
    chain_set = set(parent.get("chain") or [])
    owned = parent.get("owned_atoms")
    out: list[dict] = []
    for a in info.get("amines") or []:
        one = _one_amino(info, a, chain_set, owned)
        if one is not None and one["attach_idx"] not in principal:
            out.append(one)
    return out
