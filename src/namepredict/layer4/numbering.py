"""L4 编号入口：定向编号并组装最终 result 包。"""
from __future__ import annotations
from namepredict.constants import C, HYDRO_MULT_N, MULT_EN, MULT_ZH
from namepredict.layer4.indicated_hydrogen import indicated_hydrogen_prefix, saturated_ring_atoms
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.layer4.locant_calc import _pack, _with_locants
from namepredict.layer4.locant_calc import locant_key


def _fallback_hydro_atoms(parent: dict) -> frozenset:
    """无保留模板可比的稠环（未注册母体）的加氢位回退：取环系内仅以单键连邻环原子、带氢且非芳香位的 sp3 位（P-31.2.2）。 """
    mol = parent.get("mol")
    chain = list(parent.get("chain") or ())
    if mol is None or not chain or parent.get("fused_tree") is None:
        return frozenset()
    out = frozenset(i for i in saturated_ring_atoms(mol, set(chain))
                    if not mol.GetAtomWithIdx(i).GetIsAromatic())
    carbons = frozenset(a for a in out if mol.GetAtomWithIdx(a).GetAtomicNum() == C)
    return carbons if len(carbons) != len(out) and len(carbons) in HYDRO_MULT_N else out


def _lowest_extra_to_indicated(packed: dict, labels, hydro: frozenset) -> frozenset:
    """把加氢位中最低位次改用指示氢表达（P-31.2.2、P-58.2.1.2、P-58.2.2.2：指示氢优先得低位次，其余仍是 hydro 位）。

    与现行指示氢位次最高者互换，只改名次归属、不改变两边计数，故 hydro 倍增前缀的合法性与奇数回退都不受影响（如 1,4-dihydro-2H- → 2,4-dihydro-1H-）。无 hydro 位或无指示氢位时不重排。
    """
    if not hydro:
        return hydro
    mol = packed.get("mol")
    chain = list(packed.get("chain") or ())
    if mol is None or not chain:
        return hydro
    chain_set = set(chain)
    use_labels = bool(labels) and len(labels) == len(chain)

    def _key(a: int):
        return locant_key(labels[chain.index(a)] if use_labels else str(chain.index(a) + 1))

    hydro = set(hydro)
    sats = {a for a in set(saturated_ring_atoms(mol, chain_set)) | hydro if a in chain_set}
    indicated = sats - hydro  # 现行指示氢位（上游按元素切分后的剩余）
    if not indicated:
        return frozenset(hydro)
    lowest = min(sats, key=_key)
    if lowest in hydro:  # 最低位次落在 hydro 里 -> 与指示氢位中位次最高者互换
        hydro = (hydro - {lowest}) | {max(indicated, key=_key)}
    return frozenset(hydro)


def _odd_hydro_to_indicated(packed: dict, labels, hydro: frozenset) -> frozenset:
    """加氢位数不成倍增（奇数）时把其中最低位次改用指示氢表达（P-51.1.1.4、P-58.2.1.2）。

    饱和环上出现无氢饱和位（偕二甲基季碳等）使其不计入 hydro，环内带 H 的饱和位遂成奇数：hydro 只取偶数个，余下最低位次那个由指示氢承载，故 4,4,7-三甲基-2,3-二氢-1H-萘 取 hydro=2,3 / 指示氢=1，而非整体放弃加氢描述退成「4,4,7-三甲基萘」。已合法（在倍增表内）时不改动。
    """
    from namepredict.constants import HYDRO_MULT_N

    n = len(hydro or ())
    if not n or n in HYDRO_MULT_N:
        return hydro
    chain = list(packed.get("chain") or ())
    if not chain or any(a not in chain for a in hydro):  # 位次表达不全，本层补救不了
        return hydro
    use_labels = bool(labels) and len(labels) == len(chain)
    lowest = min(hydro, key=lambda a: locant_key(
        labels[chain.index(a)] if use_labels else str(chain.index(a) + 1)))
    return frozenset(hydro - {lowest})


def hydro_prefix(chain, labels, hydro_atoms) -> tuple[str, str]:
    """返回加氢前缀 (en, zh)：位次取加氢原子 locant（按位次升序），倍增数为加氢原子数；环原子全加氢时省略位次（P-14.3.4.5，如 decahydronaphthalene）。无加氢或数量超表返回空串对。"""
    chain = list(chain or ())
    hydro = set(hydro_atoms or ())
    if not chain or not hydro:
        return "", ""
    n = len(hydro)
    if n not in HYDRO_MULT_N:
        return "", ""
    if set(chain) <= hydro:  # 完全氢化：省略全部位次（P-14.3.4.5）
        return f"{MULT_EN[n]}hydro", f"{MULT_ZH[n]}氢"
    use_labels = bool(labels) and len(labels) == len(chain)
    locs = [labels[chain.index(a)] if use_labels else str(chain.index(a) + 1)
            for a in hydro if a in chain]
    if len(locs) != n:  # 加氢原子不全在编号链内：位次无法完整表达，放弃而非给错名
        return "", ""
    locs.sort(key=locant_key)
    loc = ",".join(locs)
    return f"{loc}-{MULT_EN[n]}hydro", f"{loc}-{MULT_ZH[n]}氢"


def number(parent: dict, substituents: list) -> dict:
    """对 parent 定向编号，校验编号骨架事实后组装位次结果（含指示氢前缀 P-58.2.1）。"""
    chain = orient_numbering(parent, substituents)
    if parent.get("numbering_scaffold_required") and not parent.get("numbering_scaffold"):
        raise ValueError("numbering_scaffold facts required for selected scaffold")
    oriented = {**parent, "chain": chain}
    result = _pack(oriented, _with_locants(chain, substituents, oriented.get("numbering_scaffold")))
    packed = result.get("parent") or {}
    labels = (packed.get("numbering_scaffold") or {}).get("labels")
    hydro = packed.get("hydro_atoms") or frozenset()
    if not hydro:  # 未注册稠环无保留模板可比 -> hydro_atoms 缺失，回退由分子自身饱和环位推导
        hydro = _fallback_hydro_atoms(packed)
    else:  # 保留模板只写了一个 Kekulé 式，稠合/角位可能漏计（cyclopenta[a]phenanthrene 型甾体的 tetradecahydro 被写成 dodecahydro-1H,2H）；分子自身推得的饱和环位为合法倍增数且是模板集的真超集时改用它，否则保留模板集。
        fb = _fallback_hydro_atoms(packed)
        if len(fb) in HYDRO_MULT_N and set(hydro) < fb:
            hydro = fb
    hydro = _lowest_extra_to_indicated(packed, labels, hydro)
    hydro = _odd_hydro_to_indicated(packed, labels, hydro)
    pre = hydro_prefix(packed.get("chain"), labels, hydro)
    if not pre[0]:  # hydro 位次表达不出（奇数值/超表/不在链内）则整体退回指示氢，不产半截名
        hydro, pre = frozenset(), ("", "")
    extra = _extra_indicated(packed)
    packed["indicated_h"] = indicated_hydrogen_prefix(
        packed.get("mol"), packed.get("chain"), labels, hydro, extra)
    packed["indicated_h_forced"] = bool(extra)
    packed["hydro_prefix"] = pre
    return result


def _extra_indicated(packed: dict) -> frozenset:
    """保留母体名未隐含、须显式标指示氢的环位（P-58.2.1）：全芳香环系唯一的指示氢来源。"""
    from namepredict.layer2.ring_scaffold import extra_indicated_atoms

    sid, match, mol = packed.get("scaffold_id"), packed.get("scaffold_match"), packed.get("mol")
    return extra_indicated_atoms(mol, sid, match) if sid and match and mol is not None else frozenset()
