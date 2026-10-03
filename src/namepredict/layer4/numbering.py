"""L4 编号入口：定向编号并组装最终 result 包。"""
from __future__ import annotations
from namepredict.constants import C, HYDRO_MULT_N, MULT_EN, MULT_ZH
from namepredict.layer4.indicated_hydrogen import indicated_hydrogen, saturated_ring_atoms
from namepredict.layer4.numbering_engine import orient_numbering
from namepredict.layer4.locant_calc import _pack, _with_locants
from namepredict.layer4.locant_calc import atom_locant, locant_key


def _fallback_hydro_atoms(parent: dict) -> frozenset:
    """稠环加氢位回退：取环内单键连邻环的饱和 sp3 位（P-31.2.2）。"""
    mol = parent.get("mol")
    chain = list(parent.get("chain") or ())
    if mol is None or not chain or parent.get("fused_tree") is None:
        return frozenset()
    out = frozenset(saturated_ring_atoms(mol, set(chain)))  # 芳香稠环上的 NH 位也是加氢位（P-31.2.2）
    carbons = frozenset(a for a in out if mol.GetAtomWithIdx(a).GetAtomicNum() == C)
    return carbons if len(carbons) != len(out) and len(carbons) in HYDRO_MULT_N else out


def _locant_of(chain, labels):
    """链原子 → 排序键闭包（按整体编号标签定位次）。"""
    facts = {"labels": labels}
    return lambda a: locant_key(str(atom_locant(chain, a, facts)))


def _lowest_extra_to_indicated(packed: dict, labels, hydro: frozenset) -> frozenset:
    """把加氢位中最低位次改用指示氢表达 """
    if not hydro:
        return hydro
    mol = packed.get("mol")
    chain = list(packed.get("chain") or ())
    if mol is None or not chain:
        return hydro
    chain_set = set(chain)
    _key = _locant_of(chain, labels)
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
    """加氢位为奇数时最低者改用指示氢（P-51.1.1.4、P-58.2.1.2）。"""
    n = len(hydro or ())
    if not n or n in HYDRO_MULT_N:
        return hydro
    chain = list(packed.get("chain") or ())
    if not chain or any(a not in chain for a in hydro):  # 位次表达不全，本层补救不了
        return hydro
    lowest = min(hydro, key=_locant_of(chain, labels))
    return frozenset(hydro - {lowest})


def hydro_prefix(chain, labels, hydro_atoms) -> tuple[str, str]:
    """返回加氢前缀 (en, zh)：全加氢时省略位次（P-14.3.4.5）。"""
    chain = list(chain or ())
    hydro = set(hydro_atoms or ())
    if not chain or not hydro:
        return "", ""
    n = len(hydro)
    if n not in HYDRO_MULT_N:
        return "", ""
    if set(chain) <= hydro:  # 完全氢化：省略全部位次（P-14.3.4.5）
        return f"{MULT_EN[n]}hydro", f"{MULT_ZH[n]}氢"
    facts = {"labels": labels}
    locs = [str(atom_locant(chain, a, facts)) for a in hydro if a in chain]
    if len(locs) != n:  # 加氢原子不全在编号链内：位次无法完整表达，放弃而非给错名
        return "", ""
    locs.sort(key=locant_key)
    loc = ",".join(locs)
    return f"{loc}-{MULT_EN[n]}hydro", f"{loc}-{MULT_ZH[n]}氢"


def number(parent: dict, substituents: list) -> dict:
    """对 parent 定向编号并组装位次结果（含指示氢前缀 P-58.2.1）。"""
    chain = orient_numbering(parent, substituents)
    if chain is None:  # 定向失败（如桥环候选不可判定）：显式失败，不让 None 顺流触发 TypeError
        raise ValueError("numbering_failed")
    oriented = {**parent, "chain": chain}
    result = _pack(oriented, _with_locants(chain, substituents, oriented.get("numbering_scaffold")))
    packed = result.get("parent") or {}
    labels = (packed.get("numbering_scaffold") or {}).get("labels")
    hydro = packed.get("hydro_atoms") or frozenset()
    fb = _fallback_hydro_atoms(packed)  # 未注册稠环推导
    if not hydro or (set(hydro) < fb and (len(fb) in HYDRO_MULT_N or len(fb) - 1 in HYDRO_MULT_N)):
        hydro = fb  # 奇数个饱和位也可回退：随后 _odd_hydro_to_indicated 会把最低位次改由指示氢表达
    # print(hydro)
    hydro_all = hydro
    hydro = _lowest_extra_to_indicated(packed, labels, hydro)
    # print(hydro)
    hydro = _odd_hydro_to_indicated(packed, labels, hydro)
    # print(hydro)
    pre = hydro_prefix(packed.get("chain"), labels, hydro)
    # 加氢位被搬进指示氢且前缀不再提它：仅当该位是模板未隐含的加氢位（保留名已含的 CH2 不算）时，
    # 才要求 L5 写出指示氢，否则饱和度整段丢失
    moved = hydro_all - hydro
    packed["hydro_fallback"] = bool(moved) and bool(
        moved & _extra_hydrogenated(packed))
    if not pre[0]:  # hydro 位次表达不出（奇数值/超表/不在链内）则整体退回指示氢，不产半截名
        hydro, pre = frozenset(), ("", "")
    from namepredict.layer2.hantzsch_widman import is_hw_scaffold
    extra = _extra_indicated(packed)
    packed["lambda_locants"] = _lambda_locants(packed.get("mol"), packed.get("chain"), labels)
    packed["indicated_h_locants"] = indicated_hydrogen(
        packed.get("mol"), packed.get("chain"), labels, hydro, extra,
        packed.get("scaffold_id"))
    packed["indicated_h_forced"] = bool(extra) or is_hw_scaffold(packed.get("scaffold_id"))
    packed["ind_h_carbon_ok"] = _ind_h_carbon_scaffold_ok(packed.get("scaffold_id"))
    packed["hydro_prefix"] = pre
    # print("layer4result",result,"\n\n\n\n\n")
    return result


def _lambda_locants(mol, chain, labels) -> tuple[tuple[str, int, int], ...]:
    """整体编号下需 λ/δ 标注的骨架原子 → ((位次, 键数|0, δ 双键数|0), …)。

    环内 λ 只用于 HW_RING_LAMBDA_Z 元素：硫族/氮氧的环内高价由 dioxo/oxide 前缀表达，
    名称里不复标 λ（与 P-22.2.7 同一口径）。δ 计该原子直接相连的连续双键数（P-25.7.2）。
    """
    from namepredict.layer2.hantzsch_widman import HW_RING_LAMBDA_Z
    from namepredict.tools.lambda_notation import bonding_number, is_lambda_marked

    if mol is None or not chain or not labels:
        return ()
    ring = set(chain)
    out: list[tuple[str, int, int]] = []
    for idx, loc in zip(chain, labels):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            continue
        n = bonding_number(atom) if atom.GetAtomicNum() in HW_RING_LAMBDA_Z and is_lambda_marked(atom) else 0
        delta = sum(1 for b in atom.GetBonds()  # 只数骨架内的连续双键，=O/=S 等环外双键不计
                    if b.GetBondTypeAsDouble() == 2.0 and not b.GetIsAromatic()
                    and b.GetOtherAtomIdx(idx) in ring)
        delta = delta if delta >= 2 else 0  # 单个双键常规表达，δ 只标连续形式（P-25.7.2）
        if n or delta:
            out.append((str(loc), n, delta))
    return tuple(out)


def _ind_h_carbon_scaffold_ok(scaffold_id: str | None) -> bool:
    """母体是否为可写碳位指示氢的保留杂环名（单环杂芳/稠合杂芳；苯并二氧戊环、金刚烷等除外）。"""
    if not scaffold_id:
        return False
    from namepredict.layer2.ring_scaffold import get_spec

    sp = get_spec(scaffold_id)
    # benzodioxole 的 CH2、adamantane/mono_carbo 的饱和位由保留名本身表达，不得补指示氢（P-58.2.1）
    return sp is not None and getattr(sp, "naming_class", "") not in (
        "benzodioxole", "adamantane", "mono_carbo")


def _extra_indicated(packed: dict) -> frozenset:
    """保留母体名未隐含、须显式标指示氢的环位（P-58.2.1）。"""
    from namepredict.layer2.ring_scaffold import extra_indicated_atoms

    sid, match, mol = packed.get("scaffold_id"), packed.get("scaffold_match"), packed.get("mol")
    return extra_indicated_atoms(mol, sid, match) if sid and match and mol is not None else frozenset()


def _extra_hydrogenated(packed: dict) -> frozenset:
    """氢数多于保留母体模板的环位（P-58.2.1 加氢指示氢）。"""
    from namepredict.layer2.ring_scaffold import extra_hydrogenated_atoms

    sid, match, mol = packed.get("scaffold_id"), packed.get("scaffold_match"), packed.get("mol")
    return extra_hydrogenated_atoms(mol, sid, match) if sid and match and mol is not None else frozenset()
