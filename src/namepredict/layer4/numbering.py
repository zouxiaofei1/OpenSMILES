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
    out = frozenset(i for i in saturated_ring_atoms(mol, set(chain))
                    if not mol.GetAtomWithIdx(i).GetIsAromatic())
    carbons = frozenset(a for a in out if mol.GetAtomWithIdx(a).GetAtomicNum() == C)
    return carbons if len(carbons) != len(out) and len(carbons) in HYDRO_MULT_N else out


def _lowest_extra_to_indicated(packed: dict, labels, hydro: frozenset) -> frozenset:
    """把加氢位中最低位次改用指示氢表达 """
    if not hydro:
        return hydro
    mol = packed.get("mol")
    chain = list(packed.get("chain") or ())
    if mol is None or not chain:
        return hydro
    chain_set = set(chain)
    facts = {"labels": labels}

    def _key(a: int):
        return locant_key(str(atom_locant(chain, a, facts)))

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
    from namepredict.constants import HYDRO_MULT_N

    n = len(hydro or ())
    if not n or n in HYDRO_MULT_N:
        return hydro
    chain = list(packed.get("chain") or ())
    if not chain or any(a not in chain for a in hydro):  # 位次表达不全，本层补救不了
        return hydro
    facts = {"labels": labels}
    lowest = min(hydro, key=lambda a: locant_key(str(atom_locant(chain, a, facts))))
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
    oriented = {**parent, "chain": chain}
    result = _pack(oriented, _with_locants(chain, substituents, oriented.get("numbering_scaffold")))
    packed = result.get("parent") or {}
    labels = (packed.get("numbering_scaffold") or {}).get("labels")
    hydro = packed.get("hydro_atoms") or frozenset()
    fb = _fallback_hydro_atoms(packed)  # 未注册稠环推导
    if not hydro or (len(fb) in HYDRO_MULT_N and set(hydro) < fb):
        hydro = fb
    # print(hydro)
    hydro = _lowest_extra_to_indicated(packed, labels, hydro)
    # print(hydro)
    hydro = _odd_hydro_to_indicated(packed, labels, hydro)
    # print(hydro)
    pre = hydro_prefix(packed.get("chain"), labels, hydro)
    if not pre[0]:  # hydro 位次表达不出（奇数值/超表/不在链内）则整体退回指示氢，不产半截名
        hydro, pre = frozenset(), ("", "")
    extra = _extra_indicated(packed)
    packed["indicated_h_locants"] = indicated_hydrogen(
        packed.get("mol"), packed.get("chain"), labels, hydro, extra)
    packed["indicated_h_forced"] = bool(extra)
    packed["hydro_prefix"] = pre
    # print("layer4result",result,"\n\n\n\n\n")
    return result


def _extra_indicated(packed: dict) -> frozenset:
    """保留母体名未隐含、须显式标指示氢的环位（P-58.2.1）。"""
    from namepredict.layer2.ring_scaffold import extra_indicated_atoms

    sid, match, mol = packed.get("scaffold_id"), packed.get("scaffold_match"), packed.get("mol")
    return extra_indicated_atoms(mol, sid, match) if sid and match and mol is not None else frozenset()
