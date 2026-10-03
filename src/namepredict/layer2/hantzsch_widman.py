"""P-22.2.2 Hantzsch-Widman 杂单环命名：3-10 元环的词干、位次与名称组装。

未命中保留模板的杂单环由本模块生成词干（饱和/mancude 两形态），
位次按 P-22.2.2.1.3 的元素优先序排列，省略判定见 P-22.2.2.1.7。
"""
from __future__ import annotations

import re
from functools import lru_cache
from itertools import combinations, permutations

from rdkit import Chem

from namepredict.constants import (
    As, B, Bi, C, Ge, HW_CLASS, HW_COMPONENT_PREFIX, HW_ID, HW_MAX_VALENCE,
    HW_PREFIX_EN, HW_PREFIX_ZH, HW_SAT_SIX, HW_SAT_TAIL, HW_SIX_A, HW_SIX_B,
    HW_UNSAT_SIX, HW_UNSAT_TAIL, HW_VOWELS, HW_ZH_RING, HW_ZH_SHORT, I,
    MULT_EN, MULT_ZH, N, P, P145_SENIOR, Pb, Sb, Si, Sn, zh_numeral,
)
from namepredict.tools.lambda_notation import bonding_number, is_lambda_marked, lambda_mark

_PT = Chem.GetPeriodicTable()


def ring_heteros(zs) -> tuple[int, ...]:
    """环序元素表中的杂原子序列（保序）。"""
    return tuple(z for z in zs if z != C)


def six_group(zs) -> str:
    """六元环词干分组：由优先性最低的杂原子（名称紧邻词干者）所属组决定（P-22.2.2.1.6）。"""
    least = max(ring_heteros(zs), key=P145_SENIOR.index)
    if least in HW_SIX_A:
        return "A"
    return "B" if least in HW_SIX_B else "C"


def unsat_stem(n: int, zs) -> str | None:
    """不饱和（mancude）词干（Table 2.5、P-22.2.2.1.5.1）。"""
    if n == 3:
        return "irine" if set(ring_heteros(zs)) == {N} else "irene"
    if n == 4:
        return "ete"
    if n == 5:
        return "ole"
    if n == 6:
        return HW_UNSAT_SIX[six_group(zs)]
    return HW_UNSAT_TAIL.get(n)


def sat_stem(n: int, zs) -> str | None:
    """饱和词干（Table 2.5、P-22.2.2.1.5.2）。"""
    has_n = N in zs
    if n == 3:
        return "iridine" if has_n else "irane"
    if n == 4:
        return "etidine" if has_n else "etane"
    if n == 5:
        return "olidine" if has_n else "olane"
    if n == 6:
        return HW_SAT_SIX[six_group(zs)]
    return HW_SAT_TAIL.get(n)


def zh_stem(n: int, zs, n_db: int) -> str | None:
    """中文基干：含氮五/六元环取唑/嗪系，其余按环大小与实际环内双键数。"""
    has_n = N in zs
    if n == 5 and has_n:
        return "唑" if n_db else "唑烷"
    if n == 6 and has_n:
        return "嗪" if n_db else "嗪烷"
    ring = HW_ZH_RING.get(n)
    if ring is None:
        return None
    return ring + _zh_unsat_tail(n_db)


def _zh_unsat_tail(n_db: int) -> str:
    """中文双键词尾：0 烷、1 烯、n 二烯…"""
    if n_db <= 0:
        return "烷"
    if n_db == 1:
        return "烯"
    return f"{zh_numeral(n_db)}烯"


def elide_a(terms: list[str]) -> str:
    """按 P-22.2.2.1.1 拼接：前项末尾 'a' 在后项起首元音前省略。"""
    out = terms[0]
    for term in terms[1:]:
        if out.endswith("a") and term[:1] in HW_VOWELS:
            out = out[:-1]
        out += term
    return out


def locant_map(zs) -> dict[int, list[int]]:
    """按 P-22.2.2.1.3 引用顺序分组的 {元素: 位次升序}（zs 为已编号的环序）。"""
    out: dict[int, list[int]] = {}
    for i, z in enumerate(zs):
        if z != C:
            out.setdefault(z, []).append(i + 1)
    return {z: out[z] for z in P145_SENIOR if z in out}


def locant_string(by_z: dict[int, list[int]]) -> str:
    """位次串按引用顺序排列（非数值升序），如 1,6,2-dioxazepane 的 '1,6,2'。"""
    return ",".join(str(loc) for locs in by_z.values() for loc in locs)


def locant_string_lambda(by_z: dict[int, list[int]], lambda_at: dict[int, int] | None, *,
                         en: bool) -> str:
    """带 λ 的位次串（P-22.2.7.1：λn 紧跟杂原子位次之后，如 1λ6,3λ5）。"""
    marks = lambda_at or {}
    return ",".join(f"{loc}{lambda_mark(marks[loc], en=en)}" if loc in marks else str(loc)
                    for locs in by_z.values() for loc in locs)


HW_RING_LAMBDA_Z = frozenset({P, As, Sb, Bi, B, I, Si, Ge, Sn, Pb})  # 环内 λ 适用元素：硫族/氮氧的高价由 dioxo/oxide 前缀表达，不复标 λ


def ring_lambda_atoms(ordered, mol) -> dict[int, int]:
    """已编号环序 → {1 起位次: 键数}，仅收键数偏离标准值的环内杂原子（P-14.1.3）。"""
    out: dict[int, int] = {}
    for i, idx in enumerate(ordered):
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() in HW_RING_LAMBDA_Z and is_lambda_marked(atom):
            out[i + 1] = bonding_number(atom)
    return out


def _hetero_terms(n: int, zs, by_z: dict[int, list[int]], *, en: bool) -> list[str]:
    """元素前缀项（含倍数词）；中文在含氮五/六元环里把 N 折入唑/嗪基干，故中文跳过 N。"""
    has_n = N in zs
    folded = (not en) and has_n and n in (5, 6)
    short = has_n and n in (5, 6)  # 噁/噻短式只在唑/嗪系（含氮环）启用
    terms = []
    for z, locs in by_z.items():
        if folded and z == N:
            continue
        mult = (MULT_EN if en else MULT_ZH)[len(locs)]
        if en:
            word = HW_PREFIX_EN[z]
        else:
            word = (HW_ZH_SHORT.get((n, z)) if short else None) or HW_PREFIX_ZH[z]
        terms.append(elide_a([mult, word]) if mult else word)  # tetra+aza → tetraza
    return terms


def _zh_folded_n(n: int, zs, by_z: dict[int, list[int]], zh_base: str) -> str:
    """折入唑/嗪的氮数量词：三唑、二嗪烷（仅五/六元含氮环的唑嗪系）。"""
    if not (N in zs and n in (5, 6)):
        return zh_base
    count = len(by_z.get(N, ()))
    return f"{MULT_ZH[count]}{zh_base}" if count else zh_base


def hw_name_from_cycle(zs, n_db: int, lambda_at: dict[int, int] | None = None) -> tuple[str, str, str] | None:
    """已编号环序元素表 + 环内双键数 → (英文名, 中文名, 英文位次前缀)。

    传入的 zs 必须是编号后的顺序：第 i 个原子的位次为 i+1。
    lambda_at 为 {1 起位次: 键数}，命中者在位次后带 λn（P-22.2.7.1）。
    """
    if not 3 <= len(zs) <= 10:
        return None
    by_z = locant_map(zs)
    if not by_z:
        return None
    stem = unsat_stem(len(zs), zs) if n_db else sat_stem(len(zs), zs)
    zh_base = zh_stem(len(zs), zs, n_db)
    if stem is None or zh_base is None:
        return None
    zh_base = _zh_folded_n(len(zs), zs, by_z, zh_base)
    if omit_locants(len(zs), zs):
        locs_en = locs_zh = ""
    else:
        locs_en = f"{locant_string_lambda(by_z, lambda_at, en=True)}-"
        locs_zh = f"{locant_string_lambda(by_z, lambda_at, en=False)}-"
    en = f"{locs_en}{elide_a(_hetero_terms(len(zs), zs, by_z, en=True) + [stem])}"
    zh = f"{locs_zh}{''.join(_hetero_terms(len(zs), zs, by_z, en=False))}{zh_base}"
    return en, zh, locs_en


def omit_locants(n: int, zs) -> bool:
    """P-22.2.2.1.7：单杂原子、或元素组成在环上排布唯一时省略全部位次。"""
    heteros = ring_heteros(zs)
    if not heteros:
        return False
    if len(heteros) == 1:
        return True
    return _arrangement_classes(n, tuple(sorted(heteros))) == 1


@lru_cache(maxsize=None)
def _arrangement_classes(n: int, heteros: tuple[int, ...]) -> int:
    """杂原子多重集在 n 元环上的循环排布数（模旋转与翻转）。"""
    k = len(heteros)
    classes = set()
    for spots in combinations(range(1, n), k - 1):
        pos = (0,) + spots
        for perm in set(permutations(heteros)):
            seq = [C] * n
            for p, z in zip(pos, perm):
                seq[p] = z
            rev = seq[::-1]
            classes.add(min(min(tuple(seq[i:] + seq[:i]) for i in range(n)),
                            min(tuple(rev[i:] + rev[:i]) for i in range(n))))
    return len(classes)


# ── 骨架身份与词干入口（L2 接线）──────────────

def is_hw_scaffold(sid: str | None) -> bool:
    """sid 是否为生成式 HW 杂单环骨架。"""
    return sid == HW_ID


def in_scope(zs) -> bool:
    """环内元素是否都在 HW 词表内（词表外元素不得臆造前缀）。"""
    return all(z == C or z in HW_PREFIX_EN for z in zs)


def isolated_ring(mol, atoms) -> bool:
    """原子集是否为不与他环稠合的单环（稠合/桥环另走 P-25/P-23）。"""
    from namepredict.layer1.ring_systems import sssr_rings

    ring = set(atoms)
    rings = sssr_rings(mol)
    if not any(set(r) == ring for r in rings):
        return False
    return all(sum(1 for r in rings if a in r) == 1 for a in ring)


def ring_double_bonds(mol, chain) -> int:
    """环内 Kekulé 双键数（芳香环按 Kekulé 视图计）。"""
    from namepredict.layer1.ring_systems import kekulized

    kek = kekulized(mol) or mol
    atoms = set(chain)
    return sum(1 for b in kek.GetBonds()
               if b.GetBondType() == Chem.BondType.DOUBLE
               and b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms)


def ring_numbering(mol, chain) -> list[int] | None:
    """P-22.2.2.1.3 环编号：复用 L4 窄化引擎，保证名中位次与 L4 编号同源。"""
    from namepredict.layer4.numbering_engine import (
        _narrow_hetero_ring, _ring_cands, _to_chain,
    )

    cands = _narrow_hetero_ring(_ring_cands(list(chain)), mol, list(chain), False)
    return _to_chain(cands[0]) if cands else None


def parent_names(sid: str | None, mol, chain) -> tuple[str, str] | None:
    """生成式 HW 母体双语名（含位次前缀）；非 HW 骨架或词表外元素返回 None。"""
    if not is_hw_scaffold(sid) or mol is None or not chain:
        return None
    ordered = ring_numbering(mol, chain)
    if ordered is None:
        return None
    zs = tuple(mol.GetAtomWithIdx(a).GetAtomicNum() for a in ordered)
    if not in_scope(zs):
        return None
    names = hw_name_from_cycle(zs, ring_double_bonds(mol, ordered), ring_lambda_atoms(ordered, mol))
    return (names[0], names[1]) if names else None


def ring_order(mol, atoms) -> list[int] | None:
    """按环连接走出原子集的环序；原子集不是单环或输入顺序无关。"""
    rset = set(atoms)
    adj = {a: [b.GetOtherAtomIdx(a) for b in mol.GetAtomWithIdx(a).GetBonds()
               if b.GetOtherAtomIdx(a) in rset] for a in rset}
    if any(len(v) != 2 for v in adj.values()):
        return None
    start = min(rset)
    out = [start, adj[start][0]]
    while len(out) < len(rset):
        nxt = [x for x in adj[out[-1]] if x != out[-2]]
        if not nxt:
            return None
        out.append(nxt[0])
    return out


def _max_matching_mask(n: int, capable) -> int:
    """环上最大匹配（边 i 连接 i 与 i+1）；n<=10 直接枚举边子集。"""
    best, best_mask = -1, 0
    for mask in range(1 << n):
        if mask & (mask >> 1) or ((mask & 1) and (mask >> (n - 1)) & 1):
            continue  # 相邻边或首尾边不可同选
        if all(capable[i] and capable[(i + 1) % n] for i in range(n) if mask >> i & 1):
            cnt = bin(mask).count("1")
            if cnt > best:
                best, best_mask = cnt, mask
    return best_mask


def effective_valence(atom) -> int | None:
    """有效价：Table 2.4 键数 + 形式电荷（N+ 视作 4 价、O- 视作 1 价）。"""
    base = HW_MAX_VALENCE.get(atom.GetAtomicNum())
    return None if base is None else base + atom.GetFormalCharge()


def mancude_hydrogens(mol, chain) -> tuple[dict[int, int], int] | None:
    """mancude 参照：({环原子: 参照氢数}, 参照环内双键数)（P-22.2.2.1.1/P-31.2）。"""
    from namepredict.layer1.ring_systems import kekulized

    order = ring_order(mol, chain)
    if order is None:
        return None
    kek = kekulized(mol) or mol
    rset = set(order)
    n = len(order)
    outer: dict[int, float] = {}
    caps: dict[int, int] = {}
    capable = []
    for a in order:
        atom = kek.GetAtomWithIdx(a)
        cap = effective_valence(atom)
        if cap is None:
            return None
        caps[a] = cap
        outer[a] = sum(b.GetBondTypeAsDouble() for b in atom.GetBonds()
                       if b.GetOtherAtomIdx(a) not in rset
                       and kek.GetAtomWithIdx(b.GetOtherAtomIdx(a)).GetAtomicNum() != 1)
        capable.append(outer[a] + 3 <= cap)  # 环内两键 + 一个 π 键不超价
    mask = _max_matching_mask(n, capable)
    ref: dict[int, int] = {}
    for i, a in enumerate(order):
        matched = bool((mask >> i) & 1) or bool((mask >> ((i - 1) % n)) & 1)
        ref[a] = int(caps[a] - (3 if matched else 2) - outer[a])
    return ref, bin(mask).count("1")


def hydro_atoms(mol, chain) -> frozenset[int]:
    """加氢位：实际氢数多于 mancude 参照的环原子（P-31.2.2 / P-54.4.1）。

    环内无重键时取饱和词干，加氢前缀不再使用（否则饱和环会被误标 dihydro）。
    """
    if ring_double_bonds(mol, chain) == 0:
        return frozenset()
    ref = mancude_hydrogens(mol, chain)
    if ref is None:
        return frozenset()
    h_ref, _ = ref
    return frozenset(a for a, h in h_ref.items()
                     if mol.GetAtomWithIdx(a).GetTotalNumHs() > h)


def _parse_key(sid: str) -> tuple[int, ...] | None:
    """'hw:SiCCO' → 环序元素表（按编号顺序）。"""
    body = sid[len(HW_COMPONENT_PREFIX):]
    zs = []
    for m in re.finditer(r"[A-Z][a-z]?", body):
        z = _PT.GetAtomicNumber(m.group(0))
        if z <= 0 or (z != C and z not in HW_PREFIX_EN):
            return None
        zs.append(z)
    return tuple(zs) or None


def _mancude_db(zs) -> int:
    """mancude 环内双键数（纯拓扑最大匹配：价 >= 3 者才可承担 π 键）。"""
    capable = [HW_MAX_VALENCE[z] >= 3 for z in zs]
    return bin(_max_matching_mask(len(zs), capable)).count("1")


def component_key(mol, atoms) -> str | None:
    """稠合组分键 'hw:OCOCC'：仅 3-10 元、含杂、mancude 的孤立环（P-25.2.2.1.1）。"""
    order = ring_order(mol, atoms)
    if order is None or not 3 <= len(order) <= 10:
        return None
    zs = tuple(mol.GetAtomWithIdx(a).GetAtomicNum() for a in order)
    if not ring_heteros(zs) or not in_scope(zs):
        return None
    ref = mancude_hydrogens(mol, order)
    if ref is None or ring_double_bonds(mol, order) != ref[1]:
        return None  # 非 mancude 环不作组分（饱和/部分饱和仍走模板路径）
    ordered = ring_numbering(mol, order)
    if ordered is None:
        return None
    seq = tuple(mol.GetAtomWithIdx(a).GetAtomicNum() for a in ordered)
    return HW_COMPONENT_PREFIX + "".join(_PT.GetElementSymbol(z) for z in seq)


def component_names(sid: str) -> tuple[str, str] | None:
    """'hw:…' → 稠合母体组分词干 (en, zh)（P-25.3.2.1.2）。

    附加组分前缀由 fused_namer 的「去尾 e 加 o」通用式给出（P-25.3.2.2.2）。
    """
    zs = _parse_key(sid)
    if zs is None or not 3 <= len(zs) <= 10:
        return None
    names = hw_name_from_cycle(zs, _mancude_db(zs))
    if names is None:
        return None
    en, zh, locs = names
    if locs:  # 稠合组分位次须加方括号（P-25.3.2.1.2）
        bare = locs[:-1]
        en, zh = f"[{bare}]{en[len(locs):]}", f"[{bare}]{zh[len(locs):]}"
    return en, zh


def identity(info: dict, skeleton):
    """未命中保留模板的 3-10 元孤立含杂单环 → 生成式 HW 骨架身份（P-22.2.2）。"""
    from namepredict.layer2.ring_scaffold import ScaffoldIdentity

    mol = info.get("mol")
    atoms = tuple(skeleton.atom_ids)
    if mol is None or not 3 <= len(atoms) <= 10:
        return None
    zs = tuple(mol.GetAtomWithIdx(a).GetAtomicNum() for a in atoms)
    if not ring_heteros(zs) or not in_scope(zs):
        return None  # 全碳环走通用 carbocycle；词表外元素不接手
    if not isolated_ring(mol, atoms):
        return None
    return ScaffoldIdentity(HW_ID, HW_CLASS, 1, "hetero")
