"""L4 定向编号引擎：按 P-14.4 规则筛出链/环原子顺序候选。"""
from __future__ import annotations

from namepredict.tools import memo
from namepredict.layer1.ring_systems import sssr_rings

# P-25.3.3：这些保留骨架按传统编号，不走 P-25.3.3.1 的优选取向自动编号。
# xanthene 及其硫属类似物（xanthene/thioxanthene）与 cyclopenta[a]phenanthrene（甾体 1-17）已按
# _TEMPLATES 的 standard 字段登记传统编号，故一并列入。
_TRADITIONAL_NUMBERING_IDS = frozenset({
    "anthracene", "phenanthrene", "acridine", "carbazole", "purine",
    "xanthene", "thioxanthene", "cyclopenta[a]phenanthrene",
})


# ── 候选生成 ──────────────────────────────────────────────────

def _numbered(chain: list[int]) -> dict[int, int]:
    """把链原子顺序映射为 {原子: 位次}。"""
    return {a: i + 1 for i, a in enumerate(chain)}


def _to_chain(cand: dict[int, int]) -> list[int]:
    """按位次升序还原原子顺序列表。"""
    return sorted(cand, key=cand.get)


def _chain_cands(chain: list[int]) -> list[dict[int, int]]:
    """生成线性链正反两个方向的编号候选。"""
    return [_numbered(chain), _numbered(list(reversed(chain)))]


def _ring_cands(chain: list[int]) -> list[dict[int, int]]:
    """生成环的全部旋转/翻转编号候选。"""
    n = len(chain)
    out: list[dict[int, int]] = []
    for i in range(n):
        fwd = chain[i:] + chain[:i]
        out.append(_numbered(fwd))
        out.append(_numbered(list(reversed(fwd))))
    return out


# ── 位次集合键 ───────────────────────────────────────────────────────

def _locant_set(cand: dict[int, int], atoms: list[int]) -> tuple[int, ...] | None:
    """计算原子集合在候选编号下的排序位次元组。"""
    locs = sorted(cand[a] for a in atoms if a in cand)
    return tuple(locs) if locs else None


def _bond_locants(cand: dict[int, int], bonds) -> tuple[int, ...] | None:
    """计算各多重键沿编号方向占据的边位置（seam 感知，每条键一个数）。"""
    if not bonds:
        return None
    order = sorted(cand, key=cand.get)  # locant 升序 → 环/链遍历顺序
    n = len(order)
    pos = {a: i for i, a in enumerate(order)}
    locs = []
    for a, b in bonds:
        if a not in pos or b not in pos:
            continue
        ia, ib = pos[a], pos[b]
        if ia > ib:
            ia, ib = ib, ia
        if ib - ia == 1:
            locs.append(ia + 1)              # 普通相邻边
        elif ia == 0 and ib == n - 1:
            locs.append(n)                   # seam 闭合边：跨编号首尾，记 n
        else:
            return None                      # 端点不沿编号相邻：判据不适用
    return tuple(sorted(locs)) if len(locs) == len(bonds) else None  # print(cand, bonds)


def _narrow(cands: list[dict], key_fn) -> list[dict]:
    """保留位次集合键最小的候选；得到一个即提前停止。"""
    if len(cands) <= 1:
        return cands
    keys = [key_fn(c) for c in cands]
    if any(k is None for k in keys):
        return cands  # 特征全部缺失 → 规则不适用
    best = min(keys)
    return [c for c, k in zip(cands, keys) if k == best]


# ── P-14.4(j)：CIP 立体描述符平局破（低位次给 R/M/r，优先于 S/P/s）─────

_RS_HI = frozenset({"R", "M", "r"})   # 编号优先级高的 CIP 描述符（P-91.2）
_RS_LO = frozenset({"S", "P", "s"})   # 与之成对、取较低位次的次位描述符


def assign_cip(mol) -> None:
    """强制重算分子 CIP（隐式 H 手性碳先补显式 H）；编号期与 L5 stereo 打印共用此唯一实现。

    对同一个 mol 在一次命名内只算一次：结果是 mol 的确定函数，重复调用只是重写同样的
    `_CIPCode`。调用点分布在不同层，键仍按 mol 对象身份，故记忆放在本函数内而非各调用点。
    """
    memo.by_mol("cip", _assign_cip_uncached, mol)


def _assign_cip_uncached(mol) -> None:
    """实际执行 CIP 重算（无记忆版本，见 assign_cip）。"""
    from rdkit import Chem
    from rdkit.Chem import ChiralType, rdCIPLabeler

    r = Chem.RWMol(mol)
    for a in list(r.GetAtoms()):
        if a.GetChiralTag() != ChiralType.CHI_UNSPECIFIED and a.GetTotalNumHs() == 0 and a.GetDegree() < 4:
            r.AddBond(a.GetIdx(), r.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
    m = r.GetMol()
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)
    rdCIPLabeler.AssignCIPLabels(m)
    for a in mol.GetAtoms():
        if a.HasProp("_CIPCode"):
            a.ClearProp("_CIPCode")
        b = m.GetAtomWithIdx(a.GetIdx())
        if b.HasProp("_CIPCode"):
            a.SetProp("_CIPCode", b.GetProp("_CIPCode"))


def _chain_rs_codes(mol, chain: list[int]) -> dict[int, str]:
    """取母体链原子上的 CIP 代码映射；键上 E/Z 与无双键歧义不入此表（(j) 原子级破局只需 R/S）。"""
    if mol is None or not chain:
        return {}
    assign_cip(mol)
    codes: dict[int, str] = {}
    for idx in chain:
        a = mol.GetAtomWithIdx(int(idx))
        if a.HasProp("_CIPCode") and a.GetProp("_CIPCode") in _RS_HI | _RS_LO:
            codes[int(idx)] = a.GetProp("_CIPCode")
    return codes


def _rs_locant_key(codes: dict[int, str], chain: list[int]) -> tuple:
    """P-14.4(j) 排序键：R/M/r 描述符位次升序在前、S/P/s 位次在后，字典序小者取低位次。"""
    hi, lo = [], []
    for loc, idx in enumerate(chain, 1):
        code = codes.get(int(idx))
        if code in _RS_HI:
            hi.append(loc)
        elif code in _RS_LO:
            lo.append(loc)
    return tuple(hi), tuple(lo)


# ── 从 parent dict 提取 P-14.4 特征 ────────────────────────

def _principal_atoms(parent: dict) -> list[int]:
    """P-14.4(c)：principal characteristic group 的附着原子。"""
    atoms = parent.get("principal_attachment_atoms")
    if atoms:
        return list(atoms)
    facts = parent.get("principal_expression_facts")
    if facts is not None:
        attach = getattr(facts, "attachment_atoms", None)
        if attach:
            return sorted(attach)
    return []


def _unsat_bonds(parent: dict) -> tuple[list, list]:
    """（全部多重键、双键）端点对：单数 double_bond 亦进 all_bonds（混合烯炔整体最小化）。"""
    all_bonds, doubles = [], []
    for key in ("double_bond", "triple_bond"):
        v = parent.get(key)
        if v:
            pair = (v[0], v[1])
            all_bonds.append(pair)
            if key == "double_bond":
                doubles.append(pair)
    for key in ("double_bonds", "triple_bonds"):
        for b in parent.get(key) or []:
            pair = (b[0], b[1])
            all_bonds.append(pair)
            if key == "double_bonds":
                doubles.append(pair)
    return all_bonds, doubles


def _is_ring(parent: dict) -> bool:
    """按 scaffold_id 判断 parent 是否为环系。"""
    return bool(parent.get("scaffold_id"))


# P-14.4(a)：parent dict 中标定必须为 locant 1 的原子的字段（环外羰基连接、自由基中心）。
# 杂环起点已改由 _narrow_hetero_ring 按元素序决定，不在此列。
_FIXED_START_KEYS = (
    "ring_attach_idx", "n_idx", "nh_idx", "hetero_idx", "radical_c_idx", "acyl_c_idx",
)


def _fixed_start(parent: dict) -> int | None:
    """取碳环/链固定 locant 1 起点原子（P-14.4(a) FG 锚/自由基字段）；杂环不落入此函数。"""
    for key in _FIXED_START_KEYS:
        v = parent.get(key)
        if v is not None:
            return v
    return None


def _narrow_hetero_ring(cands: list[dict], mol, chain: list[int], float_hetero: bool) -> list[dict]:
    """杂环编号 P-22.2.2.1.3/P-25.3.3.1.2(b)：(a) 全杂原子集最低位次→(b) 按 F>…>O>S>…>N>… 元素序逐元素收窄→(c) 同元素 N 中带 H/3 价取代者得低位（唑 NH=1）；float_hetero（稠合组分）跳过 (c)。"""
    from namepredict.constants import P145_SENIOR
    heteros = [a for a in chain if mol.GetAtomWithIdx(a).GetAtomicNum() != 6]
    cands = _narrow(cands, lambda c: _locant_set(c, heteros))            # (a)
    by_z: dict[int, list[int]] = {}
    for a in heteros:
        by_z.setdefault(mol.GetAtomWithIdx(a).GetAtomicNum(), []).append(a)
    for z in P145_SENIOR:                                               # (b)
        atoms = by_z.get(z)
        if atoms:
            cands = _narrow(cands, lambda c, at=sorted(atoms): _locant_set(c, at))
    if not float_hetero:                                                # (c)
        n_active = [a for a in heteros
                    if mol.GetAtomWithIdx(a).GetAtomicNum() == 7
                    and (mol.GetAtomWithIdx(a).GetTotalNumHs() > 0
                         or mol.GetAtomWithIdx(a).GetDegree() == 3)]
        if n_active:
            cands = _narrow(cands, lambda c: _locant_set(c, sorted(n_active)))
    return cands


def _template_matches(mol, sid: str, atoms: frozenset) -> list:
    """模板 atom_ids 全覆盖环系的全部匹配：先直接匹配，失败则退到氢化骨架匹配（原子序不变，映射可回原分子）。"""
    from namepredict.layer2.ring_scaffold import _Q, _Q_H, _hydrogenated
    q = _Q.get(sid)
    matches = [m for m in mol.GetSubstructMatches(q, uniquify=False) if set(m) == atoms]
    if matches:
        return matches
    qh = _Q_H.get(sid)
    mol_h = memo.by_mol("hydrogenated", _hydrogenated, mol)
    if qh is None or mol_h is None:
        return []
    return [m for m in mol_h.GetSubstructMatches(qh, uniquify=False) if set(m) == atoms]


def _fixed_numbering(parent: dict, chain: list[int], substituents: list | None = None) -> list[int] | None:
    """P-14.4(a)：fused 环经模板 standard_path 映射固定编号（P-25.4 起点方向不随取代基变，避免如喹啉 10-氯 vs 2-氯 算错）；对称 scaffold 用全自同构等价链按取代基位次最小化。"""
    sid = parent.get("scaffold_id")
    mol = parent.get("mol")
    if not sid or mol is None:
        return None
    from namepredict.layer2.ring_scaffold import _Q, standard_chain
    if _Q.get(sid) is None:
        return None
    atoms = frozenset(chain)
    chains = []
    for m in _template_matches(mol, sid, atoms):
        std = standard_chain(sid, tuple(m))
        if std is not None and set(std) == atoms:
            chains.append(std)
    if not chains:
        return None
    if len(chains) == 1:
        return chains[0]
    subs = [s["attach_idx"] for s in (substituents or []) if s.get("attach_idx") in chain]
    if parent.get("radical_c_idx") in chain:  # 自由基主基团：自由价连接点（radical_c_idx）按 P-14.4 最低化到 locant 1，避免咔唑类对称 scaffold 取镜像首方向把自由价标成 8（应 1）
        subs.append(parent["radical_c_idx"])
    if not subs:
        return chains[0]
    return min(chains, key=lambda std: tuple(sorted(std.index(a) + 1 for a in subs)))


def _fused_numbering(parent: dict, chain: list[int],
                     substituents: list | None = None) -> list[int] | None:
    """P-25.3.3 稠环编号（护栏：P-25.3.3 传统编号例外骨架保持固定编号，其余芳香或未注册稠环走优选取向+外周编号）。"""
    sid = parent.get("scaffold_id")
    if sid in _TRADITIONAL_NUMBERING_IDS:
        return None
    mol = parent.get("mol")
    from namepredict.layer2.ring_scaffold import get_spec
    spec = get_spec(parent.get("scaffold_id") or "")
    registered_fused = bool(spec and spec.n_rings >= 2)  # 稠环保留模板（含氢化衍生物）：编号随骨架拓扑，与芳香性无关（fused_numbering 本身不依赖 aromatic）
    chain_arom = bool(mol is not None and chain) and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in chain)
    # 未注册稠环（如 fused_hetero 身份的饱和/部分饱和稠环）：scaffold_id 不在 _TEMPLATES，registered_fused 与 chain_arom 均假。
    # 若在此拦下会回落到 P-14.4 单环通用枚举：桥头得数字位而非字母位(4a/8a)、起点方向随输入原子序漂移。
    # 故按「chain 覆盖 ≥2 环」放行走 P-25.3.3 外周编号；单环/开链覆盖 ≤1 环仍返回 None，保持原路径不变。
    chain_fused = bool(mol is not None and chain) and sum(
        1 for r in sssr_rings(mol) if set(r) <= set(chain)) >= 2
    if mol is None or not chain or not (chain_arom or registered_fused or chain_fused):
        return None
    from namepredict.layer1.ring_systems import build_ring_systems
    from namepredict.layer4.fused_orientation import preferred_orientations
    from namepredict.layer4.fused_numbering import number_fused_system
    systems = [s for s in build_ring_systems(mol) if (s.get("atom_ids") or []) == sorted(set(chain))]
    if not systems:
        return None
    system = systems[0]
    sssr = list(system.get("sssr_indices") or [])
    if len(sssr) < 2:
        return None  # 单环走 P-14.4 通用枚举
    atom_rings = list(sssr_rings(mol))  # 仅取稠合系统自身环（全分子 AtomRings 会把取代基无关环算进 layout 致 orientation 失败，如 chebi-300）；rings 过滤后索引重排，fusion_edges 的 SSSR 索引须同步重映射，否则 horizontal_rows KeyError 退化成只命名侧链。
    rings = [atom_rings[i] for i in sssr]
    idx_map = {orig: new for new, orig in enumerate(sssr)}
    fusion_edges = [(idx_map[i], idx_map[j], sh) for i, j, sh in (system.get("fusion_edges") or [])
                    if i in idx_map and j in idx_map]
    orients = preferred_orientations(mol, rings, fusion_edges)  # print("riings",rings) / print(orients)
    if not orients:
        return None
    chain_set = set(chain)  # 环外附着原子按优先级分层：纯碳环上准则 (a)-(d) 全平局，靠它逐层收窄镜像取向
    layers: list[list[int]] = []
    radical = parent.get("radical_c_idx")
    if radical in chain_set:
        layers.append([radical])  # P-29: 取代基游离价连接点优先得最低位次
    principal_atoms = sorted(a for a in (parent.get("principal_attachment_atoms") or [])
                             if a in chain_set)
    if principal_atoms:
        layers.append(principal_atoms)  # P-14.4(c): principal 特征基团优先于取代基
    sub_atoms = sorted(s["attach_idx"] for s in (substituents or [])
                       if s.get("attach_idx") in chain_set)
    if sub_atoms:
        layers.append(sub_atoms)  # P-14.4(f): 取代基位次集合最小化
    hydro_atoms = sorted(a for a in (parent.get("hydro_atoms") or ()) if a in chain_set)
    if hydro_atoms:
        layers.append(hydro_atoms)  # P-31.2.2: hydro 加氢位次最低（取代基之后、指示氢之前）
    from namepredict.layer3.substituent_extractor import alkyl_alpha_key
    alpha_subs = [(alkyl_alpha_key(s.get("en") or ""), s["attach_idx"])
                  for s in (substituents or []) if s.get("attach_idx") in chain_set]
    result = number_fused_system(mol, rings, [o.coord_dict() for o in orients], layers, alpha_subs)
    if result is None:
        return None
    fused_chain, labels = result
    parent["numbering_scaffold"] = {
        "scaffold_id": "fused", "labels": tuple(labels), "relative_stereo": None,
    }
    return fused_chain #外环原子顺序


# ── 入口 ─────────────────────────────────────────────────────────────────

def orient_numbering(parent: dict, substituents: list, *, float_hetero: bool = False) -> list[int] | None:
    """返回 P-14.4 定向后的原子顺序；不适用返回 None（杂环走元素序窄化；float_hetero 见 _narrow_hetero_ring，仅稠合组分路径开启）。"""
    chain = parent.get("chain") or []
    if not chain:
        return None
    fixed = _fixed_numbering(parent, chain, substituents)
    if fixed is not None:
        return fixed
    fused = _fused_numbering(parent, chain, substituents)
    if fused is not None:
        return fused
    mol = parent.get("mol")
    anchor_as_principal = None  # 固定 locant 1 失败（链中部自由基/锚点）时把该原子并入 principal 竞争最低位次；杂环分支无此回退
    if mol is not None and _is_ring(parent) and any(
            mol.GetAtomWithIdx(a).GetAtomicNum() != 6 for a in chain):
        cands = _narrow_hetero_ring(_ring_cands(chain), mol, chain, float_hetero)  # 杂环：P-22.2.2.1.3 元素序窄化先于 principal（唑类 N 必须 1,3/1,2、吡啶甲酸 N=1 后由 principal 定方向），不走 FG 锚点/自由基字段，避免醛基环碳抢占 locant 1
    else:
        cands = _ring_cands(chain) if _is_ring(parent) else _chain_cands(chain)  # 碳环/链：P-14.4(a) 固定 locant 1（FG 锚/自由基字段）锚定后退化
        start = _fixed_start(parent)
        if start is not None:
            forced = [c for c in cands if c.get(start) == 1]
            if forced:
                cands = forced
            else:
                anchor_as_principal = start  # 固定起点在链候选中不可能为 locant 1（如链中部自由基/锚点）：不原样保留链序（否则自由价/双键/取代基位次全不最小化、随上游原子序漂移），回退全候选并并入 P-14.4(c) principal 竞争最低位次
    principal = _principal_atoms(parent)
    if anchor_as_principal is not None:
        principal = sorted(set(principal) | {anchor_as_principal})
    if principal:
        cands = _narrow(cands, lambda c: _locant_set(c, principal))
    bonds, doubles = _unsat_bonds(parent)
    if bonds:
        cands = _narrow(cands, lambda c: (_bond_locants(c, bonds), _bond_locants(c, doubles)))
    subs = [s["attach_idx"] for s in substituents if s.get("attach_idx") in chain]

    if subs:
        cands = _narrow(cands, lambda c: _locant_set(c, subs))
    if len(cands) > 1 and substituents:
        from namepredict.layer4._chain_orient import _stem_loc_pairs  # P-14.4(f) 平局：最低位次集合已相同 → 把最低位次给字母序最前的取代基（stem-alpha 对，P-14.5）
        cands = _narrow(cands, lambda c: _stem_loc_pairs(_to_chain(c), substituents))
    if len(cands) > 1:
        codes = _chain_rs_codes(mol, chain)  # P-14.4(j) 立体平局：镜像/反向等价编号仍等优时按 CIP 立体描述符定方向——较低位次赋予 R/M/r（优先于对应的 S/P/s），避免随输入原子序漂移（如内消旋环 1R/7S vs 1S/7R）
        if codes:
            cands = _narrow(cands, lambda c: _rs_locant_key(codes, _to_chain(c)))
    return _to_chain(cands[0])


def _component_labels(parent: dict, chain: list[int]) -> list[str]:
    """orient_numbering 只返回链, 补并行 labels: 优先标准 locant 标签, 否则纯数字。"""
    scaffold = parent.get("numbering_scaffold")
    labels = scaffold.get("labels") if scaffold else None
    if labels and len(labels) == len(chain):
        return list(labels)
    from namepredict.layer2.ring_scaffold import _STANDARD_LABELS
    labels = _STANDARD_LABELS.get(parent.get("scaffold_id") or "")
    if labels and len(labels) == len(chain):
        return list(labels)
    return [str(i + 1) for i in range(len(chain))]


def fused_component_numbering(mol, scaffold_id, sub_rings, shared=None, sub_edges=None):
    """稠合组分自身编号（P-25.4/P-25.3.3），稠合点 shared 作取代基最小化位次；返回 (chain, labels) 或 (None, None)。"""
    if not sub_rings:
        return None, None
    subs = [{"attach_idx": a} for a in (shared or ())] if shared else []
    chain0 = sorted(set().union(*sub_rings))
    from namepredict.layer2.ring_scaffold import _STANDARD_ORDERS
    if len(sub_rings) == 1 or (scaffold_id and scaffold_id in _STANDARD_ORDERS):
        if not scaffold_id:
            ring = list(sub_rings[0])
            return ring, [str(i + 1) for i in range(len(ring))]
        parent = {"mol": mol, "scaffold_id": scaffold_id, "chain": chain0}
        res = orient_numbering(parent, subs, float_hetero=bool(shared))  # float_hetero: 对称等价杂环(嘧啶双 N 等)经元素序窄化后保留两镜像方向，locant 1 交给稠合原子位次最小化决定，使碱环取向与规范稠合字母(d 侧)一致（见 _narrow_hetero_ring 的 (c)-skip）
        if not res:
            return None, None
        return res, _component_labels(parent, res)
    from namepredict.layer4.fused_orientation import preferred_orientations  # 多环无固定编号: P-25.3.3 通用编号(fused_numbering)
    from namepredict.layer4.fused_numbering import number_fused_system
    orients = preferred_orientations(mol, sub_rings, sub_edges or [])
    if not orients:
        return None, None
    result = number_fused_system(mol, sub_rings, [o.coord_dict() for o in orients])
    if result is None:  # print(result)
        return None, None
    return result[0], result[1]
