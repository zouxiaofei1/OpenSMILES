"""L4 编号结构图: SMILES → 带 locant 标注的分子 SVG(供 /name/locants-svg)。"""
from __future__ import annotations

from typing import Any

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

# 主题色是绿(app.css --color-accent),白底 SVG 上用深绿保证可读。
_LOCANT_FILL = "#16a34a"

# preferred_orientation 模板坐标的环边长为 1.0, 缩放到 rdkit 标准 2D 键长 1.5,
# 使稠环水平行摆放与其余原子(取代基等)的 rdkit 布局尺度一致。
_ORIENT_SCALE = 1.5


def _orientation_coords(mol, chain: list[int]) -> dict[int, tuple[float, float]] | None:
    """稠环系统若存在优选取向(P-25.3.2.3), 返回 {原子: (x, y)} 水平行坐标。

    仅对全芳香多环系统应用(与 L4 编号 _fused_numbering 的护栏一致), 否则
    None 表示退回 rdkit 默认 2D 布局。坐标已缩放到 rdkit 键长尺度。
    """
    from namepredict.layer1.ring_systems import build_ring_systems
    from namepredict.layer4.fused_orientation import preferred_orientation

    try:
        chain_set = set(chain)
        systems = [
            s for s in build_ring_systems(mol)
            if (s.get("atom_ids") or []) == sorted(chain_set)
        ]
        if not systems:
            return None
        sys0 = systems[0]
        # 注意：L1 的 _system_dict 不产出 is_aromatic_mancude，故 not None 恒真、
        # 此处恒退。要启用本路径须先在 L1 按 P-21.1 补出该字段。
        if not sys0.get("is_aromatic_mancude") or len(sys0.get("sssr_indices") or []) < 2:
            return None
        rings = list(mol.GetRingInfo().AtomRings())
        orient = preferred_orientation(mol, rings, sys0["fusion_edges"])
        if orient is None:
            return None
        return {
            a: (x * _ORIENT_SCALE, y * _ORIENT_SCALE)
            for a, (x, y) in orient.coord_dict().items()
        }
    except Exception:
        return None


def _chain_locants(chain: list[int], labels) -> dict[int, str]:
    """映射 {原子: locant};labels 与 chain 同位时优先采用,否则用位置 + 1。"""
    out: dict[int, str] = {}
    for i, atom in enumerate(chain):
        out[atom] = str(labels[i]) if labels is not None else str(i + 1)
    return out


def _draw_with_locants(mol, locants: dict[int, str],
                       orient_coords: dict[int, tuple[float, float]] | None = None) -> str:
    """计算 2D 坐标、绘制结构图并叠加编号文本,返回完整 SVG 字符串。

    orient_coords 非空时,稠环原子按优选取向水平行坐标摆放(P-25.3.2.3),
    其余原子(取代基)保留 rdkit 布局;为 None 时全部用 rdkit 默认布局。
    """
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    if orient_coords:
        for atom, (x, y) in orient_coords.items():
            conf.SetAtomPosition(atom, Point3D(x, y, 0))
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    w = max(int(max(xs) - min(xs) + 140), 320)
    h = max(int(max(ys) - min(ys) + 140), 320)
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().padding = 0.15
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    texts = []
    for atom, loc in locants.items():
        x, y = d.GetDrawCoords(atom)
        texts.append(f'<text x="{x + 6:.1f}" y="{y - 8:.1f}">{loc}</text>')
    group = (
        '<g font-family="Arial, Helvetica, sans-serif" font-size="12" '
        'text-anchor="middle" fill="' + _LOCANT_FILL + '" stroke="#ffffff" '
        'stroke-width="2.5" stroke-linejoin="round" paint-order="stroke">'
        + "".join(texts) + "</g>"
    )
    idx = svg.rfind("</svg>")
    if idx == -1:
        raise ValueError("no svg root")
    return svg[:idx] + group + svg[idx:]


def build_locants_svg(smiles: str, orient: bool = True) -> dict[str, Any] | None:
    """SMILES → 带 L4 编号的结构图 SVG;任何失败返回 None。

    orient=True 时稠环按优选取向水平行摆放;False 退回 rdkit 默认布局。
    namepredict 在函数内导入(且置于 try 外): src/ 坏掉时报 500 而非静默返回 None。
    """
    from namepredict.layer0.preprocessor import preprocess
    from namepredict.layer1.analyzer import analyze
    from namepredict.layer2.parent_ownership import finalize_parent_ownership
    from namepredict.layer2.parent_selector import select_parent
    from namepredict.layer3.substituent_extractor import extract_substituents
    from namepredict.layer4.numbering import number

    try:
        mol = preprocess(smiles)
        if mol is None:
            return None
        info = analyze(mol)
        selected = select_parent(info)
        if selected is None:
            return None
        parent = finalize_parent_ownership(selected, mol)
        subst = extract_substituents(info, parent)
        chain_set = set(parent.get("chain") or [])
        subs_numbered = [s for s in subst if s.get("attach_idx") in chain_set]
        numbered = number(parent, subs_numbered)
    except Exception:
        return None

    p = numbered.get("parent") or {}
    chain = p.get("chain") or []
    if not chain:
        return None

    labels = None
    scaf = p.get("numbering_scaffold")
    if scaf:
        lab = scaf.get("labels")
        if lab and len(lab) == len(chain):
            labels = lab

    locants = _chain_locants(chain, labels)
    try:
        orient_coords = _orientation_coords(mol, chain) if orient else None
        svg = _draw_with_locants(mol, locants, orient_coords)
    except Exception:
        return None
    return {"svg": svg}
