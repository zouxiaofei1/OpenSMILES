"""Atom ID + SSSR ring ID 结构图: SMILES → 标注 SVG(供 /name/atom-ids-svg)。"""
from __future__ import annotations

from typing import Any

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# 原子索引红色、SSSR 环徽章蓝色,白底 SVG 上可读。
_ATOM_FILL = "#ef4444"
_RING_FILL = "#0284c7"


def _draw_with_ids(mol) -> str:
    """计算 2D 坐标、绘制结构图并叠加原子/环索引,返回完整 SVG 字符串。"""
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    w = max(int(max(xs) - min(xs) + 140), 320)
    h = max(int(max(ys) - min(ys) + 140), 320)
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().padding = 0.15
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()

    atom_texts = []
    for i in range(mol.GetNumAtoms()):
        x, y = d.GetDrawCoords(i)
        atom_texts.append(f'<text x="{x + 5:.1f}" y="{y - 5:.1f}">{i}</text>')
    atom_group = (
        '<g font-family="Arial, Helvetica, sans-serif" font-size="11" '
        'text-anchor="middle" fill="' + _ATOM_FILL + '" stroke="#ffffff" '
        'stroke-width="2" stroke-linejoin="round" paint-order="stroke">'
        + "".join(atom_texts) + "</g>"
    )

    ring_texts = []
    for k, ring in enumerate(mol.GetRingInfo().AtomRings(), start=1):
        pts = [d.GetDrawCoords(a) for a in ring]
        cx = sum(p.x for p in pts) / len(pts)
        cy = sum(p.y for p in pts) / len(pts)
        ring_texts.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="10" fill="{_RING_FILL}" '
            f'stroke="#ffffff" stroke-width="1"/>'
            f'<text x="{cx:.1f}" y="{cy + 3.5:.1f}">R{k}</text>'
        )
    ring_group = (
        '<g font-family="Arial, Helvetica, sans-serif" font-size="10" '
        'font-weight="bold" text-anchor="middle" fill="#ffffff">'
        + "".join(ring_texts) + "</g>"
    )

    idx = svg.rfind("</svg>")
    if idx == -1:
        raise ValueError("no svg root")
    return svg[:idx] + atom_group + ring_group + svg[idx:]


def build_atom_ids_svg(smiles: str) -> dict[str, Any] | None:
    """SMILES → 标注原子/环索引的结构图 SVG;任何失败返回 None。

    namepredict 在函数内导入(且置于 try 外): src/ 坏掉时报 500 而非静默返回 None。
    """
    from namepredict.layer0.preprocessor import preprocess

    try:
        mol = preprocess(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return None
        svg = _draw_with_ids(mol)
    except Exception:
        return None
    return {"svg": svg}
