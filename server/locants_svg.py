"""L4 编号结构图: SMILES → 带 locant 标注的分子 SVG(供 /name/locants-svg)。"""
from __future__ import annotations

from typing import Any

from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from namepredict.layer0.preprocessor import preprocess
from namepredict.layer1.analyzer import analyze
from namepredict.layer2.parent_ownership import finalize_parent_ownership
from namepredict.layer2.parent_selector import select_parent
from namepredict.layer3.substituent_extractor import extract_substituents
from namepredict.layer4.numbering import number

# 主题色是绿(app.css --color-accent),白底 SVG 上用深绿保证可读。
_LOCANT_FILL = "#16a34a"


def _chain_locants(chain: list[int], labels) -> dict[int, str]:
    """映射 {原子: locant};labels 与 chain 同位时优先采用,否则用位置 + 1。"""
    out: dict[int, str] = {}
    for i, atom in enumerate(chain):
        out[atom] = str(labels[i]) if labels is not None else str(i + 1)
    return out


def _draw_with_locants(mol, locants: dict[int, str]) -> str:
    """计算 2D 坐标、绘制结构图并叠加编号文本,返回完整 SVG 字符串。"""
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


def build_locants_svg(smiles: str) -> dict[str, Any] | None:
    """SMILES → 带 L4 编号的结构图 SVG;任何失败返回 None。"""
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
        svg = _draw_with_locants(mol, locants)
    except Exception:
        return None
    return {"svg": svg}
