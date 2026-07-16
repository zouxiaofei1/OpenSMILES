"""L4 orientation and assembled stem facts for unsaturated polycarboxylic parents."""
from __future__ import annotations
from rdkit.Chem import BondStereo
from namepredict.layer4.polyene import ene_locants

_ALKANE_NAMES = {
    1: ("methane", "甲烷"), 2: ("ethane", "乙烷"), 3: ("propane", "丙烷"),
    4: ("butane", "丁烷"), 5: ("pentane", "戊烷"), 6: ("hexane", "己烷"),
    7: ("heptane", "庚烷"), 8: ("octane", "辛烷"), 9: ("nonane", "壬烷"),
    10: ("decane", "癸烷"),
}
_MULT_EN = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
_MULT_ZH = {1: "", 2: "二", 3: "三", 4: "四", 5: "五", 6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}

def _bond_locs(chain: list[int], bonds) -> list[int]:
    return sorted(min(chain.index(a) + 1, chain.index(b) + 1) for a, b in bonds)

def unsat_locant_key(chain: list[int], parent: dict) -> tuple[int, ...]:
    return tuple(sorted([*_bond_locs(chain, parent.get("double_bonds") or []), *_bond_locs(chain, parent.get("triple_bonds") or [])]))

def orient_polycarboxylic(chain: list[int], parent: dict, fallback) -> list[int]:
    """Prefer the complete ene/yne locant set before carboxyl locants."""
    reverse = list(reversed(chain))
    forward_key, reverse_key = unsat_locant_key(chain, parent), unsat_locant_key(reverse, parent)
    if forward_key != reverse_key:
        return chain if forward_key < reverse_key else reverse
    return fallback(chain, parent, "cooh_c_idxs", [])

def _ez_prefix(parent: dict) -> str:
    mol, chain = parent.get("mol"), parent.get("chain") or []
    parts = []
    for a, b in parent.get("double_bonds") or []:
        bond = None if mol is None else mol.GetBondBetweenAtoms(a, b)
        tag = None if bond is None else bond.GetStereo()
        letter = "E" if tag == BondStereo.STEREOE else "Z" if tag == BondStereo.STEREOZ else ""
        if letter: parts.append((min(chain.index(a) + 1, chain.index(b) + 1), letter))
    return f"({','.join(f'{loc}{tag}' for loc, tag in sorted(parts))})-" if parts else ""

def _en_stem(base: str, enes: list[int], ynes: list[int], ez: str) -> str:
    stem = base[:-3]
    ene = f"{','.join(map(str, enes))}-{_MULT_EN.get(len(enes), '')}ene" if enes else ""
    yne = f"{','.join(map(str, ynes))}-{_MULT_EN.get(len(ynes), '')}yne" if ynes else ""
    if enes and len(enes) > 1: stem += "a"
    return f"{ez}{stem}-{'-'.join(filter(None, (ene, yne)))}" if (ene or yne) else base

def _zh_stem(base: str, enes: list[int], ynes: list[int], ez: str) -> str:
    parts = [f"{','.join(map(str, enes))}-{_MULT_ZH.get(len(enes), '')}烯" if enes else "", f"{','.join(map(str, ynes))}-{_MULT_ZH.get(len(ynes), '')}炔" if ynes else ""]
    suffix = "-".join(filter(None, parts))
    return f"{ez}{base[:-1] if base.endswith('烷') else base}-{suffix}" if suffix else base

def polycarboxylic_facts(oriented: dict) -> dict:
    """Materialize assembly-only stem and unsaturation facts after orientation."""
    base = _ALKANE_NAMES.get(oriented.get("n_carbons", 0))
    if not base: return {}
    chain = oriented.get("chain") or []
    enes, ynes = _bond_locs(chain, oriented.get("double_bonds") or []), _bond_locs(chain, oriented.get("triple_bonds") or [])
    ez = _ez_prefix(oriented)
    return {"ene_locants": enes, "yne_locants": ynes, "stereo_prefix": ez, "stem_en": _en_stem(base[0], enes, ynes, ez), "stem_zh": _zh_stem(base[1], enes, ynes, ez)}
