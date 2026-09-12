# IUPAC: P-65.3 / architecture layer purity
# Layer: L2,L5
"""L5 special FG name modules must not import L2 private APIs.

Aryl / cycloalkyl stems are precomputed into parent side dicts at L2 pack time.
L5 only reads side['en']/side['zh'] (and existing kind/mode fields).
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

from namepredict.tools.re import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

_ROOT = Path(__file__).resolve().parents[2]
_L5 = _ROOT / "src" / "namepredict" / "layer5"

# L5 special-FG name modules that previously walked via L2 private APIs.


def _layer2_import_hits(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            if node.module == "namepredict.layer2" or node.module.startswith(
                "namepredict.layer2."
            ):
                hits.append(f"from {node.module} (L{node.lineno})")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "namepredict.layer2" or alias.name.startswith(
                    "namepredict.layer2."
                ):
                    hits.append(f"import {alias.name} (L{node.lineno})")
    return hits


@pytest.mark.parametrize("path", sorted(_L5.glob("*.py")), ids=lambda p: p.name)
def test_all_l5_modules_have_no_layer2_import(path: Path) -> None:
    hits = _layer2_import_hits(path)
    assert hits == [], f"{path.name} must not import layer2: {hits}"


# End-to-end bilingual gold (behavior must stay flat after cut).
_E2E = [
    # sulfonamide aryl
    ("c1ccc(S(=O)(=O)N)cc1", "benzenesulfonamide", "苯磺酰胺"),
    # urea aryl
    ("c1ccc(NC(=O)N)cc1", "phenylurea", "苯基脲"),
    # carbonate aryl
    (
        "C(OC)(OC1=C(C=C(C=C1)C)C)=O",
        "methyl 2,4-dimethylphenyl carbonate",
        "2,4-二甲苯基甲基碳酸酯",
    ),
    # alkyl sulfonamide still correct (non-aryl path)
    ("CS(=O)(=O)N", "methanesulfonamide", "甲磺酰胺"),
]
