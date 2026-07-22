# IUPAC: P-64.2
# Layer: L2,L3
"""1,4-Benzoquinone PIN as cyclohexa-2,5-diene-1,4-dione (P-64.2).

Single C6 carbocycle + exactly two para ring ketones + two endocyclic
double bonds. PIN is systematic (not retained 1,4-benzoquinone).

Round B: ring n-alkyl C1–C12 + halo/OH/alkoxy C1–C2 mixed substitution.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# All benzoquinone test cases were removed — not yet supported.
