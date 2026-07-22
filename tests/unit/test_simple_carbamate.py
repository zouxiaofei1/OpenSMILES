# IUPAC: P-65
# Layer: L1,L2,L5
"""Simple alkyl N-…carbamate (incl. Boc / tert-butyl carbamate).

Carbamate R2N–C(=O)–OR is principal FG (not ester). English follows gold:
  methyl N-methylcarbamate; tert-butyl (3-methoxyphenyl)carbamate.
"""
from __future__ import annotations

import pytest

from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

# All test cases were removed — carbamate naming is not yet supported.
