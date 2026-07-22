# tests/unit/test_recursive_benzamide_north_star.py
# IUPAC: P-66.1 / P-29 / P-25 / P-65
# Layer: L2–L5
from namepredict.constants import normalize_en, normalize_zh
from namepredict.namer import SMILESNNamer

NORTH = (
    "COc1ccc(C(=O)Nc2nc3c(s2)CCC3C(=O)N2CCCC2)cc1",
    "4-methoxy-N-(4-(pyrrolidine-1-carbonyl)-5,6-dihydro-4H-cyclopenta[d]thiazol-2-yl)benzamide",
    "4-甲氧基-N-(4-(吡咯烷-1-羰基)-5,6-二氢-4H-环戊并[d]噻唑-2-基)苯甲酰胺",
)
