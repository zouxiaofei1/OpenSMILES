from __future__ import annotations
import re

# ── Atomic Numbers ──────────────────────────────────────────────
H  = 1
Li = 3
B  = 5
C  = 6
N  = 7
O  = 8
F  = 9
Na = 11
Al = 13
Si = 14
P  = 15
S  = 16
Cl = 17
K  = 19
Ga = 31
Ge = 32
As = 33
Se = 34
Br = 35
In = 49
Sn = 50
Sb = 51
Te = 52
I  = 53
Tl = 81
Pb = 82
Bi = 83

# ── Common Sets ─────────────────────────────────────────────────
HALO_Z = frozenset({F, Cl, Br, I})
HALO_EN = {F: "fluoro", Cl: "chloro", Br: "bromo", I: "iodo"}
HALO_ZH = {F: "氟", Cl: "氯", Br: "溴", I: "碘"}

# ── Multiplicity Prefixes ────────────────────────────────────────
MULT_EN = {1: "", 2: "di", 3: "tri", 4: "tetra", 5: "penta",
           6: "hexa", 7: "hepta", 8: "octa", 9: "nona", 10: "deca"}
MULT_ZH = {1: "", 2: "二", 3: "三", 4: "四", 5: "五",
           6: "六", 7: "七", 8: "八", 9: "九", 10: "十"}

# ── Text Normalization ──────────────────────────────────────────
_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",").replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    return (name or "").strip()
