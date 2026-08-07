"""L5 names for sulfur functional parents.

Covers sulfonamide / sulfonate / sulfone / sulfonic acid / sulfonyl chloride /
sulfoxide (P-65.3, P-63.3). Sides are precomputed into parent dicts at L2 pack
time; L5 only reads side['en']/side['zh'] (and existing kind/mode fields).
"""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, to_benzene, wrap_aryl
from namepredict.layer5.stems import (
    SULFIDE_ALKYL_EN,
    SULFIDE_ALKYL_ZH,
    ester_alkyl_en,
    ester_alkyl_zh,
)

# ---------------------------------------------------------------------------
# sulfonamide (P-65.3)
# ---------------------------------------------------------------------------

_ALKYL_SAM_EN = {1: "methanesulfonamide", 2: "ethanesulfonamide",
                 3: "propanesulfonamide", 4: "butanesulfonamide"}
_ALKYL_SAM_ZH = {1: "甲磺酰胺", 2: "乙磺酰胺", 3: "丙磺酰胺", 4: "丁磺酰胺"}


def _sa_alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SAM_EN.get(n), _ALKYL_SAM_ZH.get(n)
    return (en, zh) if en and zh else None


def _sa_alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _sa_alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _sa_aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = to_benzene(*side_en_zh(s))
    return f"{en}sulfonamide", f"{zh}磺酰胺"


def _n_aryl_alkyl(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    stem = _sa_alkyl_stem(int(s.get("n") or 0))
    if stem is None or n.get("kind") != "aryl":
        return None
    ar_en, ar_zh = wrap_aryl(*side_en_zh(n))
    return f"N-{ar_en}{stem[0]}", f"N-{ar_zh}{stem[1]}"


def _split_pref(en: str, zh: str) -> tuple[str, str, str, str]:
    """Split '2,4-dichlorobenzene' → ('2,4-dichloro', 'benzene', ...)."""
    if en.endswith("benzene"):
        return en[: -len("benzene")], "benzene", zh[: -len("苯")], "苯"
    return "", en, "", zh


def _with_n_pref(pe: str, pz: str, cy_en: str, cy_zh: str, be: str, bz: str):
    head_en = f"{pe}-N-" if pe else "N-"
    head_zh = f"{pz}-N-" if pz else "N-"
    return f"{head_en}{cy_en}{be}sulfonamide", f"{head_zh}{cy_zh}{bz}磺酰胺"


def _n_cyclo_aryl(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    if s.get("kind") != "aryl" or n.get("kind") != "cyclo":
        return None
    cy_en, cy_zh = side_en_zh(n, ("", ""))
    if not cy_en:
        return None
    pe, be, pz, bz = _split_pref(*to_benzene(*side_en_zh(s)))
    return _with_n_pref(pe, pz, cy_en, cy_zh, be, bz)


_N_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_N_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}


def _n_alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _N_EN.get(n), _N_ZH.get(n)
    return (en, zh) if en and zh else None


def _n_alkyl_names(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    stem = _sa_alkyl_stem(int(s.get("n") or 0))
    n_stem = _n_alkyl_stem(int(n.get("n") or 0))
    if stem is None or n_stem is None:
        return None
    return f"N-{n_stem[0]}{stem[0]}", f"N-{n_stem[1]}{stem[1]}"


def _n_n_name_en_zh(ns, n1, n2, stem):
    if ns[0] == ns[1]:
        return f"N,N-di{n1[0]}{stem[0]}", f"N,N-二{n1[1]}{stem[1]}"
    return f"N-{n1[0]}-N-{n2[0]}{stem[0]}", f"N-{n1[1]}-N-{n2[1]}{stem[1]}"


def _n_n_dialkyl_names(parent: dict) -> tuple[str, str] | None:
    s, n = parent.get("s_side") or {}, parent.get("n_side") or {}
    stem = _sa_alkyl_stem(int(s.get("n") or 0))
    ns = n.get("ns") or []
    if len(ns) != 2 or stem is None:
        return None
    n1, n2 = _n_alkyl_stem(ns[0]), _n_alkyl_stem(ns[1])
    if n1 is None or n2 is None:
        return None
    return _n_n_name_en_zh(ns, n1, n2, stem)


_SA_MODE = {
    "alkyl": _sa_alkyl_names,
    "aryl": _sa_aryl_stem,
    "n_aryl_alkyl": _n_aryl_alkyl,
    "n_cyclo_aryl": _n_cyclo_aryl,
    "n_alkyl": _n_alkyl_names,
    "n_n_dialkyl": _n_n_dialkyl_names,
}


def sulfonamide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonamide":
        return None
    fn = _SA_MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# sulfonate ester (P-65.3.2)
# ---------------------------------------------------------------------------

_ALKYL_SO_EN = {1: "methanesulfonate", 2: "ethanesulfonate",
                3: "propanesulfonate", 4: "butanesulfonate"}
_ALKYL_SO_ZH = {1: "甲磺酸", 2: "乙磺酸", 3: "丙磺酸", 4: "丁磺酸"}


def _alkyl_o_pair(n: int) -> tuple[str, str] | None:
    """O-side alkyl: EN methyl…; ZH 甲 / 乙 / 丙 / 正丁."""
    en = ester_alkyl_en(n)
    if not en:
        return None
    if n == 4:
        return en, "正丁"
    zh = ester_alkyl_zh(n)
    return (en, zh) if zh else None


def _toluene_retained(en: str, zh: str) -> tuple[str, str]:
    """Gold prefers 对甲苯磺酸 over 4-甲基苯磺酸 for p-tosylates."""
    if en == "4-methylbenzene":
        return en, "对甲苯"
    if en == "2-methylbenzene":
        return en, "邻甲苯"
    if en == "3-methylbenzene":
        return en, "间甲苯"
    return en, zh


def _so3_alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SO_EN.get(n), _ALKYL_SO_ZH.get(n)
    return (en, zh) if en and zh else None


def _so3_aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = _toluene_retained(*to_benzene(*side_en_zh(s)))
    return f"{en}sulfonate", f"{zh}磺酸"


def _so3_s_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") == "aryl":
        return _so3_aryl_stem(parent)
    if s.get("kind") == "alkyl":
        return _so3_alkyl_stem(int(s.get("n") or 0))
    return None


def _so3_join(o: tuple[str, str], stem: tuple[str, str]) -> tuple[str, str]:
    return f"{o[0]} {stem[0]}", f"{stem[1]}{o[1]}酯"


def sulfonate_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonate":
        return None
    o = _alkyl_o_pair(int((parent.get("o_side") or {}).get("n") or 0))
    stem = _so3_s_stem(parent)
    return _so3_join(o, stem) if o and stem else None


# ---------------------------------------------------------------------------
# sulfone (P-65.3.1.2)
# ---------------------------------------------------------------------------

_SONE_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
_SONE_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}
_SULFONE_EN = "sulfone"
_SULFONE_ZH = "砜"


def _sone_alkyl(n: int) -> tuple[str, str] | None:
    en, zh = _SONE_ALKYL_EN.get(n), _SONE_ALKYL_ZH.get(n)
    return (en, zh) if en and zh else None


def _sulfone_symmetric(a1, a2):
    en = f"di{a1[0]} {_SULFONE_EN}"
    zh = f"二{a1[1]}{_SULFONE_ZH}"
    return en, zh


def _sulfone_asymmetric(a1, a2):
    en = f"{a1[0]} {a2[0]} {_SULFONE_EN}"
    zh = f"{a1[1]}{a2[1]}{_SULFONE_ZH}"
    return en, zh


def _sone_names(parent: dict) -> tuple[str, str] | None:
    ns = parent.get("alkyl_ns") or (0, 0)
    n1, n2 = ns[0], ns[1]
    if n1 < 1 or n2 < 1:
        return None
    a1 = _sone_alkyl(n1)
    a2 = _sone_alkyl(n2)
    if a1 is None or a2 is None:
        return None
    return _sulfone_symmetric(a1, a2) if n1 == n2 else _sulfone_asymmetric(a1, a2)


def sulfone_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfone":
        return None
    return _sone_names(parent)


# ---------------------------------------------------------------------------
# sulfonic acid (P-65.3)
# ---------------------------------------------------------------------------

_ALKYL_SOA_EN = {1: "methanesulfonic acid", 2: "ethanesulfonic acid",
                 3: "propanesulfonic acid", 4: "butanesulfonic acid"}
_ALKYL_SOA_ZH = {1: "甲磺酸", 2: "乙磺酸", 3: "丙磺酸", 4: "丁磺酸"}


def _soh_alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SOA_EN.get(n), _ALKYL_SOA_ZH.get(n)
    return (en, zh) if en and zh else None


def _soh_alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _soh_alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _soh_aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = to_benzene(*side_en_zh(s))
    return f"{en}sulfonic acid", f"{zh}磺酸"


_SOH_MODE = {
    "alkyl": _soh_alkyl_names,
    "aryl": _soh_aryl_stem,
}


def sulfonic_acid_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonic_acid":
        return None
    fn = _SOH_MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# sulfonyl chloride (P-65.3)
# ---------------------------------------------------------------------------

_ALKYL_SC_EN = {1: "methanesulfonyl chloride", 2: "ethanesulfonyl chloride",
                3: "propanesulfonyl chloride", 4: "butanesulfonyl chloride"}
_ALKYL_SC_ZH = {1: "甲磺酰氯", 2: "乙磺酰氯", 3: "丙磺酰氯", 4: "丁磺酰氯"}


def _scl_alkyl_stem(n: int) -> tuple[str, str] | None:
    en, zh = _ALKYL_SC_EN.get(n), _ALKYL_SC_ZH.get(n)
    return (en, zh) if en and zh else None


def _scl_alkyl_names(parent: dict) -> tuple[str, str] | None:
    return _scl_alkyl_stem(int((parent.get("s_side") or {}).get("n") or 0))


def _scl_aryl_stem(parent: dict) -> tuple[str, str] | None:
    s = parent.get("s_side") or {}
    if s.get("kind") != "aryl":
        return None
    en, zh = to_benzene(*side_en_zh(s))
    return f"{en}sulfonyl chloride", f"{zh}磺酰氯"


_SCL_MODE = {
    "alkyl": _scl_alkyl_names,
    "aryl": _scl_aryl_stem,
}


def sulfonyl_chloride_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "sulfonyl_chloride":
        return None
    fn = _SCL_MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# sulfoxide (P-63.3)
# ---------------------------------------------------------------------------

# Symmetric: dimethyl sulfoxide / 二甲基亚砜 (gold uses 基)
SULFOXIDE_SYM_EN = {
    1: "dimethyl sulfoxide", 2: "diethyl sulfoxide",
    3: "dipropyl sulfoxide", 4: "dibutyl sulfoxide",
}
SULFOXIDE_SYM_ZH = {
    1: "二甲基亚砜", 2: "二乙基亚砜", 3: "二丙基亚砜", 4: "二丁基亚砜",
}


def _sym_names(n: int) -> tuple[str, str] | None:
    en, zh = SULFOXIDE_SYM_EN.get(n), SULFOXIDE_SYM_ZH.get(n)
    return (en, zh) if en and zh else None


def _asym_names(n1: int, n2: int) -> tuple[str, str] | None:
    en1, en2 = SULFIDE_ALKYL_EN.get(n1), SULFIDE_ALKYL_EN.get(n2)
    zh1, zh2 = SULFIDE_ALKYL_ZH.get(n1), SULFIDE_ALKYL_ZH.get(n2)
    if not en1 or not en2 or not zh1 or not zh2:
        return None
    a, b = sorted([(en1, zh1), (en2, zh2)], key=lambda x: x[0])
    return f"{a[0]} {b[0]} sulfoxide", f"{a[1]}{b[1]}亚砜"


def sulfoxide_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    ns = parent.get("alkyl_ns")
    if not ns or len(ns) != 2:
        return None
    n1, n2 = int(ns[0]), int(ns[1])
    if n1 == n2:
        return _sym_names(n1)
    return _asym_names(n1, n2)
