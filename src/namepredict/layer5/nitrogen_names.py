"""L5 names for nitrogen functional parents.

Covers urea / guanidine / hydrazine / isocyanate (P-66, P-68). Sides are
precomputed into parent dicts at L2 pack time; L5 only reads side['en']/
side['zh'] (and existing kind/mode fields).
"""
from __future__ import annotations

from namepredict.layer5.aryl_helpers import side_en_zh, tolyl_swap, wrap_aryl
from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh

# ---------------------------------------------------------------------------
# urea (P-66.1.6.1.1)
# ---------------------------------------------------------------------------

_UREA_BARE = frozenset({"phenyl", "p-tolyl", "o-tolyl", "m-tolyl"})


def _urea_aryl_label(side: dict) -> tuple[str, str]:
    return tolyl_swap(*side_en_zh(side))


def _urea_unsub_names() -> tuple[str, str]:
    return "urea", "脲"


def _urea_mono_aryl_names(parent: dict) -> tuple[str, str]:
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = wrap_aryl(*_urea_aryl_label(ar), bare=_UREA_BARE)
    return f"{en}urea", f"{zh}脲"


def _urea_me2_aryl_names(parent: dict) -> tuple[str, str]:
    """Benchmark style: 3-(4-chlorophenyl)-1,1-dimethylurea."""
    n1, n3 = parent.get("n1") or {}, parent.get("n3") or {}
    ar = n3 if n3.get("kind") == "aryl" else n1
    en, zh = _urea_aryl_label(ar)
    return (
        f"3-({en})-1,1-dimethylurea",
        f"3-({zh})-1,1-二甲基脲",
    )


def _urea_pattern(parent: dict) -> str:
    kinds = {(parent.get("n1") or {}).get("kind"), (parent.get("n3") or {}).get("kind")}
    if kinds == {"h"}:
        return "unsub"
    if kinds == {"aryl", "h"}:
        return "mono_aryl"
    if kinds == {"dialkyl", "aryl"}:
        return "me2_aryl"
    return ""


_UREA_PAT = {
    "unsub": lambda p: _urea_unsub_names(),
    "mono_aryl": _urea_mono_aryl_names,
    "me2_aryl": _urea_me2_aryl_names,
}


def urea_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "urea":
        return None
    fn = _UREA_PAT.get(_urea_pattern(parent))
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# guanidine (P-66.4.1.2.1)
# ---------------------------------------------------------------------------

def _gua_unsub_names() -> tuple[str, str]:
    return "guanidine", "胍"


def _gua_mono_aryl_names(parent: dict) -> tuple[str, str]:
    en, zh = wrap_aryl(*side_en_zh(parent.get("sub") or {}))
    return f"1-{en}guanidine", f"1-{zh}胍"


def _gua_arylsulfonyl_names(parent: dict) -> tuple[str, str]:
    en, zh = side_en_zh(parent.get("sub") or {})
    return (
        f"1-({en}sulfonyl)guanidine",
        f"1-({zh}磺酰基)胍",
    )


_GUA_MODE = {
    "unsub": lambda p: _gua_unsub_names(),
    "aryl": _gua_mono_aryl_names,
    "arylsulfonyl": _gua_arylsulfonyl_names,
}


def guanidine_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "guanidine":
        return None
    fn = _GUA_MODE.get(parent.get("mode") or "")
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# hydrazine (P-68.3.1.2)
# ---------------------------------------------------------------------------

def _hyd_unsub_names() -> tuple[str, str]:
    return "hydrazine", "肼"


def _hyd_phenyl_names() -> tuple[str, str]:
    return "phenylhydrazine", "苯肼"


def _hyd_alkyl_label(n: int) -> tuple[str, str] | None:
    en, zh = ester_alkyl_en(n), ester_alkyl_zh(n)
    return (en, zh) if en and zh else None


def _hyd_sym_dialkyl_names(ns: list[int]) -> tuple[str, str] | None:
    if len(ns) != 2 or ns[0] != ns[1]:
        return None
    lab = _hyd_alkyl_label(ns[0])
    if lab is None:
        return None
    en_a, zh_a = lab
    return f"1,1-di{en_a}hydrazine", f"1,1-二{zh_a}肼"


def _hyd_pattern(parent: dict) -> str:
    n1 = (parent.get("n1") or {}).get("kind")
    n2 = (parent.get("n2") or {}).get("kind")
    kinds = {n1, n2}
    if kinds == {"h"}:
        return "unsub"
    if kinds == {"phenyl", "h"}:
        return "phenyl"
    if kinds == {"dialkyl", "h"}:
        return "dialkyl"
    return ""


def _hyd_dialkyl_names(parent: dict) -> tuple[str, str] | None:
    n1 = parent.get("n1") or {}
    n2 = parent.get("n2") or {}
    side = n1 if n1.get("kind") == "dialkyl" else n2
    return _hyd_sym_dialkyl_names(list(side.get("ns") or []))


_HYD_PAT = {
    "unsub": lambda _p: _hyd_unsub_names(),
    "phenyl": lambda _p: _hyd_phenyl_names(),
    "dialkyl": _hyd_dialkyl_names,
}


def hydrazine_names(numbered: dict) -> tuple[str, str] | None:
    parent = numbered.get("parent") or {}
    if parent.get("kind") != "hydrazine":
        return None
    fn = _HYD_PAT.get(_hyd_pattern(parent))
    return fn(parent) if fn else None


# ---------------------------------------------------------------------------
# isocyanate / isothiocyanate (P-61.9)
# ---------------------------------------------------------------------------

def _iso_alkyl_pair(n: int) -> tuple[str, str] | None:
    en, zh = ester_alkyl_en(n), ester_alkyl_zh(n)
    return (en, zh) if en and zh else None


def _iso_pair(n: int, en_tail: str, zh_head: str) -> tuple[str, str] | None:
    alk = _iso_alkyl_pair(n)
    if alk is None:
        return None
    en_a, zh_a = alk
    return f"{en_a} {en_tail}", f"{zh_head}{zh_a}酯"


def isocyanate_names(n: int) -> tuple[str, str] | None:
    return _iso_pair(n, "isocyanate", "异氰酸")


def isothiocyanate_names(n: int) -> tuple[str, str] | None:
    return _iso_pair(n, "isothiocyanate", "异硫氰酸")


def iso_kind_names(kind: str, n: int) -> tuple[str, str] | None:
    if kind == "isocyanate":
        return isocyanate_names(n)
    if kind == "isothiocyanate":
        return isothiocyanate_names(n)
    return None
