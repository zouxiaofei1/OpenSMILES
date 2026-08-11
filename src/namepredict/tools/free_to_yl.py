
from __future__ import annotations

import re

# free_en → (yl_en, yl_zh) when locant is conventional / omitted
_RETAINED: dict[str, tuple[str, str]] = {
    "benzene": ("phenyl", "苯基"),
}


def _drop_terminal_e(en: str) -> str:
    return en[:-1] if en.endswith("e") else en


def _hetero_locant_prefix(en: str) -> str:
    """Leading heteroatom locants from free EN: '1,3-thiazole' → '1,3-'."""
    i = 0
    n = len(en)
    while i < n and en[i].isdigit():
        i += 1
        while i < n and en[i] == ",":
            i += 1
            while i < n and en[i].isdigit():
                i += 1
    return en[: i + 1] if i and i < n and en[i] == "-" else ""


# ── FG suffix → prefix conversion ──

def _alkoxy_en(en: str) -> str | None:
    """alcohol → alkoxy: methanol→methoxy, propan-1-ol→propoxy, propan-2-ol→propan-2-yloxy."""
    if not en.endswith("ol") or "diol" in en or "triol" in en:
        return None
    if en in ("methanol", "ethanol"):
        return en[:-4] + "oxy"
    m = re.match(r"^(.+)an-1-ol$", en)
    if m:
        return m.group(1) + "oxy"
    m = re.match(r"^(.+)an-(\d+)-ol$", en)
    if m:
        return f"{m.group(1)}an-{m.group(2)}-yloxy"
    # heptan-1-ol etc. (C7+ with full "ane" stem)
    m = re.match(r"^(.+)ane-1-ol$", en)
    if m:
        return m.group(1) + "oxy"
    return None


def _alkoxy_zh(zh: str) -> str | None:
    """alcohol → alkoxy ZH: 甲醇→甲氧基, 丙-1-醇→丙氧基, 丙-2-醇→丙-2-基氧基."""
    if not zh.endswith("醇") or "二酚" in zh:
        return None
    if zh in ("甲醇", "乙醇"):
        return zh[:-1] + "氧基"
    m = re.match(r"^(.+)-1-醇$", zh)
    if m:
        return m.group(1) + "氧基"
    m = re.match(r"^(.+)-(\d+)-醇$", zh)
    if m:
        return f"{m.group(1)}-{m.group(2)}-基氧基"
    return None


def _sulfanyl_en(en: str) -> str | None:
    """thiol → sulfanyl: methanethiol→methylsulfanyl, ethanethiol→ethylsulfanyl, propane-2-thiol→propane-2-ylsulfanyl."""
    if not en.endswith("thiol"):
        return None
    # propane-2-thiol → propane-2-ylsulfanyl
    m = re.match(r"^(.+)ane-(\d+)-thiol$", en)
    if m:
        return f"{m.group(1)}ane-{m.group(2)}-ylsulfanyl"
    # methanethiol → methylsulfanyl
    m = re.match(r"^(.+)anethiol$", en)
    if m:
        return m.group(1) + "ylsulfanyl"
    return None


def _sulfanyl_zh(zh: str) -> str | None:
    """thiol → sulfanyl ZH: 甲硫醇→甲硫基, 乙硫醇→乙硫基, 丙-2-硫醇→丙-2-基硫基."""
    if zh in ("甲硫醇", "乙硫醇"):
        return zh[:-1] + "基"
    m = re.match(r"^(.+)-(\d+)-硫醇$", zh)
    if m:
        return f"{m.group(1)}-{m.group(2)}-基硫基"
    return None


def _amino_en(en: str) -> str | None:
    """1° amine → amino: methanamine→methylamino, ethanamine→ethylamino, propan-2-amine→propan-2-ylamino."""
    if not en.endswith("amine") or "N-" in en:
        return None
    # propan-2-amine → propan-2-ylamino
    m = re.match(r"^(.+)ane-(\d+)-amine$", en)
    if m:
        return f"{m.group(1)}ane-{m.group(2)}-ylamino"
    # methanamine → methylamino
    m = re.match(r"^(.+)anamine$", en)
    if m:
        return m.group(1) + "ylamino"
    return None


def _amino_zh(zh: str) -> str | None:
    """1° amine → amino ZH: 甲胺→甲氨基, 乙胺→乙氨基, 丙-2-胺→丙-2-基氨基."""
    if zh in ("甲胺", "乙胺"):
        return zh[:-1] + "氨基"
    m = re.match(r"^(.+)-(\d+)-胺$", zh)
    if m:
        return f"{m.group(1)}-{m.group(2)}-基氨基"
    return None


# ── yl form builders ──

def _try_fg_prefix(en: str, zh: str) -> tuple[str, str] | None:
    """Try to convert FG suffix name to substituent prefix form."""
    for en_fn, zh_fn in [
        (_alkoxy_en, _alkoxy_zh),
        (_sulfanyl_en, _sulfanyl_zh),
        (_amino_en, _amino_zh),
    ]:
        en_out = en_fn(en)
        zh_out = zh_fn(zh) if en_out else None
        if en_out and zh_out:
            return en_out, zh_out
    return None


def _yl_en(en: str, k: int) -> str:
    hit = _RETAINED.get(en)
    return hit[0] if hit is not None else f"{_drop_terminal_e(en)}-{k}-yl"


def _mirror_azole_locants_zh(zh: str, en: str) -> str:
    """Insert 1,3- into ZH azole stems when EN has 1,3- but free ZH omits it."""
    pre = _hetero_locant_prefix(en)
    if pre and not zh.startswith(pre):
        return pre + zh
    if "1,3-thiazole" in en and "1,3-噻唑" not in zh and zh.endswith("噻唑"):
        return zh[: -len("噻唑")] + "-1,3-噻唑" if zh != "噻唑" else "1,3-噻唑"
    return zh


def _yl_zh(zh: str, k: int, en: str) -> str:
    hit = _RETAINED.get(en)
    if hit is not None:
        return hit[1]
    return f"{_mirror_azole_locants_zh(zh, en)}-{k}-基"


def free_to_yl(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool]:
    """Convert free parent name to P-29 -yl dual form at attach locant.

    P-63.2.2 alcohol→alkoxy, P-63.2.1 thiol→alkylsulfanyl, P-62.2 1° amine→alkylamino.
    """
    fg = _try_fg_prefix(en, zh)
    if fg is not None:
        # P-29.3.6: compound prefixes (methylamino=CH3-NH-, not plain amino)
        # need parentheses to distinguish from two separate substituents.
        need_paren = fg[0].endswith("amino") and fg[0] != "amino"
        return fg[0], fg[1], need_paren
    return _yl_en(en, attach_locant), _yl_zh(zh, attach_locant, en), paren
