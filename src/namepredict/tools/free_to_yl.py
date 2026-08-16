"""free 母体名 → P-29 -yl 取代基形式的双语转换（含官能团前缀化）。"""

from __future__ import annotations

import re

# free_en → (yl_en, yl_zh)（当位次为常规/省略时）
_RETAINED: dict[str, tuple[str, str]] = {
    "benzene": ("phenyl", "苯基"),
}


def _drop_terminal_e(en: str) -> str:
    """去掉英文名末尾的 e（用于拼接 -yl 词干）。"""
    return en[:-1] if en.endswith("e") else en


def _hetero_locant_prefix(en: str) -> str:
    """取 free EN 的前导杂原子位次：'1,3-thiazole' → '1,3-'。"""
    i = 0
    n = len(en)
    while i < n and en[i].isdigit():
        i += 1
        while i < n and en[i] == ",":
            i += 1
            while i < n and en[i].isdigit():
                i += 1
    return en[: i + 1] if i and i < n and en[i] == "-" else ""


# ── 官能团后缀 → 前缀转换 ──

_MULT_OL_SUF = re.compile(r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca)ol\b")


def _alkoxy_en(en: str) -> str | None:
    """醇 → 烷氧基：methanol→methoxy, propan-1-ol→propoxy, propan-2-ol→propan-2-yloxy。"""
    if not en.endswith("ol") or _MULT_OL_SUF.search(en):
        return None  # 多醇（diol…decaol）不能去氢成烷氧基
    if en in ("methanol", "ethanol"):
        return en[:-4] + "oxy"
    m = re.match(r"^(.+)an-1-ol$", en)
    if m:
        return m.group(1) + "oxy"
    m = re.match(r"^(.+)an-(\d+)-ol$", en)
    if m:
        return f"{m.group(1)}an-{m.group(2)}-yloxy"
    # heptan-1-ol 等（C7+，带完整 "ane" 词干）
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
    """硫醇 → 硫基：methanethiol→methylsulfanyl, ethanethiol→ethylsulfanyl, propane-2-thiol→propane-2-ylsulfanyl。"""
    if not en.endswith("thiol"):
        return None
    # 例：propane-2-thiol → propane-2-ylsulfanyl
    m = re.match(r"^(.+)ane-(\d+)-thiol$", en)
    if m:
        return f"{m.group(1)}ane-{m.group(2)}-ylsulfanyl"
    # 例：methanethiol → methylsulfanyl
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
    """1° 胺 → 氨基：methanamine→methylamino, ethanamine→ethylamino, propan-2-amine→propan-2-ylamino。"""
    if not en.endswith("amine") or "N-" in en:
        return None
    # 例：propan-2-amine → propan-2-ylamino
    m = re.match(r"^(.+)ane-(\d+)-amine$", en)
    if m:
        return f"{m.group(1)}ane-{m.group(2)}-ylamino"
    # 例：methanamine → methylamino
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


# 单核母体氢化物（P-15.4.1 表 2.1）→ 去氢取代基名（表 1.5 'a' 前缀体系）。
# 组装名 "ethyl-oxidane" → "ethyloxy" / "乙基-氧化烷" → "乙氧基"。
_MONONUCLEAR = (("oxidane", "氧化烷", "oxy", "氧基"),
                ("azane", "氮烷", "amino", "氨基"),
                ("sulfane", "硫烷", "sulfanyl", "硫基"))


def _mononuclear_en(en: str) -> str | None:
    """单核氢化物 free 名 → 去氢取代基名：ethyl-oxidane → ethyloxy。"""
    for en_suf, _, yl, _ in _MONONUCLEAR:
        if en.endswith("-" + en_suf):
            return en[: -len(en_suf) - 1] + yl
    return None


def _mononuclear_zh(zh: str) -> str | None:
    """中文组装名去氢：乙基-氧化烷 → 乙氧基、乙基-氮烷 → 乙氨基。"""
    for _, zh_suf, _, zy in _MONONUCLEAR:
        if zh.endswith("-" + zh_suf):
            base = zh[: -len(zh_suf) - 1]
            return base[: -1] + zy if base.endswith("基") else base + zy
    return None


# ── yl 形式构建器 ──

def _try_fg_prefix(en: str, zh: str) -> tuple[str, str] | None:
    """尝试将官能团后缀名转换为取代基前缀形式。"""
    for en_fn, zh_fn in [
        (_alkoxy_en, _alkoxy_zh),
        (_sulfanyl_en, _sulfanyl_zh),
        (_amino_en, _amino_zh),
        (_mononuclear_en, _mononuclear_zh),
    ]:
        en_out = en_fn(en)
        zh_out = zh_fn(zh) if en_out else None
        if en_out and zh_out:
            return en_out, zh_out
    return None


def _yl_en(en: str, k: int) -> str:
    """生成英文 -yl 形式；保留名用映射，否则按位次拼接。"""
    hit = _RETAINED.get(en)
    return hit[0] if hit is not None else f"{_drop_terminal_e(en)}-{k}-yl"


def _mirror_azole_locants_zh(zh: str, en: str) -> str:
    """当 EN 含 1,3- 而 free ZH 省略时，向 ZH 唑类词干插入 1,3-。"""
    pre = _hetero_locant_prefix(en)
    if pre and not zh.startswith(pre):
        return pre + zh
    if "1,3-thiazole" in en and "1,3-噻唑" not in zh and zh.endswith("噻唑"):
        return zh[: -len("噻唑")] + "-1,3-噻唑" if zh != "噻唑" else "1,3-噻唑"
    return zh


def _yl_zh(zh: str, k: int, en: str) -> str:
    """生成中文 -基形式；保留名用映射，否则按位次拼接。"""
    hit = _RETAINED.get(en)
    if hit is not None:
        return hit[1]
    return f"{_mirror_azole_locants_zh(zh, en)}-{k}-基"


def free_to_yl(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool]:
    """在键合位次处将 free 母体名转换为 P-29 -yl 双语形式。

    P-63.2.2 醇→烷氧基, P-63.2.1 硫醇→烷基硫基, P-62.2 1° 胺→烷基氨基。
    """
    fg = _try_fg_prefix(en, zh)
    if fg is not None:
        # P-29.3.6：复合前缀（methylamino=CH3-NH-，非普通 amino）需括号与两个独立取代基区分。
        need_paren = fg[0].endswith("amino") and fg[0] != "amino"
        return fg[0], fg[1], need_paren
    return _yl_en(en, attach_locant), _yl_zh(zh, attach_locant, en), paren
