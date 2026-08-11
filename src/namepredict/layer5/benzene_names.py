"""Retained / multi-substituted benzene parent and prefix helpers (P-22.1.3)."""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer3.substituent_extractor import alkyl_alpha_key




def _xylene_locants(numbered: dict) -> str:
    locs = sorted(int(s["locant"]) for s in numbered.get("substituents") or [] if "locant" in s)
    return ",".join(str(x) for x in locs)


def _methoxy_subs(subs: list) -> list:
    return [s for s in subs if s.get("kind") == "alkoxy" and s.get("n_carbons") == 1]


def _is_poly_anisole(numbered: dict) -> bool:
    return None

def _anisole_reloc(loc: int, meo: int, rev: bool) -> int:
    d = (meo - loc) % 6 if rev else (loc - meo) % 6
    return d + 1


def _anisole_dir_key(others: list, meo: int, rev: bool) -> tuple:
    """Locants in alphabetical citation order (not sorted-set)."""
    ordered = sorted(
        (s for s in others if "locant" in s),
        key=lambda s: alkyl_alpha_key(s.get("en") or ""),
    )
    return tuple(_anisole_reloc(int(s["locant"]), meo, rev) for s in ordered)


def _anisole_pick_rev(others: list, meo: int) -> bool:
    return _anisole_dir_key(others, meo, True) < _anisole_dir_key(others, meo, False)


def _anisole_renum(others: list, meo_loc: int, rev: bool) -> list:
    return [{**s, "locant": _anisole_reloc(int(s["locant"]), meo_loc, rev)} for s in others]




def _xylene_prefix(numbered: dict) -> tuple[str, str]:
    locs = _xylene_locants(numbered)
    return f"{locs}-", f"{locs}-二甲基"


def benzene_prefix(numbered: dict, build_prefix) -> tuple[str, str]:
    en_pre, zh_pre = build_prefix(numbered.get("substituents") or [], 6, "benzene")
    return (en_pre, None) if _is_poly_anisole(numbered) else (en_pre, zh_pre)


def _stereo_lead(parent: str) -> tuple[str, str]:
    """Split leading stereo '(E)-' / '(9Z,11E)-' from parent stem."""
    if not parent.startswith("("):
        return "", parent
    close = parent.find(")-")
    if close < 0:
        return "", parent
    return parent[: close + 2], parent[close + 2 :]


def join_parent_name(prefix: str, parent: str) -> str:
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if stem[:1].isdigit() or stem.startswith("1H-") else f"{prefix}{stem}"
    return f"{stereo}{body}"


def _ester_alkoxy_from(numbered) -> tuple[str, str]:
    """O-side alkyl: linear 走 parent.alkoxy_n 保留名表; 特殊基团/复杂回落 o_side 取代基.

    linear (如 hexadecan-16-yl 的 o_side 命名带位次) 必须走 ester_alkyl_en/zh 的
    "hexadecyl"/"十六", 否则长链 O 侧会带 -16-yl 位次 (ester/benzoate 共用此路径).
    """
    if not numbered:
        return "", ""
    parent = numbered.get("parent") or {}
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    if len(o) <= 1 and parent.get("alkoxy_n") is not None:
        from namepredict.layer5.stems import ester_alkyl_en, ester_alkyl_zh
        en, zh = ester_alkyl_en(parent["alkoxy_n"]), ester_alkyl_zh(parent["alkoxy_n"])
        if en and zh:
            return en, zh
    if not o:
        return "", ""
    en0, zh0 = o[0].get("en") or "", (o[0].get("zh") or "").rstrip("基")
    if len(o) == 1:
        return en0, zh0
    # diester: identical arms aggregate (dimethyl / 二甲)
    if all(s.get("en") == en0 for s in o):
        me, mz = MULT_EN.get(len(o), ""), MULT_ZH.get(len(o), "")
        return f"{me}{en0}", f"{mz}{zh0}"
    return en0, zh0


def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    """Acid stem + O-side alkyl: methyl butanoate / 丁酸甲酯; stereo after alkyl."""
    en, zh = names
    alk_en, alk_zh = _ester_alkoxy_from(numbered)
    st, body = _stereo_lead(en)
    mid = f"{pre_en}{body}" if pre_en else body
    en = f"{alk_en} {st}{mid}" if alk_en else f"{st}{mid}"
    stz, bodyz = _stereo_lead(zh)
    midz = f"{stz}{pre_zh}{bodyz}" if pre_zh else f"{stz}{bodyz}"
    zh = f"{midz}{alk_zh}酯" if alk_zh else midz
    return en, zh


def join_benzoate_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    """O-side alkyl + benzoate acid stem: ethyl 4-chlorobenzoate / 4-氯苯甲酸乙酯.

    与 ester 同构: 母体只含酸部分, 烷氧基从 o_side 取代基取; zh 恒拼"酯"
    (无烷氧基时保留"苯甲酸酯", 如复杂 O-烷基场景).
    """
    en, zh = names
    alk_en, alk_zh = _ester_alkoxy_from(numbered)
    body = f"{pre_en}{en}" if pre_en else en
    en = f"{alk_en} {body}" if alk_en else body
    midz = f"{pre_zh}{zh}" if pre_zh else zh
    return en, f"{midz}{alk_zh}酯"


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str]:
    if kind in ("ester", "diester"):
        return join_ester_name(pre[0], pre[1], names, numbered)
    if kind == "benzoate":
        return join_benzoate_name(pre[0], pre[1], names, numbered)
    en = join_parent_name(pre[0], names[0])
    zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    return en, zh


def zh_1h_parent(en_parent: str, zh_parent: str, prefix: str) -> str:
    """Prefix Chinese retained 1H-parents with 1H- when ring is substituted."""
    if not prefix or not en_parent.startswith("1H-") or zh_parent.startswith("1H-"):
        return zh_parent
    return f"1H-{zh_parent}"
