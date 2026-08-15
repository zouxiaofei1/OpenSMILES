"""保留名 / 多取代苯母体与前缀辅助（P-22.1.3）。"""
from __future__ import annotations

from namepredict.constants import MULT_EN, MULT_ZH
from namepredict.layer5.stereo import _split_stereo_lead as _stereo_lead

def join_parent_name(prefix: str, parent: str) -> str:
    if not prefix:
        return parent
    stereo, stem = _stereo_lead(parent)
    body = f"{prefix}-{stem}" if stem[:1].isdigit() or stem.startswith("1H-") else f"{prefix}{stem}"
    return f"{stereo}{body}"

def join_ester_name(pre_en: str, pre_zh: str, names: tuple[str, str], numbered=None) -> tuple[str, str]:
    en, zh = names
    o = [s for s in (numbered.get("substituents") or []) if s.get("o_side")]
    alk_en, alk_zh = o[0].get("en") or "", (o[0].get("zh") or "").rstrip("基")
    st, body = _stereo_lead(en)
    mid = f"{pre_en}{body}" if pre_en else body
    en = f"{alk_en} {st}{mid}" if alk_en else f"{st}{mid}"
    stz, bodyz = _stereo_lead(zh)
    midz = f"{stz}{pre_zh}{bodyz}" if pre_zh else f"{stz}{bodyz}"
    zh = f"{midz}{alk_zh}酯"
    return en, zh


def join_kind_name(
    kind: str | None, pre: tuple[str, str], names: tuple[str, str],
    numbered=None,
) -> tuple[str, str]:
    if kind in ("ester"):
        return join_ester_name(pre[0], pre[1], names, numbered)
    en = join_parent_name(pre[0], names[0])
    zh = join_parent_name(pre[1], zh_1h_parent(names[0], names[1], pre[1]))
    return en, zh


def zh_1h_parent(en_parent: str, zh_parent: str, prefix: str) -> str:
    """环被取代时，中文保留名 1H- 母体补加 1H- 前缀。"""
    if not prefix or not en_parent.startswith("1H-") or zh_parent.startswith("1H-"):
        return zh_parent
    return f"1H-{zh_parent}"
