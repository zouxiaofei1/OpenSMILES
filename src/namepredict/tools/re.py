"""命名文本剥除工具：P-14.5 字母数字序键 `alkyl_alpha_key` 及其前缀剥除辅助函数。"""
from __future__ import annotations

import re


def _strip_ital_prefix(stem: str) -> str:
    """去掉 sec-/tert- 前缀。"""
    if stem.startswith("tert-") or stem.startswith("sec-"):
        return stem[stem.index("-") + 1 :]
    return stem


def _strip_n_prefix(stem: str) -> str:
    """去掉 N- 或 N, 前缀。"""
    if stem.startswith("N,"):
        return stem.split("-")[-1] if "-" in stem else stem
    return stem[2:] if stem.startswith("N-") else stem


def _strip_lead_locant(stem: str) -> str:
    """去掉一个前导位次集：'4-'、'1,3-'、'1,1,1-'、'1H-'（P-14.5）。"""
    i = 0
    n = len(stem)
    while i < n and stem[i].isdigit():
        i += 1
        while i < n and stem[i] == ",":
            i += 1
            while i < n and stem[i].isdigit():
                i += 1
    if i and i + 1 < n and stem[i] == "H" and stem[i + 1] == "-":  # 指示氢前缀：1H-、2H-、3H-（P-14.5 / P-65.3.2.5）
        i += 1
    return stem[i + 1 :] if i and i < n and stem[i] == "-" else stem


def _strip_outer_parens(stem: str) -> str:
    """去掉外层括号；词干以开括号开头即剥（含剥 locant 后残留的前缀括号，如 '(phenyloxy)methyl'）。"""
    while stem.startswith("("):
        stem = stem[1:]
        if stem.endswith(")"):
            stem = stem[:-1]
    return stem


_STEREO_LEAD_RE = re.compile(
    r"^\((?:\d*[a-zA-Z]*[EeZzRrSs])(?:,(?:\d*[a-zA-Z]*[EeZzRrSs]))*\)-"
)


def _strip_lead_stereo(stem: str) -> str:
    """剥除词干最前的立体描述符组 '(2S,3R)-' / '(E)-' 及其后的连字符（P-14.5：字母数字序不含立体描述符）。"""
    m = _STEREO_LEAD_RE.match(stem)
    return stem[m.end():] if m else stem


def _strip_lead_bracket(stem: str) -> str:
    """剥一个前导 '['：复合前缀整体被方括号包时，开括号本身不参与字母序，须剥到其后的实质词干（P-14.5）。"""
    return stem[1:] if stem.startswith("[") else stem


def alkyl_alpha_key(stem: str) -> str:
    """字母数字序键：忽略 sec-/tert-/N-/括号/前导位次/前导立体组（P-14.5）；交替剥括号、立体组与前导位次直到稳定，避免残留开括号/位次数字按 ASCII 错误排最前。"""
    s = _strip_n_prefix(_strip_ital_prefix(stem))
    while True:
        s2 = _strip_lead_locant(_strip_outer_parens(_strip_lead_stereo(_strip_lead_bracket(s))))
        if s2 == s:
            return s
        s = s2
