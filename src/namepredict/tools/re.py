"""命名文本辅助函数。"""
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
    """去掉一个前导位次集，含 1H- 指示氢（P-14.5）。"""
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
    """去掉词干开头的外层括号，含剥位次后残留的括号。"""
    while stem.startswith("("):
        stem = stem[1:]
        if stem.endswith(")"):
            stem = stem[:-1]
    return stem


_STEREO_LEAD_RE = re.compile(
    r"^\((?:\d*[a-zA-Z]*[EeZzRrSs])(?:,(?:\d*[a-zA-Z]*[EeZzRrSs]))*\)-"
)


def _strip_lead_stereo(stem: str) -> str:
    """剥除词干最前的立体描述符组及其连字符（P-14.5）。"""
    m = _STEREO_LEAD_RE.match(stem)
    return stem[m.end():] if m else stem


def alkyl_alpha_key(stem: str) -> str:
    """字母数字序键：忽略斜体前缀/括号/位次/立体组（P-14.5），循环剥到稳定。"""
    s = _strip_n_prefix(_strip_ital_prefix(stem))
    while True:
        s2 = _strip_lead_locant(_strip_outer_parens(_strip_lead_stereo(s[1:] if s.startswith("[") else s)))  # P-14.5
        if s2 == s:
            return s
        s = s2


# ── 名称文本规范化（判分口径） ──────────────────────
_WS = re.compile(r"\s+")


def normalize_en(name: str) -> str:
    """规范化英文名：小写、去重空白并统一连字符/逗号/括号。"""
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = s.replace("[", "(").replace("]", ")")
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    """规范化中文名：去首尾空白并统一括号种类（方/圆等价，仅括注外观不同不判分）。"""
    s = (name or "").strip()
    s = s.replace("[", "(").replace("]", ")")
    return s


def nospace(name: str) -> str:
    """删去全部空白字符。"""
    return "".join((name or "").split())
