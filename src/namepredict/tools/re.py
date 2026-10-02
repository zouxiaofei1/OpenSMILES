"""命名文本辅助函数。

alpha_order_key 系列供管线内 P-14.5 排序；normalize_en/normalize_zh/nospace 是
测试与 benchmark 的判分口径，管线本身不调用它们。
"""
from __future__ import annotations

import re





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

# 立体描述符剥除（P-14.5 排序用）：允许无位次的 (E)-、(R)-，也允许 1aR 型。与
# layer5.assembler_prefixes._STEREO_LEAD_ENCLOSE_RE 不同——那份要求每个 token 带位次，
# 用于判"是否须整体围栏"，两者语义不同，勿互换。
_STEREO_LEAD_STRIP_RE = re.compile(
    r"^\((?:\d*[a-zA-Z]*[EeZzRrSs])(?:,(?:\d*[a-zA-Z]*[EeZzRrSs]))*\)-"
)

# 取代基名额外的位次段（复合环名/复合臂名判据）：段首或 -([ 之后起，排除母体词干内烯/炔位次。
# layer5.assembler 与 tools.anchored_table 共用这一份，勿各自内联。
SUB_LOCANT_RE = re.compile(r"(?:^|[-(\[])\d+(?:,\d+)*[a-z]?-(?!(?:en|yn|an|in))")


def _strip_lead_stereo(stem: str) -> str:
    """剥除词干最前的立体描述符组及其连字符（P-14.5）。"""
    m = _STEREO_LEAD_STRIP_RE.match(stem)
    return stem[m.end():] if m else stem


def alkyl_alpha_key(stem: str) -> str:
    """字母数字序键：忽略斜体前缀/括号/位次/立体组（P-14.5），循环剥到稳定。"""
    s = stem
    while True:
        s2 = _strip_lead_locant(_strip_lead_stereo(s[1:] if s.startswith("[") else s))  # P-14.5
        if s2 == s:
            return s
        s = s2


_MID_ITALIC_RE = re.compile(r"(?<![A-Za-z])(?:tert|sec)-")  # 词中斜体前缀 tert-/sec-
_LOWER_LETTER_RE = re.compile(r"[a-z]")  # 非斜体罗马字母（斜体 R/S/E/Z/H/N 为大写，不参与首轮比较）
_DIGIT_RUN_RE = re.compile(r"\d+")


def _nonitalic_letters(stem: str) -> str:
    """取词干中的非斜体罗马字母：位次/连字符/括号/立体描述符一律不参与比较。"""
    return "".join(_LOWER_LETTER_RE.findall(_MID_ITALIC_RE.sub("", stem)))


def _lead_locants(stem: str) -> tuple[int, ...]:
    """首个非斜体罗马字母之前的位次集（P-14.5.4 平局判据），无位次最优先。"""
    s = _MID_ITALIC_RE.sub("", stem)
    m = _LOWER_LETTER_RE.search(s)
    head = s[: m.start()] if m else s
    return tuple(sorted(int(d) for d in _DIGIT_RUN_RE.findall(head)))


def alpha_order_key(stem: str) -> tuple:
    """P-14.5 字母数字序键：先比字母序列（忽略位次/连字符/斜体），再比首字母前位次。"""
    s = alkyl_alpha_key(stem)
    return (_nonitalic_letters(s), _lead_locants(stem))


# ── 名称文本规范化（判分口径） ──────────────────────
_WS = re.compile(r"\s+")


def _unify_brackets(s: str) -> str:
    """方括号统一为圆括号（方/圆仅括注外观不同，不判分）。"""
    return s.replace("[", "(").replace("]", ")")


def normalize_en(name: str) -> str:
    """规范化英文名：小写、去重空白并统一连字符/逗号/括号。"""
    s = (name or "").strip().lower()
    s = s.replace("–", "-").replace("—", "-")
    s = _unify_brackets(s)
    s = _WS.sub(" ", s)
    s = s.replace(" ,", ",")
    return s


def normalize_zh(name: str) -> str:
    """规范化中文名：去首尾空白并统一括号种类（方/圆等价，仅括注外观不同不判分）。"""
    return _unify_brackets((name or "").strip())


def nospace(name: str) -> str:
    """删去全部空白字符。"""
    return "".join((name or "").split())
