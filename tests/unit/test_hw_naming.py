# IUPAC: P-22.2.2
# Layer: L2
"""Hantzsch-Widman 生成式命名核心：词干、位次引用序与位次省略。

用例元素序列均为「已编号环序」：第 i 项的位次为 i+1。
"""
from __future__ import annotations

import pytest

from namepredict.constants import (
    As, B, C, N, O, P, S, Se, Si,
)
from namepredict.layer2.hantzsch_widman import hw_name_from_cycle, omit_locants

# (环序元素表, 环内双键数, 期望英文名, 期望中文名)
CASES = [
    # ── 词干表 Table 2.5 ──
    ((O, C, C), 1, "oxirene", "氧杂环丙烯"),
    ((N, C, C), 1, "azirine", "氮杂环丙烯"),          # 仅含氮 → irine
    ((O, N, C), 1, "oxazirene", "氧杂氮杂环丙烯"),
    ((O, C, C, C), 0, "oxetane", "氧杂环丁烷"),
    ((N, C, C, C), 0, "azetidine", "氮杂环丁烷"),
    ((O, S, C, C, C), 0, "1,2-oxathiolane", "1,2-氧杂硫杂环戊烷"),
    ((O, N, C, C, P), 2, "1,2,5-oxazaphosphole", "1,2,5-噁磷杂唑"),
    ((S, C, C, C, C, C, C), 3, "thiepine", "硫杂环庚三烯"),
    ((O, C, C, C, C, C, C, C), 0, "oxocane", "氧杂环辛烷"),
    ((S, S, S, S, S, S, S, S), 0, "octathiocane", "八硫杂环辛烷"),
    # ── P-22.2.2.1.5.2：环中有氮则取 etidine/olidine ──
    ((O, S, N, C, C), 0, "1,2,3-oxathiazolidine", "1,2,3-噁噻唑烷"),
    # ── P-22.2.2.1.6：六元环词干由优先性最低的杂原子分组决定 ──
    ((O, C, C, O, C, C), 2, "1,4-dioxine", "1,4-二氧杂环己二烯"),   # 6A
    ((S, C, Se, C, C, C), 0, "1,3-thiaselenane", "1,3-硫杂硒杂环己烷"),
    ((N, C, N, C, N, C), 3, "1,3,5-triazine", "1,3,5-三嗪"),        # 6B
    ((O, C, N, C, C, C), 0, "1,3-oxazinane", "1,3-氧杂嗪烷"),
    ((Si, Si, Si, Si, Si, Si), 0, "hexasilinane", "六硅杂环己烷"),   # 6B，位次省略
    ((O, C, As, C, C, C), 0, "1,3-oxarsinane", "1,3-氧杂砷杂环己烷"),  # 6C
    ((P, B, P, B, P, B), 0, "1,3,5,2,4,6-triphosphatriborinane", "1,3,5,2,4,6-三磷杂三硼杂环己烷"),
    ((P, B, P, B, P, B), 3, "1,3,5,2,4,6-triphosphatriborinine", "1,3,5,2,4,6-三磷杂三硼杂环己三烯"),
    # ── P-22.2.2.1.3：位次按引用顺序排列，非数值升序 ──
    ((O, S, C, C, C, S, C), 0, "1,2,6-oxadithiepane", "1,2,6-氧杂二硫杂环庚烷"),
    ((O, N, C, C, C, O, C), 0, "1,6,2-dioxazepane", "1,6,2-二氧杂氮杂环庚烷"),
    ((O, N, C, C, N, C), 3, "1,2,5-oxadiazine", "1,2,5-氧杂二嗪"),
    ((O, N, C, C, C, N, C), 3, "1,2,6-oxadiazepine", "1,2,6-氧杂二氮杂环庚三烯"),  # 7 元含氮用长式
]

# (环序元素表, 是否省略位次)
OMIT_CASES = [
    ((O, C, C, C, C, C, C, C), True),        # 单杂原子：thiepine 不写 1-
    ((N, N, N, N, C), True),                 # 1H-tetrazole：4N 在 5 元环排布唯一
    ((Si, Si, Si, Si, Si, Si), True),        # hexasilinane
    ((N, C, N, C, N, C), False),             # 1,3,5-triazine：3N 有 3 种排布
    ((O, S, C, C, C), False),                # 1,2-oxathiolane：1,2 与 1,3 两种
    ((N, N, C, C, C), False),                # 1,2,4-oxadiazole 型：2N 在 5 元环有两种
]


@pytest.mark.parametrize("zs,n_db,en,zh", CASES)
def test_hw_name_from_cycle(zs, n_db, en, zh) -> None:
    r = hw_name_from_cycle(zs, n_db)
    assert r is not None
    assert r[0] == en
    assert r[1] == zh


@pytest.mark.parametrize("zs,expected", OMIT_CASES)
def test_omit_locants(zs, expected: bool) -> None:
    assert omit_locants(len(zs), zs) is expected


def test_locant_prefix_attached() -> None:
    """位次不可省时前缀带在名首（P-22.2.2.1.2）。"""
    r = hw_name_from_cycle((O, S, C, C, C), 0)
    assert r is not None
    assert r[2] == "1,2-"
    assert r[0].startswith("1,2-")


def test_single_heteroatom_has_no_locant_prefix() -> None:
    """单杂原子省略全部位次（P-22.2.2.1.7）。"""
    r = hw_name_from_cycle((O, C, C, C, C, C, C, C), 0)
    assert r is not None
    assert r[2] == ""
