# 合并自 1 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_namer_smoke.py: 
"""
from __future__ import annotations

from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# ==========================================================================
# 合并自 test_namer_smoke.py
# ==========================================================================
def test_methane():
    r = SMILESNNamer().name("C")
    assert r.success
    assert normalize_en(r.en) == "methane"
    assert normalize_zh(r.zh) == "甲烷"


def test_ethanol():
    r = SMILESNNamer().name("CCO")
    assert r.success
    assert normalize_en(r.en) == "ethanol"
    assert normalize_zh(r.zh) == "乙醇"
