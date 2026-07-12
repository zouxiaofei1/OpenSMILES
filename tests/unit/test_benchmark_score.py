from benchmarks.benchmark import score_record


def test_dual_only_en_when_no_zh():
    row = {"english_name": "ethanol", "chinese_name": "", "eval_en": True, "eval_zh": False}
    r = score_record("ethanol", "", row)
    assert r["en_ok"] is True
    assert r["zh_ok"] is None
    assert r["dual_ok"] is True


def test_dual_requires_both_when_flags():
    row = {"english_name": "ethanol", "chinese_name": "乙醇", "eval_en": True, "eval_zh": True}
    r = score_record("ethanol", "酒精", row)
    assert r["en_ok"] is True
    assert r["zh_ok"] is False
    assert r["dual_ok"] is False
