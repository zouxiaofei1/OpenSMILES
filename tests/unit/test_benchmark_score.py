from benchmarks.benchmark import score_record, _bucket_report, _empty_bucket, _tally, _print_summary


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


def test_en_normalize_case_and_space():
    row = {"english_name": "Ethanol", "chinese_name": "", "eval_en": True, "eval_zh": False}
    r = score_record("  ethanol  ", "", row)
    assert r["en_ok"] is True
    assert r["dual_ok"] is True


def test_zh_strip_exact():
    row = {"english_name": "", "chinese_name": "乙醇", "eval_en": False, "eval_zh": True}
    r = score_record("", "  乙醇  ", row)
    assert r["zh_ok"] is True
    assert r["en_ok"] is None
    assert r["dual_ok"] is True


def test_neither_flag_dual_false():
    row = {"english_name": "x", "chinese_name": "y", "eval_en": False, "eval_zh": False}
    r = score_record("x", "y", row)
    assert r["en_ok"] is None
    assert r["zh_ok"] is None
    assert r["dual_ok"] is False


def test_bucket_report_keys():
    b = _empty_bucket()
    _tally(b, {"en_ok": True, "zh_ok": False, "dual_ok": False})
    rep = _bucket_report(b)
    for k in ("acc_en", "acc_zh", "acc_dual", "n_en", "n_zh", "n_dual"):
        assert k in rep
    assert rep["n_en"] == 1 and rep["n_zh"] == 1 and rep["n_dual"] == 1
    assert rep["acc_en"] == 1.0 and rep["acc_zh"] == 0.0 and rep["acc_dual"] == 0.0


def test_print_summary_format(capsys):
    report = {
        "acc_en": 0.65, "acc_zh": 0.6, "acc_dual": 0.6,
        "n_en": 20, "n_zh": 20, "n_dual": 20,
        "ok_en": 13, "ok_zh": 12, "ok_dual": 12,
        "fails": [1, 2, 3, 4, 5, 6, 7, 8],
    }
    _print_summary(report)
    out = capsys.readouterr().out.strip()
    assert out == "en=65.0% (13/20) zh=60.0% (12/20) dual=60.0% (12/20) fails=8"
