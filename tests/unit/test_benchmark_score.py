from benchmarks.benchmark import (
    score_record,
    _bucket_report,
    _empty_bucket,
    _tally,
    _print_summary,
    _build_snapshot,
    _compare_to_previous,
    _diff_kind,
    _print_diffs,
    _handle_report,
)


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


# 括号种类(圆/方)仅是括注外观，不判分：中文 1-氯-3-(氯(苯基)甲基)苯 与
# 1-氯-3-[氯(苯基)甲基]苯 视为相同；英文 normalize_en 早已折叠，一并固化。
def test_zh_bracket_kind_ignored():
    row = {
        "english_name": "",
        "chinese_name": "1-氯-3-(氯(苯基)甲基)苯",
        "eval_en": False,
        "eval_zh": True,
    }
    r = score_record("", "1-氯-3-[氯(苯基)甲基]苯", row)
    assert r["zh_ok"] is True
    assert r["dual_ok"] is True


def test_en_bracket_kind_ignored():
    row = {
        "english_name": "1-chloro-3-[chloro(phenyl)methyl]benzene",
        "chinese_name": "",
        "eval_en": True,
        "eval_zh": False,
    }
    r = score_record("1-chloro-3-(chloro(phenyl)methyl)benzene", "", row)
    assert r["en_ok"] is True
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


def _entry(key, dual_ok, pred_en="a", pred_zh="甲", en_ok=True, zh_ok=True, smiles="C"):
    return {
        "key": key,
        "id": key.replace("id:", ""),
        "smiles": smiles,
        "english_name": "gold_en",
        "chinese_name": "金标",
        "pred_en": pred_en,
        "pred_zh": pred_zh,
        "en_ok": en_ok,
        "zh_ok": zh_ok,
        "dual_ok": dual_ok,
    }


def test_diff_kind_regress_improve_change():
    ok = {"dual_ok": True, "en_ok": True, "zh_ok": True, "pred_en": "a", "pred_zh": "甲"}
    bad = {"dual_ok": False, "en_ok": False, "zh_ok": True, "pred_en": "b", "pred_zh": "甲"}
    bad2 = {"dual_ok": False, "en_ok": True, "zh_ok": False, "pred_en": "c", "pred_zh": "乙"}
    assert _diff_kind(ok, bad) == "REGRESS"
    assert _diff_kind(bad, ok) == "IMPROVE"
    assert _diff_kind(bad, bad2) == "CHANGE"
    assert _diff_kind(ok, ok) is None


def test_compare_to_previous_detects_diffs():
    prev = {
        "items": {
            "id:1": {
                "dual_ok": True, "en_ok": True, "zh_ok": True,
                "pred_en": "ethanol", "pred_zh": "乙醇",
                "smiles": "CCO", "english_name": "ethanol", "chinese_name": "乙醇",
            },
            "id:2": {
                "dual_ok": False, "en_ok": False, "zh_ok": False,
                "pred_en": "x", "pred_zh": "y",
                "smiles": "CC", "english_name": "ethane", "chinese_name": "乙烷",
            },
            "id:3": {
                "dual_ok": False, "en_ok": True, "zh_ok": False,
                "pred_en": "methane", "pred_zh": "错",
                "smiles": "C", "english_name": "methane", "chinese_name": "甲烷",
            },
        }
    }
    report = {
        "results": [
            _entry("id:1", dual_ok=False, pred_en="wrong", pred_zh="错", en_ok=False, zh_ok=False, smiles="CCO"),
            _entry("id:2", dual_ok=True, pred_en="ethane", pred_zh="乙烷", en_ok=True, zh_ok=True, smiles="CC"),
            _entry("id:3", dual_ok=False, pred_en="methane", pred_zh="还错", en_ok=True, zh_ok=False, smiles="C"),
        ]
    }
    diffs = _compare_to_previous(report, prev)
    kinds = {d["key"]: d["kind"] for d in diffs}
    assert kinds["id:1"] == "REGRESS"
    assert kinds["id:2"] == "IMPROVE"
    assert kinds["id:3"] == "CHANGE"
    assert len(diffs) == 3


def test_compare_no_prev_empty():
    report = {"results": [_entry("id:1", True)]}
    assert _compare_to_previous(report, None) == []


def test_print_diffs_format(capsys):
    diffs = [{
        "kind": "REGRESS",
        "key": "id:1",
        "prev": {
            "dual_ok": True, "en_ok": True, "zh_ok": True,
            "pred_en": "ethanol", "pred_zh": "乙醇",
            "smiles": "CCO", "english_name": "ethanol", "chinese_name": "乙醇",
        },
        "cur": _entry(
            "id:1", dual_ok=False, pred_en="x", pred_zh="y",
            en_ok=False, zh_ok=False, smiles="CCO",
        ),
    }]
    _print_diffs(diffs)
    out = capsys.readouterr().out
    assert "REGRESS=1" in out
    assert "[REGRESS] id:1" in out
    assert "smiles=CCO" in out


def test_print_diffs_groups_improve_regress_change(capsys):
    # 三类并存时按 IMPROVE -> REGRESS -> CHANGE 分组显示，组内保持原顺序。
    def d(kind, key):
        return {
            "kind": kind, "key": key,
            "prev": _entry(key, False),
            "cur": _entry(key, True),
        }

    diffs = [
        d("CHANGE", "id:c1"),
        d("IMPROVE", "id:i1"),
        d("REGRESS", "id:r1"),
        d("IMPROVE", "id:i2"),
        d("CHANGE", "id:c2"),
    ]
    _print_diffs(diffs)
    out = capsys.readouterr().out
    assert "IMPROVE=2 REGRESS=1 CHANGE=2" in out
    assert out.index("[IMPROVE]") < out.index("[REGRESS]") < out.index("[CHANGE]")
    assert out.index("[CHANGE] id:c1") < out.index("[CHANGE] id:c2")  # 组内稳定


def test_print_diffs_none(capsys):
    _print_diffs([])
    out = capsys.readouterr().out
    assert "diff_vs_last: none" in out


def test_handle_report_saves_and_diffs(tmp_path, capsys):
    snap = tmp_path / "last.json"
    # First run: no previous snapshot → no diffs, then save.
    report1 = {
        "acc_en": 1.0, "acc_zh": 1.0, "acc_dual": 1.0,
        "n_en": 1, "n_zh": 1, "n_dual": 1,
        "ok_en": 1, "ok_zh": 1, "ok_dual": 1,
        "fails": [],
        "results": [
            _entry("id:1", dual_ok=True, pred_en="ethanol", pred_zh="乙醇", smiles="CCO"),
        ],
    }
    _handle_report(report1, as_json=False, snapshot_path=snap)
    out1 = capsys.readouterr().out
    assert "diff_vs_last: none" in out1
    assert snap.is_file()

    # Second run: regress on id:1 → printed.
    report2 = {
        "acc_en": 0.0, "acc_zh": 0.0, "acc_dual": 0.0,
        "n_en": 1, "n_zh": 1, "n_dual": 1,
        "ok_en": 0, "ok_zh": 0, "ok_dual": 0,
        "fails": [_entry("id:1", dual_ok=False, pred_en="x", pred_zh="y",
                         en_ok=False, zh_ok=False, smiles="CCO")],
        "results": [
            _entry("id:1", dual_ok=False, pred_en="x", pred_zh="y",
                   en_ok=False, zh_ok=False, smiles="CCO"),
        ],
    }
    _handle_report(report2, as_json=False, snapshot_path=snap)
    out2 = capsys.readouterr().out
    assert "REGRESS=1" in out2
    assert "[REGRESS] id:1" in out2


def test_build_snapshot_keys():
    report = {
        "ok_dual": 1, "n_dual": 1,
        "results": [_entry("id:9", dual_ok=True, smiles="C")],
    }
    snap = _build_snapshot(report)
    assert "id:9" in snap["items"]
    assert snap["items"]["id:9"]["dual_ok"] is True
    assert snap["ok_dual"] == 1
