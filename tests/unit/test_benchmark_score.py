# 合并自 1 个原测试文件（按主题分组，内容与断言未改动）。
"""
test_benchmark_score.py: 
"""
from __future__ import annotations

from benchmarks.benchmark import _bucket_report, _build_snapshot, _change_similarity, _compare_to_previous, _diff_kind, _empty_bucket, _fmt_change_sim, _handle_report, _print_diffs, _print_summary, _tally, score_record

# ==========================================================================
# 合并自 test_benchmark_score.py
# ==========================================================================
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


def benchmark_score___entry(key, dual_ok, pred_en="a", pred_zh="甲", en_ok=True, zh_ok=True, smiles="C"):
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
            benchmark_score___entry("id:1", dual_ok=False, pred_en="wrong", pred_zh="错", en_ok=False, zh_ok=False, smiles="CCO"),
            benchmark_score___entry("id:2", dual_ok=True, pred_en="ethane", pred_zh="乙烷", en_ok=True, zh_ok=True, smiles="CC"),
            benchmark_score___entry("id:3", dual_ok=False, pred_en="methane", pred_zh="还错", en_ok=True, zh_ok=False, smiles="C"),
        ]
    }
    diffs = _compare_to_previous(report, prev)
    kinds = {d["key"]: d["kind"] for d in diffs}
    assert kinds["id:1"] == "REGRESS"
    assert kinds["id:2"] == "IMPROVE"
    assert kinds["id:3"] == "CHANGE"
    assert len(diffs) == 3


def test_compare_no_prev_empty():
    report = {"results": [benchmark_score___entry("id:1", True)]}
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
        "cur": benchmark_score___entry(
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
            "prev": benchmark_score___entry(key, False),
            "cur": benchmark_score___entry(key, True),
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
            benchmark_score___entry("id:1", dual_ok=True, pred_en="ethanol", pred_zh="乙醇", smiles="CCO"),
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
        "fails": [benchmark_score___entry("id:1", dual_ok=False, pred_en="x", pred_zh="y",
                         en_ok=False, zh_ok=False, smiles="CCO")],
        "results": [
            benchmark_score___entry("id:1", dual_ok=False, pred_en="x", pred_zh="y",
                   en_ok=False, zh_ok=False, smiles="CCO"),
        ],
    }
    _handle_report(report2, as_json=False, snapshot_path=snap)
    out2 = capsys.readouterr().out
    assert "REGRESS=1" in out2
    assert "[REGRESS] id:1" in out2


def benchmark_score___change(key, gold_en, gold_zh, prev_en, prev_zh, cur_en, cur_zh,
            en_ok=False, zh_ok=False):
    """一条 CHANGE diff: dual 未翻转, 只是预测串变了。"""
    def side(pe, pz):
        return {
            "dual_ok": False, "en_ok": en_ok, "zh_ok": zh_ok,
            "pred_en": pe, "pred_zh": pz,
            "smiles": "C", "english_name": gold_en, "chinese_name": gold_zh,
        }

    return {"kind": "CHANGE", "key": key, "prev": side(prev_en, prev_zh), "cur": side(cur_en, cur_zh)}


def test_change_similarity_none_without_change():
    assert _change_similarity([]) is None
    assert _change_similarity([{"kind": "IMPROVE"}, {"kind": "REGRESS"}]) is None


def test_change_similarity_closer_to_gold():
    # 丙-2-基氧基 -> 丙-2-氧基 那类改动: 准确率不动, 但字符串更贴金标。
    d = benchmark_score___change("id:1", "ethanol", "乙醇", "ethanoll", "乙纯", "ethanol", "乙醇")
    sim = _change_similarity([d])
    assert sim["n"] == 1
    assert sim["en"]["prev"] < sim["en"]["cur"] == 1.0
    assert sim["zh"]["prev"] < sim["zh"]["cur"] == 1.0
    assert "+" in _fmt_change_sim(sim)


def test_change_similarity_signed_delta_when_worse():
    d = benchmark_score___change("id:1", "ethanol", "乙醇", "ethanol", "乙醇", "ethanoll", "乙纯")
    sim = _change_similarity([d])
    assert sim["en"]["cur"] < sim["en"]["prev"]
    assert "-" in _fmt_change_sim(sim)


def test_change_similarity_skips_unevaluated_language():
    # 中文未考核的行不按打分口径外的标准评判: 均值和 n 都不含它们。
    d = benchmark_score___change("id:1", "ethanol", "乙醇", "ethanoll", "乙纯", "ethanol", "乙醇")
    d["prev"]["zh_ok"] = None
    d["cur"]["zh_ok"] = None
    sim = _change_similarity([d])
    assert "zh" not in sim
    assert sim["en"]["n"] == 1


def test_print_diffs_prints_change_sim_after_summary(capsys):
    _print_diffs([benchmark_score___change("id:1", "ethanol", "乙醇", "zzz", "甲", "ethanol", "乙醇")])
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith("diff_vs_last:")
    assert lines[1] == (
        "change_sim: 1 CHANGE  "
        "en 0.0%->100.0% (+100.00pt, n=1)  zh 0.0%->100.0% (+100.00pt, n=1)"
    )


def test_print_diffs_no_change_sim_line_without_change(capsys):
    diffs = [{
        "kind": "REGRESS", "key": "id:1",
        "prev": {"dual_ok": True, "en_ok": True, "zh_ok": True,
                 "pred_en": "ethanol", "pred_zh": "乙醇",
                 "smiles": "CCO", "english_name": "ethanol", "chinese_name": "乙醇"},
        "cur": benchmark_score___entry("id:1", dual_ok=False, pred_en="x", pred_zh="y",
                      en_ok=False, zh_ok=False, smiles="CCO"),
    }]
    _print_diffs(diffs)
    out = capsys.readouterr().out
    assert "diff_vs_last:" in out
    assert "change_sim" not in out


def test_build_snapshot_keys():
    report = {
        "ok_dual": 1, "n_dual": 1,
        "results": [benchmark_score___entry("id:9", dual_ok=True, smiles="C")],
    }
    snap = _build_snapshot(report)
    assert "id:9" in snap["items"]
    assert snap["items"]["id:9"]["dual_ok"] is True
    assert snap["ok_dual"] == 1
