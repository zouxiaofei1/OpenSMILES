import time

from benchmarks.benchmark_parallel import (
    _build_parser,
    _print_timeouts,
    _run_row_with_timeout,
)


def _slow_score(_row):
    time.sleep(1)


def test_timeout_defaults_to_one_second():
    args = _build_parser().parse_args(["--data", "benchmark.json"])

    assert args.timeout == 1.0


def test_run_row_timeout_returns_index_without_printing(capsys):
    row = {
        "smiles": "SLOW",
        "english_name": "expected",
        "eval_en": True,
        "eval_zh": False,
    }

    result, timeout_index = _run_row_with_timeout(7, row, 0.05, score_func=_slow_score)

    assert result[1:] == ("", "", row)
    assert result[0]["en_ok"] is False
    assert timeout_index == 7
    assert capsys.readouterr().err == ""


def test_print_timeouts_combines_indexes_on_one_line(capsys):
    _print_timeouts([1033, 3268, 2944], 0.5)

    assert capsys.readouterr().err == "timeout index=1033,3268,2944 after=0.5s\n"
