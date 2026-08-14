# Layer: L4
"""fg_locants 稀疏产出契约: principal FG 位次为结构化列表 [{kind, locants, omit}],
只产实际存在的 FG; 烯/炔不进列表(扁平字段独立); cooh 不产(死字段)."""
from __future__ import annotations

from namepredict.layer4.locant_calc import _fg_locants


def test_single_alcohol_record():
    fg = _fg_locants({"kind": "alcohol", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "oh_c_idx": 1, "oh_c_idxs": [1]})
    assert fg == [{"kind": "oh", "locants": [2], "omit": False}]


def test_methanol_omit_flag():
    fg = _fg_locants({"kind": "alcohol", "n_carbons": 1, "chain": [0],
                      "oh_c_idx": 0, "oh_c_idxs": [0]})
    assert fg == [{"kind": "oh", "locants": [1], "omit": True}]


def test_diol_single_oh_record_multi_locants():
    fg = _fg_locants({"kind": "diol", "n_carbons": 7, "chain": [9, 8, 6, 5, 3, 1, 0],
                      "oh_c_idxs": [1, 9]})
    assert fg == [{"kind": "oh", "locants": [1, 6], "omit": False}]


def test_ketone_record():
    fg = _fg_locants({"kind": "ketone", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "ketone_c_idx": 1, "ketone_c_idxs": [1]})
    assert fg == [{"kind": "ketone", "locants": [2], "omit": False}]


def test_dione_single_record():
    fg = _fg_locants({"kind": "dione", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "ketone_c_idxs": [1, 3]})
    assert fg == [{"kind": "ketone", "locants": [2, 4], "omit": False}]


def test_amine_record():
    fg = _fg_locants({"kind": "amine", "n_carbons": 4, "chain": [0, 1, 2, 3],
                      "amine_c_idx": 1, "amine_c_idxs": [1]})
    assert fg == [{"kind": "amine", "locants": [2], "omit": False}]


def test_acid_produces_no_record():
    """cooh 是死字段: 单/多酸位次隐含, 不产记录."""
    assert _fg_locants({"kind": "acid", "n_carbons": 5, "chain": [0, 1, 2, 3, 4],
                        "cooh_c_idx": 4, "cooh_c_idxs": [4]}) == []


def test_unsaturation_stays_out_of_records():
    """烯/炔不进 fg_locants: 不饱和度是扁平字段, 独立于 principal FG."""
    fg = _fg_locants({"kind": "alkene", "n_carbons": 2, "chain": [0, 1],
                      "double_bond": (0, 1)})
    assert fg == []
