"""L2 磷酸无机母体候选构建（kind=phosphate，P 中心单原子骨架）。"""
from __future__ import annotations


def build_phosphate_parent(det: dict, salt: dict | None) -> dict | None:
    """由 det（{n_oh,n_om,n_arms,p_idx}）构造磷酸母体 dict；不可命名的形态返回 None。

    盐门控：n_om>0 时若有碱金属须同数配对（可含烷基酯臂，如 (=O)([O-])([O-])OR 的
    磷酸单酯二钠盐）；完全无抗衡金属的游离磷酸根/磷酸酯阴离子也放行（如 (=O)([O-])
    ([O-])O-己基）。中性酸/酯不允许带金属；金属数目不足的欠中和形态不命名。
    """
    n_oh, n_om, n_arms = int(det.get("n_oh", 0)), int(det.get("n_om", 0)), int(det.get("n_arms", 0))
    p_idx = det.get("p_idx")
    if n_oh + n_om + n_arms != 3 or p_idx is None:
        return None
    salt = salt or {}
    metal = salt.get("metal")
    if n_om > 0:
        if metal and int(salt.get("n_metal") or 0) != n_om:
            return None
    elif metal:
        return None  # 中性形态带金属 = 电荷不平衡输入
    return {
        "kind": "phosphate", "chain": [int(p_idx)], "n_carbons": 0,
        "p_idx": int(p_idx), "n_oh": n_oh, "n_om": n_om, "n_arms": n_arms,
        "salt_meta": dict(salt) if salt else None,
        "principal_group_count": 1, "covered_principal_ids": (),
    }
