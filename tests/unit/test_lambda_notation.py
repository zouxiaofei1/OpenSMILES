"""λ 约定基础设施：键数计算、非标准判定、记号渲染与判分折叠（P-14.1）。

Layer: L1（lambda_atoms 注入）/ tools
"""
from __future__ import annotations

from rdkit import Chem

from namepredict.tools.lambda_notation import (
    bonding_number, delta_mark, is_nonstandard, lambda_mark,
    locant_lambda_str, nonstandard_bonding, put_lambda,
)
from namepredict.tools.re import normalize_en, normalize_zh


def _mol(smiles):
    """解析 SMILES 并返回 mol。"""
    return Chem.MolFromSmiles(smiles)


def _nonstandard_by_symbol(smiles):
    """{元素符号: 键数}，便于断言。"""
    mol = _mol(smiles)
    return {mol.GetAtomWithIdx(i).GetSymbol(): n for i, n in nonstandard_bonding(mol).items()}


# ── 键数计算（P-14.1.1）──────────────
def test_bonding_number_basic():
    """PH3 的 P 键数为 3、H2S 的 S 为 2（标准价）。"""
    mol = _mol("P")
    assert bonding_number(mol.GetAtomWithIdx(0)) == 3
    mol = _mol("S")
    assert bonding_number(mol.GetAtomWithIdx(0)) == 2


def test_aromatic_heteroatom_is_standard():
    """噻吩/呋喃/吡啶的杂原子按交替单双键计，均为标准价，不标 λ。"""
    assert nonstandard_bonding(_mol("c1ccsc1")) == {}
    assert nonstandard_bonding(_mol("c1ccoc1")) == {}
    assert nonstandard_bonding(_mol("c1ccncc1")) == {}


# ── 非标准判定（P-14.1.2 表 1.3）──────────
def test_nonstandard_common_organics_empty():
    """普通有机物无一处非标准键数。"""
    for smi in ("CSC", "CCO", "c1ccccc1", "CC(=O)O", "CSSC"):
        assert nonstandard_bonding(_mol(smi)) == {}, smi


def test_nonstandard_phosphorus_lambda5():
    """五配位磷（P-67）：环磷酰胺的 P 键数为 5。"""
    assert _nonstandard_by_symbol("O=P1(NCCCl)OCCCN1") == {"P": 5}


def test_nonstandard_sulfur_lambda4_lambda6():
    """四配位/六配位硫：DMSO 得 λ4、砜得 λ6。"""
    assert _nonstandard_by_symbol("CS(C)=O") == {"S": 4}
    assert _nonstandard_by_symbol("CS(C)(=O)=O") == {"S": 6}


def test_nonstandard_iodine_lambda3():
    """三配位碘（二氯碘苯型 λ3-碘烷）超出标准价 1。"""
    assert _nonstandard_by_symbol("ClI(Cl)c1ccccc1") == {"I": 3}


def test_nonstandard_charged_atoms_excluded():
    """带电原子归 -ium/-ide 路径，不标 λ（铵、氧鎓、叠氮的 N⁺/N⁻）。"""
    assert nonstandard_bonding(_mol("[NH4+]")) == {}
    assert nonstandard_bonding(_mol("C[O+](C)C")) == {}


def test_is_nonstandard_single_atom():
    """is_nonstandard 与 nonstandard_bonding 口径一致。"""
    mol = _mol("CS(C)(=O)=O")
    s = next(a for a in mol.GetAtoms() if a.GetSymbol() == "S")
    assert is_nonstandard(s)


# ── 记号渲染（P-15.4.1.3：位次与 λ 间无连字符）──────
def test_lambda_mark_en_zh():
    """中文用 λ5、英文用测试集口径的 lambda5。"""
    assert lambda_mark(5) == "λ5"
    assert lambda_mark(5, en=True) == "lambda5"


def test_put_lambda_no_hyphen():
    """位次与 λ 连写：1 + 6 → 1λ6。"""
    assert put_lambda("1", 6) == "1λ6"
    assert put_lambda("1", 3, en=True) == "1lambda3"


def test_locant_lambda_str():
    """位次串按升序拼接，无 λ 的位次素写。"""
    assert locant_lambda_str([(1, 5), (2, None), (3, 5)]) == "1λ5,2,3λ5"
    assert locant_lambda_str([(1, 6), (2, None)]) == "1λ6,2"
    assert locant_lambda_str([(3, 4)], en=True) == "3lambda4"


def test_delta_mark():
    """δ 记号（P-25.7.2）：连续形式双键数 → δc。"""
    assert delta_mark(2) == "δ2"


# ── 判分折叠（normalize 上标 → ASCII）──────────
def test_normalize_folds_superscript():
    """λ⁵ 与 λ5 在中英文判分口径下等价。"""
    assert normalize_zh("1λ⁵,2-苯并噻嗪") == normalize_zh("1λ5,2-苯并噻嗪")
    assert normalize_en("1lambda3,2-benziodoxol") == normalize_en("1lambda3,2-benziodoxol")
    assert normalize_zh("9λ⁶-硫杂") == "9λ6-硫杂"
