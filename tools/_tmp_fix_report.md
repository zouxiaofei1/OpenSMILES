## 小结

### 修改文件
1. `E:\dev\chem\src\namepredict\layer1\analyzer.py`
2. `E:\dev\chem\tests\unit\test_dialkyl_sulfide.py`

### 关键改动

**`_is_sulfide_sulfur`**（对齐醚逻辑）：
```python
def _is_sulfide_sulfur(atom) -> bool:
    if atom.GetAtomicNum() != 16 or atom.GetTotalNumHs() != 0:
        return False
    if atom.GetTotalDegree() != 2:
        return False
    cs = _sulfide_cs(atom)
    return len(cs) == 2 and not any(_has_double_bonded_o(c) for c in cs)
```

- `GetTotalDegree() == 2`：排除亚砜等更高配位硫（如 `CS(C)=O`）
- `not any(_has_double_bonded_o(c) for c in cs)`：排除硫酯（如 `CC(=O)SC`）

**负例测试** `NEG_NOT_SULFIDE`：
- `CS(C)=O` 不得为 dimethyl sulfide
- `CC(=O)SC` 不得含 “sulfide” / “硫醚”

### 验证结果
- **structure_lint**: ok
- **pytest**: 19 passed in 0.48s  
  (`test_dialkyl_sulfide.py` + `test_alkanethiol.py`)