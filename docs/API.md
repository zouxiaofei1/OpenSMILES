# API 参考

## SMILESNNamerV2

主命名器类，从 SMILES 生成 IUPAC 名称。

### 初始化

```python
from namepredict import SMILESNNamerV2

namer = SMILESNNamerV2(
    use_cache=True,      # 启用缓存
    timeout=5.0,         # 超时时间（秒）
)
```

### 主要方法

#### name()

生成化合物的中英文名称。

```python
result = namer.name(smiles: str) -> Dict[str, Any]
```

**参数**:
- `smiles` (str): 分子 SMILES 字符串

**返回**:
```python
{
    'en': str,           # 英文 IUPAC 名称
    'zh': str,           # 中文 IUPAC 名称
    'smiles': str,       # 标准化后的 SMILES
    'success': bool,     # 是否成功命名
    'confidence': float, # 置信度 (0-1)
    'source': str,       # 命名来源: 'iupac' | 'common' | 'cache'
    'type': str,         # 化合物类型
    'time_ms': float,    # 处理时间（毫秒）
}
```

**示例**:
```python
result = namer.name("CCO")
# {'en': 'ethanol', 'zh': '乙醇', 'success': True, ...}

result = namer.name("[Na+].[Cl-]")
# {'en': 'sodium chloride', 'zh': '氯化钠', 'success': True, ...}

result = namer.name("invalid_smiles")
# {'en': 'unknown compound', 'zh': '未知化合物', 'success': False, ...}
```

#### name_batch()

批量命名多个化合物。

```python
results = namer.name_batch(smiles_list: List[str]) -> List[Dict]
```

**参数**:
- `smiles_list` (List[str]): SMILES 字符串列表

**返回**: 字典列表，每个字典与 `name()` 返回值相同

**示例**:
```python
results = namer.name_batch(["CCO", "CC(C)O", "CCCCO"])
# [{'en': 'ethanol', ...}, {'en': 'propan-2-ol', ...}, {'en': 'pentan-1-ol', ...}]
```

### 辅助类

#### MoleculeParser

```python
from namepredict.core import MoleculeParser

parser = MoleculeParser()
mol_info = parser.parse("CCO")
# MoleculeInfo(atoms=[...], bonds=[...], rings=[], is_salt=False, ...)
```

#### ParentSelector

```python
from namepredict.core import ParentSelector

selector = ParentSelector()
result = selector.select(mol_info)
# ParentSelectionResult(parent_chain=[...], parent_type=..., ...)
```

#### FunctionalGroupAnalyzer

```python
from namepredict.core import FunctionalGroupAnalyzer

analyzer = FunctionalGroupAnalyzer()
fgs = analyzer.analyze(mol_info)
# [FunctionalGroup(type='hydroxy', atoms=[...], ...)]
```

### 常量

#### 元素名称

```python
from namepredict.constants import ELEMENT_NAMES

ELEMENT_NAMES['Na']
# {'en': 'sodium', 'zh': '钠'}

ELEMENT_NAMES['Fe']
# {'en': 'iron', 'zh': '铁'}
```

#### 官能团优先级

```python
from namepredict.constants import FG_PRIORITY

FG_PRIORITY['carboxylic_acid']  # 1 (最高)
FG_PRIORITY['alcohol']          # 5
FG_PRIORITY['ether']            # 10
```

#### 数字前缀

```python
from namepredict.constants import NUMERIC_PREFIXES

NUMERIC_PREFIXES[2]
# {'en': 'di', 'zh': '二'}

NUMERIC_PREFIXES[5]
# {'en': 'penta', 'zh': '五'}
```

## 错误处理

### 异常类型

```python
from namepredict import (
    InvalidSMILESError,  # SMILES 解析失败
    NamingError,         # 命名规则无法应用
    TimeoutError,        # 处理超时
)
```

### 示例

```python
try:
    result = namer.name(smiles)
    if not result['success']:
        print(f"无法命名: {result['en']}")
except InvalidSMILESError as e:
    print(f"无效 SMILES: {e}")
except TimeoutError:
    print("处理超时，分子可能过于复杂")
```

## 配置

### 环境变量

```bash
# 超时时间（秒）
NAMEPREDICT_TIMEOUT=10

# 缓存目录
NAMEPREDICT_CACHE_DIR=/path/to/cache

# 日志级别
NAMEPREDICT_LOG_LEVEL=DEBUG
```

### 运行时配置

```python
namer = SMILESNNamerV2()
namer.set_timeout(10.0)
namer.enable_cache(True)
namer.set_log_level('DEBUG')
```

## 性能

### 基准测试结果

| 分子类型 | 平均耗时 | 样本数 |
|---------|---------|--------|
| 简单烷烃 | 0.2 ms | 100 |
| 芳香族 | 0.5 ms | 100 |
| 杂环 | 1.2 ms | 100 |
| 配合物 | 2.5 ms | 50 |
| 复杂天然产物 | 10+ ms | 20 |

### 优化建议

1. **批量处理**: 使用 `name_batch()` 而非循环调用 `name()`
2. **缓存**: 启用缓存避免重复计算
3. **并行**: 对大量分子使用多进程

```python
from multiprocessing import Pool

def process_batch(smiles_list, n_workers=4):
    with Pool(n_workers) as pool:
        results = pool.map(namer.name, smiles_list)
    return results
```
