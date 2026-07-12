# 架构设计

## 概述

NamePredict V2 采用模块化架构，将命名过程分解为多个独立阶段：

```
SMILES 输入
    ↓
[MoleculeParser] 分子解析
    ↓
[ParentSelector] 母体选择
    ↓
[FunctionalGroupAnalyzer] 官能团分析
    ↓
[NumberingEngine] 编号
    ↓
[SubstituentAnalyzer] 取代基分析
    ↓
[NameAssembler] 名称组装
    ↓
输出: {en, zh}
```

## 核心模块

### 1. MoleculeParser（分子解析器）

**职责**: 解析 SMILES 字符串，提取分子结构信息

**输入**: SMILES 字符串
**输出**: `MoleculeInfo` 对象

```python
class MoleculeInfo:
    mol: Chem.Mol           # RDKit 分子对象
    atoms: List[AtomInfo]   # 原子信息
    bonds: List[BondInfo]   # 键信息
    rings: List[RingInfo]   # 环信息
    is_salt: bool           # 是否为盐
    components: List[str]   # 组分 SMILES
```

**关键方法**:
- `parse(smiles: str) -> MoleculeInfo`
- `_detect_salt_components()`
- `_normalize_smiles()`

### 2. ParentSelector（母体选择器）

**职责**: 根据官能团优先级选择母体结构

**输入**: `MoleculeInfo`
**输出**: `ParentSelectionResult`

```python
class ParentSelectionResult:
    parent_chain: List[int]      # 母链原子索引
    parent_type: ParentType      # 母体类型
    principal_fg: FunctionalGroup  # 主官能团
    secondary_fgs: List[FunctionalGroup]  # 次要官能团
```

**选择规则**:
1. 官能团优先级（蓝皮书 P-44）
2. 最长链原则
3. 最大取代基数原则

### 3. FunctionalGroupAnalyzer（官能团分析器）

**职责**: 识别和分类分子中的官能团

**支持类型**:
- 含氧官能团（羟基、羰基、羧基等）
- 含氮官能团（氨基、硝基、腈基等）
- 含硫官能团（巯基、磺酸基等）
- 杂环（吡啶、呋喃、噻吩等）

**优先级排序**:
```
羧酸 > 酯 > 酰卤 > 酰胺 > 腈 > 醛 > 酮 > 醇 > 胺 > 醚
```

### 4. NumberingEngine（编号引擎）

**职责**: 为母体和取代基分配编号

**编号原则**:
1. 主官能团获得最小编号
2. 取代基编号和最小
3. 多重键优先于单键

**输出**:
```python
class NumberedSubstituent:
    locants: List[int]     # 位次号
    name: str              # 取代基名称
    prefix: str            # 前缀
    suffix: str            # 后缀（如有）
```

### 5. SubstituentAnalyzer（取代基分析器）

**职责**: 识别、命名和排序取代基

**支持类型**:
- 烷基（甲基、乙基、异丙基等）
- 卤素（氟、氯、溴、碘）
- 芳基（苯基、苄基等）
- 复杂取代基（支链、环状）

### 6. NameAssembler（名称组装器）

**职责**: 组装最终的中英文 IUPAC 名称

**组装规则**:
1. 取代基按字母顺序排列（英文）
2. 中文按笔画数或拼音排序
3. 正确添加连字符和括号
4. 处理多重前缀

## 检测器 Mixin

无机物和特殊化合物通过 Mixin 类扩展功能：

```
_detect_inorganic_mixin.py    # 无机物检测
_detect_salt_mixin.py          # 盐类检测
_detect_oxygen_mixin.py        # 含氧化合物
_detect_nitrogen_mixin.py      # 含氮化合物
_detect_sulfur_mixin.py        # 含硫化合物
_detect_phosphorus_mixin.py    # 含磷化合物
_detect_cyanate_mixin.py       # 氰酸盐类
_detect_quinone_mixin.py       # 醌类
```

## 缓存系统

### CommonNameCache

存储常见化合物的通用名称：

```python
# 优先使用通用名
"toluene" > "methylbenzene"
"acetone" > "propan-2-one"
"chloroform" > "trichloromethane"
```

### ExceptionCache

处理 IUPAC 规则例外：

```python
# 特殊命名
"HCOOH" -> "formic acid"  # 甲酸
"CH3COOH" -> "acetic acid"  # 乙酸
```

## 数据流示例

### 乙醇 (CCO)

```
1. MoleculeParser: 识别 2 个碳原子、1 个羟基
2. ParentSelector: 选择乙烷为母体
3. FunctionalGroupAnalyzer: 识别羟基为唯一官能团
4. NumberingEngine: 羟基获得编号 1
5. SubstituentAnalyzer: 无取代基
6. NameAssembler: 
   - en: "ethanol"
   - zh: "乙醇"
```

### 4-氯甲苯 (Cc1ccc(Cl)cc1)

```
1. MoleculeParser: 识别苯环、甲基、氯原子
2. ParentSelector: 选择甲苯为母体
3. FunctionalGroupAnalyzer: 识别甲基、氯
4. NumberingEngine: 
   - 甲基在 1 位（甲苯母体）
   - 氯在 4 位
5. SubstituentAnalyzer: 氯为取代基
6. NameAssembler:
   - en: "4-chlorotoluene"
   - zh: "4-氯甲苯"
```

## 扩展指南

### 添加新官能团

1. 在 `FunctionalGroupAnalyzer` 中添加检测逻辑
2. 在 `constants.py` 中定义优先级
3. 在 `NameAssembler` 中添加命名规则
4. 添加单元测试

### 添加新无机物类型

1. 在 `_detect_inorganic_mixin.py` 中添加检测方法
2. 实现命名逻辑
3. 添加测试用例到 `tests-inorganic/`

## 性能优化

- **缓存**: 常见分子解析结果缓存
- **惰性计算**: 仅在需要时计算昂贵操作
- **并行处理**: 多组分分子并行处理

## 错误处理

```python
try:
    result = namer.name(smiles)
except InvalidSMILESError:
    # SMILES 解析失败
except NamingError:
    # 命名规则无法应用
except TimeoutError:
    # 超时（复杂分子）
```
