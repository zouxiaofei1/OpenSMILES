# Benchmark Ethanone命名不一致性分析

## 问题发现日期
2026-03-28

## 问题描述
在分析IUPAC P-64.2.2酮类位置编号规则时，发现benchmark数据库中存在命名格式不一致的问题。

## 具体案例

### 案例1：期望 `ethan-1-one`（带位置号）
- **ID**: 106998
- **SMILES**: BrCC(=O)C1=NC2=CC=CC=C2N=C1C
- **期望名称**: 2-bromo-1-(3-methylquinoxalin-2-yl)ethan-1-one
- **Tier**: 3

### 案例2：期望 `ethanone`（不带位置号）
- **ID**: 55163
- **SMILES**: BrCC(=O)C1CCNCC1
- **期望名称**: 2-Bromo-1-(piperidin-4-yl)ethanone
- **Tier**: 3

### 案例3：期望 `ethanone`（不带位置号）
- **ID**: 4138
- **SMILES**: BrCC(=O)C1=NC=CC(=C1)OC
- **期望名称**: 2-Bromo-1-(4-methoxypyridin-2-yl)ethanone
- **Tier**: 3

## 分析

### 共同点
- 都是乙酮衍生物
- 都在位置2有溴取代基
- 都在位置1连接杂环

### 不同点
- 案例1连接的是稠环(quinoxalin)
- 案例2连接的是单环(piperidin)
- 案例3连接的是单环(pyridin)

### 疑问
根据IUPAC P-14.4规则，当存在其他位置编号时，位置号应保留以保持一致性。但benchmark数据显示：
- 稠环连接时使用 `ethan-1-one`
- 单环连接时使用 `ethanone`

这种差异是否是IUPAC规则的要求，还是benchmark数据的不一致？

## 需要确认的IUPAC规则
1. P-14.4: 位置编号省略的具体条件
2. P-64.2.2: 酮类位置编号的具体规则
3. 是否有规则规定稠环连接时需要保留位置号

## 建议
如果这是benchmark数据的错误，需要统一命名格式。建议统一使用 `ethanone`（不带位置号），因为：
1. 乙酮的羰基碳只能位于位置1
2. 位置号是唯一的，可以省略
3. 简化命名

## 相关测试文件
- tests/test_p64_2_2_ketone_locant_format.py
- tests/test_ethanone_numbering.py
- tests/test_ketone_heterocycle_naming.py
