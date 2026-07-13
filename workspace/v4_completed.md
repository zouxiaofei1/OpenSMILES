##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]

## 其他
- 基线(本轮前): en=0.4%(16/4062) zh=2.0%(15/756) dual=0.4%(16/4062)
- 本次: en=0.5%(21/4062) zh=2.5%(19/756) dual=0.5%(21/4062)
- 改动: layer2 新增 `_arms_from`/`_join_through`，`_chain_with_oh` 取穿过 oh_c 的两最长臂
- 验证: CC(O)C→propan-2-ol；CC(O)CC→butan-2-ol；伯醇/甲醇/乙醇不变
- 下一步候选: layer3 提取甲基等简单烷基取代基 → 2-methylpropan-*-ol
