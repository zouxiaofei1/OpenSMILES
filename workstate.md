##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]
[#3352d9a][IUPAC P-29.3.1] 通用直链单烷基C1–C4侧链前缀贯通L3/L4/L5 [+10 tests, dual 0.5%→0.9%]
[#ff557f5][IUPAC P-61.3.1] 通用卤素前缀fluoro/chloro/bromo/iodo贯通L2/L3/L5 [+16 tests, dual 0.9%→1.2%]

## 其他
- 基线(本轮前): en=0.9%(35/4062) zh=4.0%(30/756) dual=0.9%(35/4062)
- 本次: en=1.2%(48/4062) zh=5.4%(41/756) dual=1.2%(48/4062) fails=4014
- 改动: L3链碳卤素检测；L5甲烷/乙烷单位次省略；L2等长链侧基数tie-break
- 验证: CCCl→chloroethane；ClCCCl→1,2-dichloroethane；BrC(C)CC(C)C→2-bromo-4-methylpentane
- 已知缺口: MULT无penta(1,1,1,2,2-pentachloroethane→1,1,1,2,2-chloroethane)
- 下一步候选: 倍数词头扩至penta/hexa；或支化烷基(isopropyl)；或卤代醇
