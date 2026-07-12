##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]
[#3352d9a][IUPAC P-29.3.1] 通用直链单烷基C1–C4侧链前缀贯通L3/L4/L5 [+10 tests, dual 0.5%→0.9%]
[#ff557f5][IUPAC P-61.3.1] 通用卤素前缀fluoro/chloro/bromo/iodo贯通L2/L3/L5 [+16 tests, dual 0.9%→1.2%]
[#pending][IUPAC P-14.2.1 / P-16.3] 倍数词头扩至penta–deca(五–十) [+9 tests, dual 1.2%(48)→1.2%(49), fails 4014→4013]

## 其他
- 基线(本轮前): en=1.2%(48/4062) zh=5.4%(41/756) dual=1.2%(48/4062)
- 本次: en=1.2%(49/4062) zh=5.4%(41/756) dual=1.2%(49/4062) fails=4013
- 改动: L5 MULT_EN/MULT_ZH 由2–4扩至5–10（penta…deca / 五…十）
- 验证: ClC(Cl)C(Cl)(Cl)Cl→1,1,1,2,2-pentachloroethane；FC(F)(F)C(F)(F)F→1,1,1,2,2,2-hexafluoroethane
- 已知缺口: 支化烷基(isopropyl)；卤代醇；更长/复杂多卤代仍缺母体与编号能力
- 下一步候选: 支化烷基(isopropyl/isobutyl)；或卤代醇；或编号引擎最低位次集边界
