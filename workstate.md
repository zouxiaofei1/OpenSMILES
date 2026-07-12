##日志
[#243ece8][IUPAC P-15.1.7 / 中文位次格式] 一元醇中文位次补连字符(丙-1-醇) [+3 zh dual, benchmark dual 0.3%→0.4%]
[#143333f][IUPAC P-44 母体链选择] 醇母体最长链穿过OH碳(propan-2-ol/butan-2-ol) [+5 dual, benchmark dual 0.4%→0.5%]
[#3352d9a][IUPAC P-29.3.1] 通用直链单烷基C1–C4侧链前缀贯通L3/L4/L5 [+10 tests, dual 0.5%→0.9%]
[#ff557f5][IUPAC P-61.3.1] 通用卤素前缀fluoro/chloro/bromo/iodo贯通L2/L3/L5 [+16 tests, dual 0.9%→1.2%]
[#4359897][IUPAC P-14.2.1 / P-16.3] 倍数词头扩至penta–deca(五–十) [+9 tests, dual 1.2%(48)→1.2%(49), fails 4014→4013]
[#9b6c9b5][IUPAC P-65.1.1] 一元羧酸：L1羧基检测+排除羧基OH；L2/L4/L5 acid母体与词干 [+8 tests, dual 1.2%(49)→1.5%(60), fails 4013→4002]
[#1bbe283][IUPAC P-64.2.1] 一元酮alkanone：L1酮羰基+L2/L4/L5 …-n-one/…-n-酮 [+11 tests, dual 1.5%(60)→1.7%(68), fails 4002→3994]
[#0a577d7][IUPAC P-66.6.1] 一元醛alkanal：L1醛羰基+L2/L4/L5 保留formaldehyde/acetaldehyde与…anal/…醛 [+10 tests, dual 1.7%(68)→1.8%(74), fails 3994→3988]
[#2b2c412][IUPAC P-31.1 / P-15.1.7.2.2] 一元烯烃alkene：L1非芳香C=C；L2穿双键母体；L4最低ene位次；L5 ane→ene/…-n-烯 [+12 tests, dual 1.8%(74)→1.9%(78), fails 3988→3984]

## 其他
- 基线(本轮前): en=1.8%(74/4062) zh=8.1%(61/756) dual=1.8%(74/4062)
- 本次: en=1.9%(78/4062) zh=8.6%(65/756) dual=1.9%(78/4062) fails=3984
- 改动: L1 double_bonds/has_alkene；L2 无更高FG且恰1个C=C时alkene母体；L4 ene_locant+C2/C3省略；L5 ethene/propene与but-n-ene
- 验证: C=C→ethene/乙烯；C=CC→propene/丙烯；CC=CC→but-2-ene/丁-2-烯；C=CCC→but-1-ene；CC(C)C=C→3-methylbut-1-ene；酸/醛/酮/醇/二烯/炔不回退
- 已知缺口: 二烯/炔/环烯；酯；支化烷基isopropyl；胺；羟基醛前缀
- 下一步候选: 一元炔烃；或酯；或支化烷基；或胺
