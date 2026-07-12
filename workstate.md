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
[#bd8354d][IUPAC P-31.1] 一元炔烃alkyne：L1 C≡C；L2穿三键且无C=C；L4最低yne位次；L5 acetylene/propyne与…-n-yne/…-n-炔 [+14 tests, dual 1.9%(78)→2.0%(82), fails 3984→3980]
[#2e8c0dd][IUPAC P-65.6] 一元酯alkyl alkanoate：L1酯检测+排除醛；L2 acid>ester>…；L5 alkyl+alkanoate/酸名+烷词干+酯 [+14 tests, dual 2.0%(82)→2.2%(90), fails 3980→3972]
[#4146028][IUPAC P-62.2.1] 一元伯胺alkanamine：L1伯胺N(排除酰胺)；L2 alcohol>amine；L4胺位次；L5 methanamine/…-n-amine与甲胺/…-n-胺 [+15 tests, dual 2.2%(90)→2.3%(95), fails 3972→3967]
[#acbcef4][IUPAC P-22.1.1] 单环烷烃cycloalkane：L1环信息；L2无取代单碳环母体；L5 cyclo+alkane/环+烷 [+18 tests, dual 2.3%(95)→2.4%(99), fails 3967→3963]

## 其他
- 基线(本轮前): en=2.3%(95/4062) zh=10.3%(78/756) dual=2.3%(95/4062)
- 本次: en=2.4%(99/4062) zh=10.8%(82/756) dual=2.4%(99/4062) fails=3963
- 改动: L1 rings/n_rings/has_ring；L2 恰1环全碳单键无侧链→cycloalkane；L5 cyclopropane…cyclodecane
- 验证: C1CC1→cyclopropane/环丙烷；C1CCCCC1→cyclohexane/环己烷；CC1CCCCC1≠cyclohexane；开链/FG不回退
- 已知缺口: 取代环烷；环醇/环酮/环烯；酰胺；仲/叔胺；支化烷基isopropyl
- 下一步候选: 酰胺；或取代环烷(甲基环己烷)；或仲胺；或支化烷基
