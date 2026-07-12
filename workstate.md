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
[#66f7f22][IUPAC P-66.1.1] 一元伯酰胺alkanamide：L1 –CONH2+排除醛；L2 acid>ester>amide>…；L5 formamide/acetamide/…amide [+16 tests, dual 2.4%(99)→2.5%(102), fails 3963→3960]
[#682a7d6][IUPAC P-66.5.1] 一元腈alkanenitrile：L1 C≡N；L2穿氰基碳；L4端基定向；L5 formonitrile/acetonitrile/…nitrile与…腈 [+15 tests, dual 2.5%(102)→2.6%(105), fails 3960→3957]
[#633ade7][IUPAC P-22.1.1 / P-14.3.4] 单取代环烷 monoalkyl-cycloalkane：L2侧链+L5省略位次 [+15 tests, dual 2.6%(105)→2.7%(109), fails 3957→3953]
[#490d44b][IUPAC P-63.1.1 / P-22.1.1] 未取代单环一元醇 cycloalkanol：L2环母体+L4省略位次+L5 cyclo…ol/环…醇 [+14 tests, dual 2.7%(109)→2.7%(111), fails 3953→3951]
[#7f8238a][IUPAC P-64.2.1 / P-22.1.1] 未取代单环一元酮 cycloalkanone：L2环母体+L4酮碳定向+L5 cyclo…one/环…酮 [+12 tests, dual 2.7%(111)→2.8%(113), fails 3951→3949]
[#5f87dd6][IUPAC P-61.3.1 / P-14.3.4 / P-22.1.1] 单卤代环烷 monohalo-cycloalkane：L2允许环上卤素+L3/L5省略位次 [+15 tests, dual 2.8%(113)→2.8%(114), fails 3949→3948]
[#7f436cb][IUPAC P-62.2.1 / P-22.1.1] 未取代单环一元伯胺 cycloalkanamine：L2环母体+L4胺碳定向+L5 cyclo…amine/环…胺 [+13 tests, dual 2.8%(114)→2.8%(115), fails 3948→3947]

## 其他
- 基线(本轮前): en=2.8%(114/4062) zh=12.7%(96/756) dual=2.8%(114/4062)
- 本次: en=2.8%(115/4062) zh=12.8%(97/756) dual=2.8%(115/4062) fails=3947
- 改动: L2 未取代单碳环+单伯胺→cycloamine；L4 胺碳定向并省略位次；L5 cyclo…amine/环…胺
- 验证: NC1CCCCC1→cyclohexanamine/环己胺；NC1CCCC1→cyclopentanamine；链胺/环醇/环酮/卤环负例保持
- 已知缺口: 多取代环烷；仲/叔胺；N-取代酰胺；支化烷基isopropyl；环烯烃
- 下一步候选: 支化烷基isopropyl；或仲胺；或多取代环烷；或环烯烃
