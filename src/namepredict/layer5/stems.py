"""Carbon-count stem tables for C1–C10 alkanes, monoalcohols, monoacids."""

from __future__ import annotations

ALKANE_EN = {
    1: "methane",
    2: "ethane",
    3: "propane",
    4: "butane",
    5: "pentane",
    6: "hexane",
    7: "heptane",
    8: "octane",
    9: "nonane",
    10: "decane",
}
ALKANE_ZH = {
    1: "甲烷",
    2: "乙烷",
    3: "丙烷",
    4: "丁烷",
    5: "戊烷",
    6: "己烷",
    7: "庚烷",
    8: "辛烷",
    9: "壬烷",
    10: "癸烷",
}
ALCOHOL_EN = {
    1: "methanol",
    2: "ethanol",
    3: "propanol",
    4: "butanol",
    5: "pentanol",
    6: "hexanol",
    7: "heptanol",
    8: "octanol",
    9: "nonanol",
    10: "decanol",
}
ALCOHOL_ZH = {
    1: "甲醇",
    2: "乙醇",
    3: "丙醇",
    4: "丁醇",
    5: "戊醇",
    6: "己醇",
    7: "庚醇",
    8: "辛醇",
    9: "壬醇",
    10: "癸醇",
}
# P-65.1.1 retained: formic/acetic; else systematic …oic acid / …酸
ACID_EN = {
    1: "formic acid",
    2: "acetic acid",
    3: "propanoic acid",
    4: "butanoic acid",
    5: "pentanoic acid",
    6: "hexanoic acid",
    7: "heptanoic acid",
    8: "octanoic acid",
    9: "nonanoic acid",
    10: "decanoic acid",
}
ACID_ZH = {
    1: "甲酸",
    2: "乙酸",
    3: "丙酸",
    4: "丁酸",
    5: "戊酸",
    6: "己酸",
    7: "庚酸",
    8: "辛酸",
    9: "壬酸",
    10: "癸酸",
}
# P-66.6.1 retained: formaldehyde/acetaldehyde; else systematic …anal / …醛
ALDEHYDE_EN = {
    1: "formaldehyde",
    2: "acetaldehyde",
    3: "propanal",
    4: "butanal",
    5: "pentanal",
    6: "hexanal",
    7: "heptanal",
    8: "octanal",
    9: "nonanal",
    10: "decanal",
}
ALDEHYDE_ZH = {
    1: "甲醛",
    2: "乙醛",
    3: "丙醛",
    4: "丁醛",
    5: "戊醛",
    6: "己醛",
    7: "庚醛",
    8: "辛醛",
    9: "壬醛",
    10: "癸醛",
}
# P-66.1.1 retained: formamide/acetamide; else systematic …amide / …酰胺
AMIDE_EN = {
    1: "formamide",
    2: "acetamide",
    3: "propanamide",
    4: "butanamide",
    5: "pentanamide",
    6: "hexanamide",
    7: "heptanamide",
    8: "octanamide",
    9: "nonanamide",
    10: "decanamide",
}
AMIDE_ZH = {
    1: "甲酰胺",
    2: "乙酰胺",
    3: "丙酰胺",
    4: "丁酰胺",
    5: "戊酰胺",
    6: "己酰胺",
    7: "庚酰胺",
    8: "辛酰胺",
    9: "壬酰胺",
    10: "癸酰胺",
}
# P-66.5.1 retained acetonitrile; C1 formonitrile; else …nitrile / …腈
NITRILE_EN = {
    1: "formonitrile",
    2: "acetonitrile",
    3: "propanenitrile",
    4: "butanenitrile",
    5: "pentanenitrile",
    6: "hexanenitrile",
    7: "heptanenitrile",
    8: "octanenitrile",
    9: "nonanenitrile",
    10: "decanenitrile",
}
NITRILE_ZH = {
    1: "甲腈",
    2: "乙腈",
    3: "丙腈",
    4: "丁腈",
    5: "戊腈",
    6: "己腈",
    7: "庚腈",
    8: "辛腈",
    9: "壬腈",
    10: "癸腈",
}
# P-65.6 functional class: alkyl alkanoate (retained formate/acetate)
ESTER_ACYL_EN = {
    1: "formate",
    2: "acetate",
    3: "propanoate",
    4: "butanoate",
    5: "pentanoate",
    6: "hexanoate",
    7: "heptanoate",
    8: "octanoate",
    9: "nonanoate",
    10: "decanoate",
}
# Alkyl stem for ester (no 基): methyl→甲, ethyl→乙, ...
ESTER_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
ESTER_ALKYL_ZH = {1: "甲", 2: "乙", 3: "丙", 4: "丁"}

# P-65.5.1 retained acetyl; C>=3 alkane stem -e + oyl chloride / …酰氯
ACYL_CHLORIDE_EN = {
    2: "acetyl chloride",
    3: "propanoyl chloride",
    4: "butanoyl chloride",
    5: "pentanoyl chloride",
    6: "hexanoyl chloride",
    7: "heptanoyl chloride",
    8: "octanoyl chloride",
    9: "nonanoyl chloride",
    10: "decanoyl chloride",
}
ACYL_CHLORIDE_ZH = {
    2: "乙酰氯",
    3: "丙酰氯",
    4: "丁酰氯",
    5: "戊酰氯",
    6: "己酰氯",
    7: "庚酰氯",
    8: "辛酰氯",
    9: "壬酰氯",
    10: "癸酰氯",
}
# P-63.2.2 symmetric dialkyl ether retained: dimethyl…dibutyl ether / 二…基醚
ETHER_SYM_EN = {
    1: "dimethyl ether", 2: "diethyl ether",
    3: "dipropyl ether", 4: "dibutyl ether",
}
ETHER_SYM_ZH = {
    1: "二甲基醚", 2: "二乙基醚", 3: "二丙基醚", 4: "二丁基醚",
}
# Asymmetric alkoxy prefix (alkoxyalkane)
ALKOXY_EN = {1: "methoxy", 2: "ethoxy", 3: "propoxy", 4: "butoxy"}
ALKOXY_ZH = {1: "甲氧基", 2: "乙氧基", 3: "丙氧基", 4: "丁氧基"}
# P-63.2.1 dialkyl sulfide: dimethyl…dibutyl sulfide / 二…硫醚
SULFIDE_SYM_EN = {
    1: "dimethyl sulfide", 2: "diethyl sulfide",
    3: "dipropyl sulfide", 4: "dibutyl sulfide",
}
SULFIDE_SYM_ZH = {
    1: "二甲硫醚", 2: "二乙硫醚", 3: "二丙硫醚", 4: "二丁硫醚",
}
# Alkyl radical for asymmetric alkyl alkyl sulfide (alphabetical EN)
SULFIDE_ALKYL_EN = {1: "methyl", 2: "ethyl", 3: "propyl", 4: "butyl"}
SULFIDE_ALKYL_ZH = {1: "甲基", 2: "乙基", 3: "丙基", 4: "丁基"}

