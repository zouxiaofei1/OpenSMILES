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
