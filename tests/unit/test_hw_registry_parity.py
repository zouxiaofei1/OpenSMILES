# IUPAC: P-22.2.2
# Layer: L2,L4,L5
"""已被生成式 Hantzsch-Widman 接管的原模板条目。

这些骨架的名字可由 P-22.2.2 生成器逐字复现（英文与中文皆同），
故从 `_TEMPLATES` 移除，避免模板与生成器两个事实来源。
"""
from __future__ import annotations

import pytest

from namepredict.layer2.ring_scaffold import _TEMPLATES
from namepredict.namer import SMILESNNamer
from namepredict.tools.re import normalize_en, normalize_zh

# (sid, 模板 SMILES, 期望英文, 期望中文)
TAKEOVER = [
    ("aziridine", "C1CN1", "aziridine", "氮杂环丙烷"),
    ("oxetane", "C1COC1", "oxetane", "氧杂环丁烷"),
    ("azetidine", "C1CNC1", "azetidine", "氮杂环丁烷"),
    ("oxepane", "O1CCCCCC1", "oxepane", "氧杂环庚烷"),
    ("oxazepane", "O1CCNCCC1", "1,4-oxazepane", "1,4-氧杂氮杂环庚烷"),
    ("thiazepane", "S1CCNCCC1", "1,4-thiazepane", "1,4-硫杂氮杂环庚烷"),
    ("oxazolidine", "C1NCCO1", "1,3-oxazolidine", "1,3-噁唑烷"),
    ("thiazolidine", "C1NCCS1", "1,3-thiazolidine", "1,3-噻唑烷"),
    ("thiadiazolidine124", "S1NCNC1", "1,2,4-thiadiazolidine", "1,2,4-噻二唑烷"),
    ("dioxaborolane132", "B1OCCO1", "1,3,2-dioxaborolane", "1,3,2-二氧杂硼杂环戊烷"),
    ("oxazinane13", "C1CNCOC1", "1,3-oxazinane", "1,3-氧杂嗪烷"),
    ("diazinane13", "C1CNCNC1", "1,3-diazinane", "1,3-二嗪烷"),
]


@pytest.mark.parametrize("sid,smiles,en,zh", TAKEOVER)
def test_generator_reproduces_removed_template(sid: str, smiles: str, en: str, zh: str) -> None:
    """移除后仍能由生成器给出同一名字（不得回退成碳环名或失败）。"""
    r = SMILESNNamer().name(smiles)
    assert r.success, r
    assert normalize_en(r.en) == normalize_en(en)
    assert normalize_zh(r.zh) == normalize_zh(zh)


@pytest.mark.parametrize("sid,smiles,en,zh", TAKEOVER)
def test_removed_templates_not_reintroduced(sid: str, smiles: str, en: str, zh: str) -> None:
    """已接管的条目不得回加到 _TEMPLATES（唯一事实来源）。"""
    assert sid not in _TEMPLATES


def test_locant_reference_order_for_oxadiazole() -> None:
    """P-22.2.2.1.3：位次 1 给引用序最前者，故 1,3,4- 而非 1,2,4-（N 不抢 O 的 1 位）。"""
    r = SMILESNNamer().name("o1cnnc1")
    assert r.success, r
    assert normalize_en(r.en) == normalize_en("1,3,4-oxadiazole")
