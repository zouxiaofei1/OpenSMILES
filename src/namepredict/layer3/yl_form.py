"""Free PIN → P-29 -yl dual form at attach locant.

Rules (fixed for recursive substituent engine):
| free EN              | attach | yl EN                          |
|----------------------|--------|--------------------------------|
| benzene              | 1      | phenyl (retained shortcut)     |
| pyridine             | k      | pyridin-k-yl                   |
| 1,3-thiazole         | k      | 1,3-thiazol-k-yl               |
| …ole / …ine etc.     | k      | drop terminal e (if any) + -k-yl |
| free name w/ prefix  | k      | keep prefix + parent → yl      |

Chinese: 吡啶 → 吡啶-2-基; 1,3-噻唑 → 1,3-噻唑-2-基; 苯 → 苯基.
"""
from __future__ import annotations

# free_en → (yl_en, yl_zh) when locant is conventional / omitted
_RETAINED: dict[str, tuple[str, str]] = {
    "benzene": ("phenyl", "苯基"),
}


def _drop_terminal_e(en: str) -> str:
    return en[:-1] if en.endswith("e") else en


def _yl_en(en: str, k: int) -> str:
    hit = _RETAINED.get(en)
    return hit[0] if hit is not None else f"{_drop_terminal_e(en)}-{k}-yl"


def _yl_zh(zh: str, k: int, en: str) -> str:
    hit = _RETAINED.get(en)
    return hit[1] if hit is not None else f"{zh}-{k}-基"


def yl_form(
    en: str, zh: str, attach_locant: int, *, paren: bool = True,
) -> tuple[str, str, bool]:
    """Convert free parent name to P-29 -yl dual form at attach locant."""
    return _yl_en(en, attach_locant), _yl_zh(zh, attach_locant, en), paren
