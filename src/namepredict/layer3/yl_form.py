"""Thin re-export of the free→yl word-stem conversion (P-29).

Word-stem knowledge lives in tools/free_to_yl.py; this module keeps the
historical layer3 import path for backward compatibility.
"""
from __future__ import annotations

from namepredict.tools.free_to_yl import free_to_yl as yl_form

__all__ = ["yl_form"]
