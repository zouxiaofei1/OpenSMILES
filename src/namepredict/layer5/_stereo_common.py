"""Shared stereo-block splitting for layer5 name assembly.

Extracted verbatim from layer5/stereo_rs.py (_strip_stereo) and
layer5/benzene_names.py (_stereo_lead), which duplicated this helper.
"""
from __future__ import annotations


def _split_stereo_lead(name: str) -> tuple[str, str]:
    """Split leading '(…)-' stereo block from a name/stem."""
    if not name.startswith("("):
        return "", name
    close = name.find(")-")
    if close < 0:
        return "", name
    return name[: close + 2], name[close + 2 :]
