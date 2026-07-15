"""L5 dispatch for special functional-parent names (thin assembler helper)."""
from __future__ import annotations

from namepredict.layer5.acyl_halide_names import acyl_halide_names
from namepredict.layer5.boronic_names import boronic_names
from namepredict.layer5.carbamate_names import carbamate_names
from namepredict.layer5.carbonate_names import carbonate_names
from namepredict.layer5.diester_names import diester_names
from namepredict.layer5.hydrazine_names import hydrazine_names
from namepredict.layer5.isocyanate_names import iso_kind_names
from namepredict.layer5.sulfonamide_names import sulfonamide_names
from namepredict.layer5.sulfonate_names import sulfonate_names
from namepredict.layer5.sulfoxide_names import sulfoxide_names
from namepredict.layer5.urea_names import urea_names


def _table(n: int, numbered: dict) -> dict:
    return {
        "carbamate": lambda: carbamate_names(numbered),
        "carbonate": lambda: carbonate_names(numbered),
        "urea": lambda: urea_names(numbered),
        "hydrazine": lambda: hydrazine_names(numbered),
        "diester": lambda: diester_names(n, numbered),
        "sulfonamide": lambda: sulfonamide_names(numbered),
        "sulfonate": lambda: sulfonate_names(numbered),
        "boronic": lambda: boronic_names(numbered),
    }


def _by_kind(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    fn = _table(n, numbered).get(kind)
    return fn() if fn else None


def special_fg_names(kind: str, n: int, numbered: dict) -> tuple[str, str] | None:
    top = _by_kind(kind, n, numbered)
    if top is not None:
        return top
    if kind in ("isocyanate", "isothiocyanate"):
        return iso_kind_names(kind, n)
    if kind in ("acyl_chloride", "acyl_bromide"):
        return acyl_halide_names(kind, n, numbered)
    return sulfoxide_names(numbered) if kind == "sulfoxide" else None
