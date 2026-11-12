"""Numerical core for aqueous boron speciation.

The analytical concentration is expressed as mol B/L. Polyborate molecular
concentrations are therefore weighted by their number of B atoms in both the
mass balance and reported fractions of total boron.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class EquilibriumConstants:
    ka: float = 10**-9.2
    kp1: float = 10**-7.29
    kp2: float = 10**-6.77
    kp3: float = 10**-14.5
    kp4: float = 10**-16.3


CONSTANTS = EquilibriumConstants()

SPECIES = (
    "B(OH)3",
    "B(OH)4-",
    "B3O3(OH)4-",
    "B5O6(OH)4-",
    "B4O5(OH)4^2-",
    "B3O3(OH)5^2-",
)

DISPLAY_NAMES = (
    "B(OH)₃",
    "B(OH)₄⁻",
    "B₃O₃(OH)₄⁻",
    "B₅O₆(OH)₄⁻",
    "B₄O₅(OH)₄²⁻",
    "B₃O₃(OH)₅²⁻",
)

BORON_NUMBERS = np.array([1.0, 1.0, 3.0, 5.0, 4.0, 3.0])

COLORS = (
    "#64748B",
    "#14B8A6",
    "#F59E0B",
    "#8B5CF6",
    "#EF4444",
    "#3B82F6",
)


def molecular_concentrations(
    boric_acid: float | np.ndarray,
    hydrogen: float | np.ndarray,
    constants: EquilibriumConstants = CONSTANTS,
) -> np.ndarray:
    """Return molecular concentrations in the order defined by ``SPECIES``."""
    x = np.asarray(boric_acid, dtype=float)
    h = np.asarray(hydrogen, dtype=float)
    return np.array(
        [
            x,
            constants.ka * x / h,
            constants.kp1 * x**3 / h,
            constants.kp2 * x**5 / h,
            constants.kp3 * x**4 / h**2,
            constants.kp4 * x**3 / h**2,
        ]
    )


def total_boron_from_x(
    boric_acid: float,
    hydrogen: float,
    constants: EquilibriumConstants = CONSTANTS,
) -> float:
    concentrations = molecular_concentrations(boric_acid, hydrogen, constants)
    return float(BORON_NUMBERS @ concentrations)


def solve_boric_acid(
    total_boron: float,
    ph: float,
    constants: EquilibriumConstants = CONSTANTS,
) -> float:
    """Solve the unique positive mass-balance root for ``[B(OH)3]``."""
    if not np.isfinite(total_boron) or total_boron <= 0:
        raise ValueError("总硼浓度必须为大于0的有限数值。")
    if not np.isfinite(ph):
        raise ValueError("pH必须为有限数值。")

    hydrogen = 10.0**(-ph)

    def residual(x: float) -> float:
        return total_boron_from_x(x, hydrogen, constants) - total_boron

    # f(0) = -C_T and f(C_T) >= 0. All coefficients are positive, so the
    # residual is strictly increasing and this root is unique.
    lower = 0.0
    upper = float(total_boron)
    # The residual is strictly increasing. Fixed-iteration bisection reaches
    # the closest representable floating-point bracket without a SciPy runtime
    # dependency, which makes the desktop build smaller and more portable.
    for _ in range(220):
        midpoint = (lower + upper) / 2.0
        if midpoint == lower or midpoint == upper:
            break
        if residual(midpoint) > 0.0:
            upper = midpoint
        else:
            lower = midpoint
    return float((lower + upper) / 2.0)


def calculate_speciation(
    total_boron: float,
    ph_values: Iterable[float],
    constants: EquilibriumConstants = CONSTANTS,
) -> dict[str, np.ndarray]:
    """Calculate molecular concentrations and B-atom fractions over a pH grid."""
    ph = np.asarray(list(ph_values), dtype=float)
    if ph.ndim != 1 or ph.size < 1:
        raise ValueError("pH网格至少需要包含一个点。")
    if not np.all(np.isfinite(ph)):
        raise ValueError("pH网格包含非有限数值。")

    x = np.array([solve_boric_acid(total_boron, value, constants) for value in ph])
    hydrogen = 10.0 ** (-ph)
    molecular = molecular_concentrations(x, hydrogen, constants)
    fractions = BORON_NUMBERS[:, None] * molecular / total_boron
    fraction_sum = fractions.sum(axis=0)
    residual = total_boron - (BORON_NUMBERS[:, None] * molecular).sum(axis=0)

    if np.any(molecular < -1e-15) or np.any(fractions < -1e-13):
        raise RuntimeError("计算得到负浓度或负物种分数。")
    if not np.allclose(fraction_sum, 1.0, rtol=2e-11, atol=2e-11):
        raise RuntimeError("总硼质量守恒检查未通过。")

    return {
        "ph": ph,
        "hydrogen": hydrogen,
        "boric_acid": x,
        "molecular": molecular,
        "fractions": fractions,
        "fraction_sum": fraction_sum,
        "mass_balance_residual": residual,
    }


def calculate_speciation_grid(
    total_boron_values: Iterable[float],
    ph_values: Iterable[float],
    constants: EquilibriumConstants = CONSTANTS,
) -> dict[str, np.ndarray]:
    """Calculate a concentration–pH speciation grid with vectorized bisection.

    ``fractions`` has shape ``(species, concentration, pH)`` and represents
    fractions of total B atoms, identical to :func:`calculate_speciation`.
    """
    concentrations = np.asarray(list(total_boron_values), dtype=float)
    ph = np.asarray(list(ph_values), dtype=float)
    if concentrations.ndim != 1 or concentrations.size < 1:
        raise ValueError("总硼浓度网格至少需要包含一个点。")
    if ph.ndim != 1 or ph.size < 1:
        raise ValueError("pH网格至少需要包含一个点。")
    if not np.all(np.isfinite(concentrations)) or np.any(concentrations <= 0):
        raise ValueError("总硼浓度网格必须全部为大于0的有限数值。")
    if not np.all(np.isfinite(ph)):
        raise ValueError("pH网格包含非有限数值。")

    concentration_grid = concentrations[:, None]
    hydrogen_grid = 10.0 ** (-ph[None, :])
    lower = np.zeros((concentrations.size, ph.size), dtype=float)
    upper = np.broadcast_to(concentration_grid, lower.shape).copy()

    for _ in range(220):
        midpoint = (lower + upper) / 2.0
        molecular = molecular_concentrations(midpoint, hydrogen_grid, constants)
        calculated = (BORON_NUMBERS[:, None, None] * molecular).sum(axis=0)
        upper = np.where(calculated > concentration_grid, midpoint, upper)
        lower = np.where(calculated > concentration_grid, lower, midpoint)

    boric_acid = (lower + upper) / 2.0
    molecular = molecular_concentrations(boric_acid, hydrogen_grid, constants)
    fractions = BORON_NUMBERS[:, None, None] * molecular / concentration_grid[None, :, :]
    fraction_sum = fractions.sum(axis=0)
    residual = concentration_grid - (BORON_NUMBERS[:, None, None] * molecular).sum(axis=0)

    if not np.allclose(fraction_sum, 1.0, rtol=2e-11, atol=2e-11):
        raise RuntimeError("三维网格总硼质量守恒检查未通过。")
    return {
        "ph": ph,
        "total_boron": concentrations,
        "boric_acid": boric_acid,
        "molecular": molecular,
        "fractions": fractions,
        "fraction_sum": fraction_sum,
        "mass_balance_residual": residual,
    }


def integrate_window(
    total_boron: float,
    ph_min: float,
    ph_max: float,
    points: int = 401,
    constants: EquilibriumConstants = CONSTANTS,
) -> dict[str, np.ndarray | float]:
    """Integrate exact B fractions over a closed pH window.

    The generated grid includes both requested boundaries exactly. The average
    fraction is dimensionless and is usually easier to interpret than the raw
    pH-integral.
    """
    if not np.isfinite(ph_min) or not np.isfinite(ph_max) or ph_min >= ph_max:
        raise ValueError("积分窗口必须满足 pH最小值 < pH最大值。")
    if points < 3:
        raise ValueError("积分点数至少为3。")
    ph = np.linspace(ph_min, ph_max, int(points))
    result = calculate_speciation(total_boron, ph, constants)
    integrals = np.trapezoid(result["fractions"], ph, axis=1)
    width = ph_max - ph_min
    averages = integrals / width
    return {
        "ph": ph,
        "integrals": integrals,
        "average_fractions": averages,
        "integral_sum": float(integrals.sum()),
        "window_width": float(width),
    }


def cumulative_charged_fractions(fractions: np.ndarray) -> np.ndarray:
    """Return the legacy Y1–Y5 cumulative charged-B curves."""
    # Legacy order: B3O3(OH)5, B4O5(OH)4, B5O6(OH)4,
    # B3O3(OH)4, B(OH)4. Neutral B(OH)3 is excluded.
    indices = (5, 4, 3, 2, 1)
    return np.cumsum(fractions[np.array(indices)], axis=0)


def result_rows(result: dict[str, np.ndarray]):
    """Yield one export-ready dictionary per pH value."""
    for column, ph in enumerate(result["ph"]):
        row = {
            "pH": ph,
            "[H+] (M)": result["hydrogen"][column],
            "[B(OH)3] solved (M)": result["boric_acid"][column],
        }
        for index, species in enumerate(SPECIES):
            row[f"[{species}] molecular concentration (M)"] = result["molecular"][index, column]
            row[f"Fraction of total B as {species}"] = result["fractions"][index, column]
        row["Sum of B fractions"] = result["fraction_sum"][column]
        row["Mass-balance residual (M B)"] = result["mass_balance_residual"][column]
        yield row