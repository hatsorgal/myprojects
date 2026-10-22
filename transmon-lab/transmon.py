"""Single-junction transmon calculations, with energies expressed as E/h in GHz."""

from dataclasses import dataclass
import math

import numpy as np
from scipy.constants import elementary_charge, h
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import least_squares


def _positive(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive.")


def charging_energy_ghz(capacitance_ff: float) -> float:
    """EC = e²/(2 CΣ), where CΣ is total island capacitance, not just a junction."""
    _positive(capacitance_ff, "Total capacitance")
    return elementary_charge**2 / (2 * capacitance_ff * 1e-15 * h) / 1e9


def capacitance_ff(ec_ghz: float) -> float:
    _positive(ec_ghz, "Charging energy")
    return elementary_charge**2 / (2 * ec_ghz * 1e9 * h) / 1e-15


def spectrum(
    ej_ghz: float,
    ec_ghz: float,
    offset_charge: float = 0.0,
    charge_cutoff: int = 25,
    levels: int = 4,
) -> np.ndarray:
    """Diagonalize 4EC(n-ng)² - EJ cos(phi) on n=-cutoff,...,+cutoff.

    ng is measured in Cooper pairs (2e). Results are GHz relative to the ground
    state; convergence with charge_cutoff remains the caller's responsibility.
    """
    _positive(ej_ghz, "Josephson energy")
    _positive(ec_ghz, "Charging energy")
    if not math.isfinite(offset_charge):
        raise ValueError("Offset charge must be finite.")
    if not isinstance(charge_cutoff, int) or charge_cutoff < 4:
        raise ValueError("Charge cutoff must be an integer of at least 4.")
    if not isinstance(levels, int) or not 2 <= levels <= 2 * charge_cutoff + 1:
        raise ValueError("Requested levels must fit inside the charge basis.")
    # Integer shifts of ng relabel the charge states. Wrapping also prevents a
    # large ng from putting the low-energy states outside the truncated basis.
    ng = (offset_charge + 0.5) % 1.0 - 0.5
    charges = np.arange(-charge_cutoff, charge_cutoff + 1, dtype=float)
    diagonal = 4.0 * ec_ghz * (charges - ng) ** 2
    off_diagonal = np.full(2 * charge_cutoff, -ej_ghz / 2.0)
    eigenvalues = eigh_tridiagonal(
        diagonal, off_diagonal, eigvals_only=True,
        select="i", select_range=(0, levels - 1),
    )
    return eigenvalues - eigenvalues[0]


@dataclass(frozen=True)
class Observables:
    f01_ghz: float
    f12_ghz: float
    alpha_mhz: float


def observables(ej_ghz: float, ec_ghz: float, offset_charge: float = 0.0,
                charge_cutoff: int = 25) -> Observables:
    energies = spectrum(ej_ghz, ec_ghz, offset_charge, charge_cutoff)
    f01 = float(energies[1] - energies[0])
    f12 = float(energies[2] - energies[1])
    return Observables(f01, f12, (f12 - f01) * 1000)


def approximation(ej_ghz: float, ec_ghz: float) -> Observables:
    _positive(ej_ghz, "Josephson energy")
    _positive(ec_ghz, "Charging energy")
    f01 = math.sqrt(8 * ej_ghz * ec_ghz) - ec_ghz
    return Observables(f01, f01 - ec_ghz, -ec_ghz * 1000)


def convergence_error_khz(ej_ghz: float, ec_ghz: float,
                          offset_charge: float = 0.0,
                          charge_cutoff: int = 25) -> float:
    """Largest transition-frequency change when adding 20 charge-basis states."""
    small = np.diff(spectrum(ej_ghz, ec_ghz, offset_charge, charge_cutoff=charge_cutoff))
    large = np.diff(spectrum(ej_ghz, ec_ghz, offset_charge, charge_cutoff=charge_cutoff + 10))
    return float(np.max(np.abs(small - large)) * 1e6)


def charge_dispersion_khz(ej_ghz: float, ec_ghz: float) -> float:
    """Difference in f01 at ng=0 and ng=0.5 (not a decoherence-rate estimate)."""
    return abs(observables(ej_ghz, ec_ghz, 0).f01_ghz
               - observables(ej_ghz, ec_ghz, 0.5).f01_ghz) * 1e6


@dataclass(frozen=True)
class Design:
    ej_ghz: float
    ec_ghz: float
    capacitance_ff: float
    predicted: Observables
    simulator_calls: int
    frequency_error_khz: float
    alpha_error_khz: float
    basis_error_khz: float


def design_transmon(target_f01_ghz: float, target_alpha_mhz: float) -> Design:
    """Fit EJ and EC to f01 and alpha at ng=0; this is numerical optimization."""
    if not math.isfinite(target_f01_ghz) or not 2 <= target_f01_ghz <= 10:
        raise ValueError("Target frequency must be between 2 and 10 GHz.")
    if not math.isfinite(target_alpha_mhz) or not -500 <= target_alpha_mhz <= -50:
        raise ValueError("Target anharmonicity must be between -500 and -50 MHz.")
    ec_guess = abs(target_alpha_mhz) / 1000
    ej_guess = (target_f01_ghz + ec_guess) ** 2 / (8 * ec_guess)
    calls = 0

    def residual(log_energies: np.ndarray) -> np.ndarray:
        nonlocal calls
        calls += 1
        ej, ec = np.exp(log_energies)
        predicted = observables(float(ej), float(ec), charge_cutoff=35)
        return np.array([
            (predicted.f01_ghz - target_f01_ghz) / target_f01_ghz,
            (predicted.alpha_mhz - target_alpha_mhz) / abs(target_alpha_mhz),
        ])

    solution = least_squares(
        residual, np.log([ej_guess, ec_guess]),
        bounds=(np.log([0.01, 0.01]), np.log([500.0, 5.0])),
        xtol=1e-11, ftol=1e-11, gtol=1e-11, max_nfev=100,
    )
    ej, ec = (float(x) for x in np.exp(solution.x))
    predicted = observables(ej, ec, charge_cutoff=35)
    frequency_error = (predicted.f01_ghz - target_f01_ghz) * 1e6
    alpha_error = (predicted.alpha_mhz - target_alpha_mhz) * 1000
    basis_error = convergence_error_khz(ej, ec, charge_cutoff=35)
    if (not solution.success or abs(frequency_error) > 10
            or abs(alpha_error) > 10 or basis_error > 1):
        raise ValueError("The solver could not verify this design. Try another target.")
    return Design(ej, ec, capacitance_ff(ec), predicted, calls,
                  frequency_error, alpha_error, basis_error)
