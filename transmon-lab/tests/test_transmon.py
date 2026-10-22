import math

import numpy as np
import pytest
from scipy.constants import elementary_charge, h

from transmon import (
    approximation, capacitance_ff, charging_energy_ghz, convergence_error_khz,
    design_transmon, observables, spectrum,
)


def test_capacitance_energy_units():
    expected = elementary_charge**2 / (2 * 80e-15 * h) / 1e9
    assert charging_energy_ghz(80) == pytest.approx(expected)
    assert 0.24 < expected < 0.25
    assert capacitance_ff(expected) == pytest.approx(80)


def test_spectrum_matches_independent_dense_hamiltonian():
    ej, ec, ng, cutoff = 17.0, 0.27, 0.23, 25
    n = np.arange(-cutoff, cutoff + 1)
    matrix = np.diag(4 * ec * (n - ng) ** 2)
    matrix += np.diag(np.full(2 * cutoff, -ej / 2), 1)
    matrix += np.diag(np.full(2 * cutoff, -ej / 2), -1)
    expected = np.linalg.eigvalsh(matrix)[:4]
    expected -= expected[0]
    assert spectrum(ej, ec, ng) == pytest.approx(expected, abs=1e-10)


def test_charge_periodicity_and_reflection():
    reference = spectrum(15, 0.3, 0.2)
    assert spectrum(15, 0.3, 1.2) == pytest.approx(reference, abs=1e-10)
    assert spectrum(15, 0.3, -0.2) == pytest.approx(reference, abs=1e-10)


def test_transmon_limit_and_basis_convergence():
    calculated = observables(20, 0.24)
    estimated = approximation(20, 0.24)
    assert calculated.f01_ghz == pytest.approx(estimated.f01_ghz, abs=0.02)
    assert -280 < calculated.alpha_mhz < -240
    assert convergence_error_khz(60, charging_energy_ghz(200), 0.5) < 1


def test_high_ratio_requires_a_larger_charge_basis():
    assert convergence_error_khz(255, 0.0494) > 1
    assert convergence_error_khz(255, 0.0494, charge_cutoff=35) < 0.01


@pytest.mark.parametrize("frequency, alpha", [(5, -250), (2, -500), (10, -50), (7.2, -180)])
def test_design_recovers_requested_transitions(frequency, alpha):
    result = design_transmon(frequency, alpha)
    # Recheck the candidate using a larger basis than the optimization uses.
    verified = observables(result.ej_ghz, charging_energy_ghz(result.capacitance_ff), charge_cutoff=45)
    assert verified.f01_ghz == pytest.approx(frequency, abs=1e-5)
    assert verified.alpha_mhz == pytest.approx(alpha, abs=0.01)
    assert result.basis_error_khz < 1


@pytest.mark.parametrize("bad", [0, -1, math.nan, math.inf])
def test_invalid_capacitance_rejected(bad):
    with pytest.raises(ValueError):
        charging_energy_ghz(bad)


@pytest.mark.parametrize("frequency, alpha", [(0, -250), (5, 250), (5, math.nan)])
def test_invalid_targets_rejected(frequency, alpha):
    with pytest.raises(ValueError):
        design_transmon(frequency, alpha)
