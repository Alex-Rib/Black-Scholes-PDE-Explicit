"""Unit tests for the explicit Black-Scholes PDE scheme."""

import numpy as np
import pytest

from Black_scholes import call_price
from Explicit_scheme import ExplicitScheme, GridConfig, OptionParams

PARAMS = OptionParams()


def test_stability_detection():
    """A coarse time grid with a fine space grid violates the CFL-type condition."""
    unstable = ExplicitScheme(PARAMS, GridConfig(M=800, N=1000))
    stable = ExplicitScheme(PARAMS, GridConfig(M=100, N=180_000))
    assert not unstable.is_stable()
    assert stable.is_stable()


def test_unstable_grid_raises():
    """Pricing on an unstable grid raises a ValueError."""
    scheme = ExplicitScheme(PARAMS, GridConfig(M=2000, N=1000))
    with pytest.raises(ValueError, match="Stability"):
        scheme.price()


def test_converges_to_closed_form():
    """PDE price converges to the closed-form Black-Scholes price."""
    pde_price = ExplicitScheme(PARAMS, GridConfig(M=800, N=180_000)).price()
    bs_price = call_price(PARAMS.S0, PARAMS.K, PARAMS.T, PARAMS.r, PARAMS.sigma)
    assert pde_price == pytest.approx(bs_price, abs=5e-3)


def test_error_decreases_with_refinement():
    """Refining the space grid (at stable N) reduces the error."""
    bs_price = call_price(PARAMS.S0, PARAMS.K, PARAMS.T, PARAMS.r, PARAMS.sigma)
    coarse = ExplicitScheme(PARAMS, GridConfig(M=50, N=180_000)).price()
    fine = ExplicitScheme(PARAMS, GridConfig(M=400, N=180_000)).price()
    assert abs(fine - bs_price) < abs(coarse - bs_price)


def test_price_positive_and_above_intrinsic():
    """Call price is positive and above the no-arbitrage lower bound."""
    price = ExplicitScheme(PARAMS, GridConfig(M=400, N=180_000)).price()
    lower_bound = PARAMS.S0 - PARAMS.K * np.exp(-PARAMS.r * PARAMS.T)
    assert price > 0.0
    assert price >= lower_bound


def test_odd_m_raises():
    """An odd number of space points is rejected (grid not centred on S0)."""
    with pytest.raises(ValueError):
        GridConfig(M=201, N=100)
