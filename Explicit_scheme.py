"""Explicit finite-difference scheme for the Black-Scholes PDE in log-price."""

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm


@dataclass(frozen=True)
class OptionParams:
    """European call option parameters."""

    S0: float = 100.0
    K: float = 100.0
    T: float = 1.0
    r: float = 0.05
    sigma: float = 0.2


@dataclass(frozen=True)
class GridConfig:
    """Space-time grid for the explicit scheme.

    The space grid is in log-price, centred on log(S0) and sized from a
    Gaussian quantile of the terminal log-price distribution, so that the
    index M // 2 corresponds exactly to S0 (M must be even).
    """

    M: int = 800
    N: int = 180_000
    confidence: float = 0.99

    def __post_init__(self) -> None:
        if self.M % 2 != 0:
            raise ValueError("M must be even so that the grid is centred on S0.")

    def space_grid(self, params: OptionParams) -> np.ndarray:
        """Return the log-price grid (M + 1 points centred on log(S0))."""
        z = norm.ppf(self.confidence)
        half_width = abs(params.r - 0.5 * params.sigma**2) * params.T + z * params.sigma * np.sqrt(
            params.T
        )
        x_min = np.log(params.S0) - half_width
        x_max = np.log(params.S0) + half_width
        return np.linspace(x_min, x_max, self.M + 1)


class ExplicitScheme:
    """Backward-in-time explicit scheme, fully vectorised in space.

    Each time step updates
    u_i^{n-1} = A u_i^n + B u_{i+1}^n + C u_{i-1}^n,
    with Dirichlet boundary conditions. The scheme is conditionally stable:
    it requires h < delta^2 / sigma^2.
    """

    def __init__(self, params: OptionParams, grid: GridConfig) -> None:
        self.params = params
        self.grid = grid

    def is_stable(self) -> bool:
        """Return True if the stability condition h < delta^2 / sigma^2 holds."""
        x = self.grid.space_grid(self.params)
        dx = x[1] - x[0]
        h = self.params.T / self.grid.N
        return h < dx**2 / self.params.sigma**2

    def price(self) -> float:
        """Return the call price at (t = 0, S = S0).

        Raises
        ------
        ValueError
            If the stability condition is violated.
        """
        if not self.is_stable():
            raise ValueError(
                "Stability condition h < delta^2 / sigma^2 violated: increase N or decrease M."
            )

        p, g = self.params, self.grid
        x = g.space_grid(p)
        dx = x[1] - x[0]
        h = p.T / g.N
        m = g.M

        u = np.maximum(np.exp(x) - p.K, 0.0)  # terminal payoff

        a = 1 - p.r * h - p.sigma**2 * h / dx**2
        b = h * p.sigma**2 / (2 * dx**2) + h * (p.r - 0.5 * p.sigma**2) / (2 * dx)
        c = h * p.sigma**2 / (2 * dx**2) - h * (p.r - 0.5 * p.sigma**2) / (2 * dx)

        u_new = np.zeros(m + 1)
        for i in range(g.N, 0, -1):
            t_new = (i - 1) * h

            # vectorised interior update
            u_new[1:-1] = a * u[1:-1] + b * u[2:] + c * u[:-2]

            # Dirichlet boundary conditions
            u_new[0] = 0.0
            u_new[m] = np.exp(x[-1]) - p.K * np.exp(-p.r * (p.T - t_new))

            u[:] = u_new

        return float(u[m // 2])
