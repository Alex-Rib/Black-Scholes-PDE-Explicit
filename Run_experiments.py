"""Convergence and stability experiments for the explicit scheme.

Run with: python run_experiments.py
"""

import matplotlib.pyplot as plt
import numpy as np

from Black_scholes import call_price
from Explicit_scheme import ExplicitScheme, GridConfig, OptionParams


def stability_and_convergence(
    params: OptionParams, m_values: np.ndarray, n_steps: int
) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute prices and errors over space resolutions, skipping unstable grids."""
    prices = np.full(len(m_values), np.nan)
    errors = np.full(len(m_values), np.nan)
    bs_price = call_price(params.S0, params.K, params.T, params.r, params.sigma)

    for i, m in enumerate(m_values):
        scheme = ExplicitScheme(params, GridConfig(M=int(m), N=n_steps))
        if not scheme.is_stable():
            print(f"Stability condition violated for M = {m}, skipped.")
            continue
        print(f"Stability condition satisfied for M = {m}.")
        prices[i] = scheme.price()
        errors[i] = abs(prices[i] - bs_price)

    return prices, errors, bs_price


def plot_convergence(m_values: np.ndarray, prices: np.ndarray, bs_price: float) -> None:
    """Plot the PDE price against the closed-form price."""
    plt.figure(figsize=(12, 6))
    plt.plot(m_values, prices, "o-", label="Explicit scheme price")
    plt.axhline(y=bs_price, color="red", linestyle="--", label="Closed-form Black-Scholes")
    plt.xlabel("Number of space points M")
    plt.ylabel("Price")
    plt.title("Convergence of the explicit scheme towards the analytical price")
    plt.legend()
    plt.grid()
    plt.gca().ticklabel_format(style="plain", axis="y", useOffset=False)
    plt.tight_layout()
    plt.show()


def plot_error(m_values: np.ndarray, errors: np.ndarray) -> None:
    """Plot the absolute error in log-log scale."""
    plt.figure(figsize=(12, 6))
    plt.plot(m_values, errors, "o-")
    plt.xlabel("Number of space points M")
    plt.ylabel("Absolute error")
    plt.title("Error of the explicit scheme as a function of M (log-log)")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.show()


def main() -> None:
    """Run the stability check and convergence study, then plot the results."""
    params = OptionParams()
    m_values = np.array([50, 100, 800, 1000, 1200, 1400, 1600, 1800, 2000])

    prices, errors, bs_price = stability_and_convergence(params, m_values, n_steps=180_000)
    plot_convergence(m_values, prices, bs_price)
    plot_error(m_values, errors)


if __name__ == "__main__":
    main()
