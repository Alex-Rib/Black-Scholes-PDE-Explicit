"""Closed-form Black-Scholes pricing for a European call."""

import numpy as np
from scipy.stats import norm


def call_price(S: float, K: float, T: float, r: float, sigma: float) -> float:
    """Return the Black-Scholes price of a European call.

    Parameters
    ----------
    S : float
        Spot price.
    K : float
        Strike price.
    T : float
        Time to maturity in years.
    r : float
        Risk-free rate.
    sigma : float
        Volatility.
    """
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    return S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
