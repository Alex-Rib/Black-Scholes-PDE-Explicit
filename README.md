# Explicit Finite-Difference Scheme for the Black-Scholes PDE

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Finance](https://img.shields.io/badge/Finance-Derivatives-green)
![Tests](https://img.shields.io/badge/Tests-pytest-purple)
![Lint](https://img.shields.io/badge/Lint-ruff-orange)
![Status](https://img.shields.io/badge/Status-Educational-orange)

## 📊 Description

Implementation of an **explicit finite-difference scheme** to solve the Black-Scholes partial differential equation (PDE) in log-price and price a European call.

## 🎯 Objectives

- Solve the Black-Scholes PDE numerically by finite differences.
- Study the convergence towards the analytical Black-Scholes price.
- Analyse the stability condition of the explicit scheme.
- Study the convergence order of the error as a function of the space step $\delta$.

## 📐 Mathematical Model

### Transformed Black-Scholes PDE

Setting $x = \ln(S)$, the option price $u(t, x)$ satisfies the PDE:

$$\frac{\partial u}{\partial t} + \left(r - \frac{1}{2}\sigma^2\right)\frac{\partial u}{\partial x} + \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2} - r u = 0$$

### Terminal and boundary conditions (Dirichlet)

- **Terminal condition**:
  $$u(T, x) = \max\left(e^x - K, 0\right)$$
- **Lower boundary** ($x \to x_{min}$): $u = 0$
- **Upper boundary** ($x \to x_{max}$): $u = e^x - K e^{-r(T-t)}$

The space grid is centred on $\ln(S_0)$ and its half-width is set from a Gaussian quantile of the terminal log-price distribution, so the grid node $M/2$ corresponds exactly to $S_0$ ($M$ even).

## 🔧 Numerical Method

### Discretisation

- **Time step**: $h = T/N$
- **Space step**: $\delta = (x_{max} - x_{min})/M$

### Explicit scheme

The value $u_i^{n-1}$ (previous time step) is computed explicitly from the values at time $n$:

$$u_i^{n-1} = A \, u_i^n + B \, u_{i+1}^n + C \, u_{i-1}^n$$

with coefficients:

$$A = 1 - rh - \frac{\sigma^2 h}{\delta^2}$$

$$B = \frac{h\sigma^2}{2\delta^2} + \frac{h(r - 0.5\sigma^2)}{2\delta}$$

$$C = \frac{h\sigma^2}{2\delta^2} - \frac{h(r - 0.5\sigma^2)}{2\delta}$$

The interior update is fully vectorised with NumPy slicing.

### Stability condition

The explicit scheme is only conditionally stable; it requires:

$$h < \frac{\delta^2}{\sigma^2}$$

The condition is checked automatically: `ExplicitScheme.price()` raises a `ValueError` on an unstable grid, and the experiment script skips those configurations.

## 📁 Project Structure

```
.
├── Black_scholes.py     # Closed-form Black-Scholes call price
├── Explicit_scheme.py   # OptionParams, GridConfig, ExplicitScheme
├── Run_experiments.py   # Stability check + convergence study
└── tests/
    └── Test_explicit_scheme.py
```

## 📈 Experiments

`Run_experiments.py` performs, for $M$ from 50 to 2000 at fixed $N = 180{,}000$:
1. **Stability check** for each grid (unstable grids are skipped).
2. **Price convergence plot** against the closed-form Black-Scholes price.
3. **Log-log error plot** showing the decay of the absolute error with $M$.

## 🚀 Usage

```bash
python Run_experiments.py
```

Minimal pricing example:

```python
from Explicit_scheme import ExplicitScheme, GridConfig, OptionParams

params = OptionParams(S0=100, K=100, T=1.0, r=0.05, sigma=0.2)
grid = GridConfig(M=800, N=180_000)
price = ExplicitScheme(params, grid).price()
```

## ✅ Tests


Six unit tests cover: stability detection, `ValueError` on unstable grids, convergence to the closed-form price, error reduction under grid refinement, no-arbitrage lower bound, and grid validation.

## 👨‍💻 Author

Alexandre R. - Université Paris Cité