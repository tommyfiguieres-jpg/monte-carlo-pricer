# European Call Option Pricer — Monte Carlo Methods

Monte Carlo pricing of a European call option under the Black-Scholes framework, with convergence analysis, variance reduction via antithetic variates, Greeks estimation by finite differences, and systematic validation against the analytical Black-Scholes solution.

---

## Overview

Under the Black-Scholes model, the stock price $S_t$ follows a **Geometric Brownian Motion (GBM)**:

$$dS_t = r S_t \, dt + \sigma S_t \, dW_t$$

which is equivalent, by Itô's lemma, to:

$$d\log S_t = \left(r - \frac{1}{2}\sigma^2\right)dt + \sigma \, dW_t$$

where $r$ is the risk-free rate, $\sigma$ the volatility, and $W_t$ a standard Brownian motion.

**Euler-Maruyama discretisation.** The continuous SDE is approximated by discretising time into $N$ steps of size $dt = T/N$ and replacing $dW_t$ by $Z_k\sqrt{dt}$, where $Z_k \sim \mathcal{N}(0,1)$:

$$S_{k+1} = S_k \left(1 + r \, dt + \sigma Z_k \sqrt{dt}\right)$$

The exact solution of the GBM provides a useful reference:

$$S_T = S_0 \exp\!\left(\left(r - \frac{\sigma^2}{2}\right)T + \sigma\sqrt{T}\, Z\right)$$

**Monte Carlo pricing.** Instead of solving the Black-Scholes PDE, we estimate the option price by averaging discounted payoffs over $N_{mc}$ simulated paths:

$$\hat{C} = e^{-rT} \frac{1}{N_{mc}} \sum_{i=1}^{N_{mc}} \max(S_T^{(i)} - K,\, 0)$$

By the Central Limit Theorem, the standard error of this estimator scales as $\sigma_N / \sqrt{N_{mc}}$, **independently of the number of underlying assets** — which is the main advantage of Monte Carlo over PDE methods, whose computational cost grows exponentially with the number of dimensions.

---

## Features

- Euler-Maruyama simulation vs. exact GBM solution
- Monte Carlo pricer with 95% confidence intervals and convergence analysis
- Variance reduction via **antithetic variates** (pairing paths $Z$ and $-Z$)
- Greeks estimation by **central finite differences** with Common Random Numbers (CRN)
- Validation against Black-Scholes closed-form values for price and all Greeks
- Bias-variance trade-off analysis for the Gamma estimator

---

## Mathematical background

**Black-Scholes price** for a European call:

$$C = S_0 \,\mathcal{N}(d_1) - K e^{-rT} \mathcal{N}(d_2)$$

$$d_1 = \frac{\ln(S_0/K) + \left(r + \frac{\sigma^2}{2}\right)T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}$$

**Analytical Greeks:**

| Greek | Formula |
|-------|---------|
| $\Delta = \partial C / \partial S_0$ | $\mathcal{N}(d_1)$ |
| $\Gamma = \partial^2 C / \partial S_0^2$ | $\mathcal{N}'(d_1) \,/\, (S_0 \sigma \sqrt{T})$ |
| $\nu = \partial C / \partial \sigma$ | $S_0 \sqrt{T}\, \mathcal{N}'(d_1)$ |
| $\Theta = -\partial C / \partial T$ | $-S_0 \mathcal{N}'(d_1)\sigma \,/\, (2\sqrt{T}) - rKe^{-rT}\mathcal{N}(d_2)$ |

The Monte Carlo Greeks are estimated by central finite differences. For example, for Delta:

$$\Delta \approx \frac{C(S_0 + \epsilon_S) - C(S_0 - \epsilon_S)}{2\epsilon_S}$$

CRN (Common Random Numbers) — using the same $Z$ draws across all perturbed paths — ensures that differences reflect parameter changes rather than Monte Carlo noise, substantially reducing estimator variance.

---

## Results

Default parameters: $S_0 = K = 100$, $r = 3\%$, $\sigma = 20\%$, $T = 1$ year.

| Quantity | MC estimate ($N_{mc} = 10^5$) | BS analytical |
|----------|-------------------------------|---------------|
| Price    | 9.42 ± 0.09                   | 9.41          |
| $\Delta$ | 0.596 ± 0.004                 | 0.599         |
| $\Gamma$ | 0.0189 ± 0.0007               | 0.0193        |
| $\nu$    | 38.2 ± 0.5                    | 38.7          |
| $\Theta$ | −5.33 ± 0.05                  | −5.38         |

Antithetic variates reduce the standard error by a factor of ~1.3 compared to standard Monte Carlo for the same number of paths.

---

## Requirements & usage

```bash
pip install numpy matplotlib scipy jupyter
jupyter notebook mc_option_pricer_v5.ipynb
```

Tested with Python 3.10+.
