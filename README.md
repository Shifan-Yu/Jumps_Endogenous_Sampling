# Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times

This repository contains the MATLAB code and simulation examples for **“Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.”**
You can find the latest draft and supplemental materials on the author’s website.

---

## What’s inside

* **Test statistic (`testLLNNY.m`)**:
Constructs the test statistic from the implied barriers $\widehat m$ and $\widehat m_{\epsilon}$ by inverting the functions $h_{2}(m)$ and $h_{2,\epsilon}(m)$.

* **Precomputed tables**:

  * `h_vec.mat`: two-column table $[m\; h_{2}(m)]$
  * `h_eps_vec.mat`: wide table across $\epsilon$ from 0.01 to 1.00 with $[m\; h_{2,\epsilon}(m)]$
  * `avar_r.mat`: wide table across $\epsilon$ with $[m\; V_{\epsilon}(m)]$

* **Simulators**:
  Heston model for the efficient price with tick-level observations. Jumps can be added with selected intensities; an option with market microstructure noise (autocorrelated Gaussian-t mixture noise and rounding errors) is included.

* **Monte Carlo examples**:
  Size and size-adjusted power across barrier widths $c$ and censoring parameters $\epsilon$.

---

## Package layout

* **`numerical_h_hbar/` — Reproduce the function tables**

  * `h_simulate.m`: compute $h_{2}$, $h_{2,\epsilon}$, first-order derivatives, and $V_{\epsilon}$ on an $m$ grid
  * `h_first_derivative.m`: local-linear slope estimator of $h'(m)$
  * `ret_delta.m`: Price duration sampled (barrier-hitting) returns $r^{(c)}$
  * `h_vec.mat`, `h_eps_vec.mat`, `avar_r.mat`: saved tables used by the test

* **`testToolbox/` — Test statistic and helper functions**

  * `testLLNNY.m`: builds the test statistic from implied barriers and $V_{\epsilon}$
  * `invFunc.m`: robust monotone inverse on tabulated $[m, Y(m)]$
  * `linearInt.m`: clamped `interp1` wrapper for table lookup
  * `ret_delta.m`: price duration (barrier-hitting) sampling
  * `finddata.m`: slice $\epsilon$-specific columns to two-column tables
  * `wb_preaveraging.m`: wild-bootstrapped preaveraging

* **`MCExamples/` — Monte Carlo Examples**

  * `simulatePrices/`

    * `simPriceEfficient.m`: Heston + Jumps
    * `simPriceNoise_autoGau_t.m`: Heston + Jumps + Microstructure Noise
    * `MA_noise.m`, `rounding.m`, `randLaplace.m`: Helper functions for microstructure noise simulation
  * `mainNoNoise.m`: empirical size & power without microstructure noise
  * `mainWithNoise.m`: empirical size & power with noise (with pre-averaging)

---

## Reproducing the tables $h_{2}$, $h_{2,\epsilon}$, $V_{\epsilon}$

1. Simulate a long Gaussian random walk and, for each grid value $m$, compute $h_{2}(m)$ and $h_{2,\epsilon}(m)$ (censoring at $m(1+\epsilon)$).
2. Estimate $h'_{2}(m)$ and $h'_{2,\epsilon}(m)$ with local-linear slopes (`h_first_derivative.m`).
3. Assemble the variance function $V_{\epsilon}(m)$ via the delta-method expressions.
4. Save `h_vec.mat`, `h_eps_vec.mat`, `avar_r.mat`.

> Assumptions for inversion: the maps are **monotone decreasing** and smooth on the grid. `invFunc.m` combines a tail quadratic (on $\mu_2$) and a local inverse quadratic (in $1/m$) with a monotone interpolation fallback.

---

## Monte Carlo usage

* **Null:** $X = X^{c}$ (diffusion; optionally with noise in levels).
* **Alternatives:** $X = X^{c} + X^{d}$ (compound Laplace jump levels).
* **Barrier choice per path:** $c_i = K \sqrt{\mathrm{Var}(\Delta X_i)}$. With noise, you may pre-average **only to compute $c_i$**.
* **Power reporting:** condition on paths with at least one jump ($N>0$).
* **Size-adjusted power:** compare to the **empirical** 95th percentile from the null **for each $K$**.

---

## Helper functions (brief)

* `testLLNNY(X, c, eps, H2_tab, H2eps_tab, Veps_tab)` → `[Z_eps, m_hat_eps, m_hat, N_c]`
  Builds $M_c$ and $M_{c,\epsilon}$, inverts to $\widehat m,\widehat m_{\epsilon}$, looks up $V_{\epsilon}(\widehat m_{\epsilon})$, and returns $Z_{\epsilon}$.

* `finddata(eps, h_eps_vec, avar_r)`
  From wide matrices across $\epsilon$, extract two-column tables $[m, h_{2,\epsilon}(m)]$ and $[m, V_{\epsilon}(m)]$.

* `invFunc(Y_trial, [m, Y(m)])`
  Robust monotone inverse of a **decreasing** tabulated map.

* `linearInt(x, [xgrid, ygrid], method)`
  Clamped `interp1` lookup (default `'linear'`).

* `ret_delta(price, m)`
  Barrier-hitting (first-passage) returns $r^{(m)}$.

* `wb_preaveraging(x, kn, 'flip',bool,'perm',bool)`
  Pre-averaged path and, if requested, wild-bootstrapped pseudo path (sign-flip + permutation).
  `pseudo.m` returns the pre-averaged levels only.
