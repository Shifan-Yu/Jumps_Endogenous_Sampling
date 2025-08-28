# Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times

This repository contains the MATLAB code and simulation examples for **“Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.”**
You can find the latest draft and supplemental materials on the author’s website.

---

## What’s inside

* **Test statistic (`testLLNNY.m`)**:
Constructs the test statistic from the empirical quantities $M_{c}$ and $M_{c,\epsilon}$ by inverting the tabulated maps $h_{2}(m)$ and $h_{2,\epsilon}(m)$ to obtain the implied barriers $\widehat m$ and $\widehat m_{\epsilon}$. The statistic uses the precomputed variance function $V_{\epsilon}(m)$ on the same barrier grid.

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
  * `ret_delta.m`: price duration sampled (barrier-hitting) returns $r^{(c)}$
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

    * `simPriceEfficient.m`: Heston + jumps
    * `simPriceNoise_autoGau_t.m`: Heston + jumps + microstructure noise
    * `MA_noise.m`, `rounding.m`, `randLaplace.m`: Helper functions for microstructure noise simulation
  * `mainNoNoise.m`: empirical size & power without microstructure noise
  * `mainWithNoise.m`: empirical size & power with noise (with pre-averaging)

---

## Workflow

1. **Choose a working barrier $c$.**
   In finite samples we set $c$ as $K$ times the standard deviation of tick-by-tick returns. Under microstructure noise, we construct the sequence of pseudo-observations with selected pre-averaging windows with `wb_preaveraging`, then set $c$ based on pre-averaged returns, e.g., `X_pa = wb_preaveraging(X, round(0.5*sqrt(n))); c = K*sqrt(var(diff(X_pa),1));`

2. **Compute $M_c$ and $M_{c,\epsilon}$.**
   Use `ret_delta(X,c)` to extract $r^{(c)}$, then form the two empirical moments in `testLLNNY`.

3. **Invert $h_2$ and $h_{2,\epsilon}$.**
   With `H2_tab = [m,h_2(m)]` and `H2eps_tab = [m,h_{2,\epsilon}(m)]` (`finddata(eps, h_eps_vec, avar_r)` returns two-column tables for the chosen $\epsilon$), `invFunc` returns $\widehat m$ and $\widehat m_{\epsilon}$.

4. **Evaluate $V_{\epsilon}$ at $\widehat m_{\epsilon}$.**
   `linearInt` performs the lookup to get $V_{\epsilon}(\widehat m_{\epsilon})$.

5. **Construct $T_{c,\epsilon}$ and test.**
   Under the null, $T_{c,\epsilon}$ is asymptotically standard normal. For Monte Carlo size-adjusted power, compare to the empirical 95th percentile from the null for each $K$.


---

## Reproducing the tables $h_{2}$, $h_{2,\epsilon}$, $V_{\epsilon}$

* **`h_simulate.m`** simulates a long Gaussian random walk (of length $10^9$) and sweeps an $m$-grid. For each $m$, it computes
  $\mu_2(m)=\mathbb{E}[(r^{(m)})^2]$ and
  $\mu_{2,\epsilon}(m)=\mathbb{E}[\{(r^{(m)}\wedge m(1+\epsilon))^2\}]$,
  then tabulates $h_{2}(m)=\mu_2(m)/m^2$ and $h_{2,\epsilon}(m)=\mu_{2,\epsilon}(m)/m^2$.

* **Derivatives.** `h_first_derivative.m` estimates $h_2'(m)$ and $h_{2,\epsilon}'(m)$ by local-linear regression around each grid point of $m$, which feed the delta-method variance.

* **Variance function.** The script also computes the needed variance and covariance terms for $(M_c,M_{c,\epsilon})$ (denoted $v$, $v_{\epsilon}$, $c_{\epsilon}$ in the code) and assembles $V_{\epsilon}(m)$. The saved file `avar_r.mat` stores $V_{\epsilon}(m)$ across the same grid.

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
