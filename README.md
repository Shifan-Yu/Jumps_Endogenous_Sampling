# 📄 Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times

This repository contains the MATLAB code and simulation examples for **“Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.”**
You can find the latest draft and supplemental materials on [the author’s website](https://www.shifanyu.com/).

---

## 🔍 What’s inside

* **Test statistic (`testLLNNY.m`)**:
Constructs the test statistic from the empirical quantities $M_{c}$ and $\overline M_{c,\epsilon}$ by inverting the tabulated maps $h_{2}(m)$ and $h_{2,\epsilon}(m)$ to obtain the implied barriers $\widehat m$ and $\widehat m_{\epsilon}$. The statistic then uses the precomputed variance function $V_{\epsilon}(m)$ on the same $m$-grid.

* **Precomputed tables**:

  * `h_vec.mat`: two-column table $[m\; h_{2}(m)]$
  * `h_eps_vec.mat`: wide table across $\epsilon \in [0.01,1.00]$ with $[m\; h_{2,\epsilon}(m)]$
  * `avar_r.mat`: wide table across $\epsilon$ with $[m\; V_{\epsilon}(m)]$

* **Simulators**:
  Heston model for the efficient price observed at tick level. Jumps can be added with user-selected intensities and scales. An option with market microstructure noise (autocorrelated Gaussian–$t$ mixture and rounding) is included.

* **Monte Carlo examples**:
  Empirical size and size-adjusted power across working barrier widths $c$ and censoring parameters $\epsilon$.

---

## 📦 Package layout

* **`numerical_h_hbar/` — Reproduce the function tables**

  * `h_simulate.m`: compute $h_{2}$, $h_{2,\epsilon}$, first-order derivatives, and $V_{\epsilon}$ on an $m$ grid
  * `h_first_derivative.m`: local-linear slope estimator of $h'(m)$
  * `ret_delta.m`: price duration (barrier-hitting) sampled returns $r^{(c)}$
  * `h_vec.mat`, `h_eps_vec.mat`, `avar_r.mat`: saved tables used by the test

* **`testToolbox/` — Test statistic and helper functions**

  * `testLLNNY.m`: builds $T_{c,\epsilon}$ from implied barriers and $V_{\epsilon}$
  * `invFunc.m`: robust monotone inverse of tabulated $[m, Y(m)]$
  * `linearInt.m`: clamped `interp1` wrapper for table lookup
  * `ret_delta.m`: price duration (barrier-hitting) sampling
  * `finddata.m`: slice $\epsilon$-specific columns to two-column tables
  * `wb_preaveraging.m`: wild-bootstrapped preaveraging

* **`MCExamples/` — Monte Carlo Examples**

  * `simulatePrices/`
  
    * `simPriceEfficient.m`: Heston + jumps
    * `simPriceNoise_autoGau_t.m`: Heston + jumps + microstructure noise
    * `MA_noise.m`, `rounding.m`, `randLaplace.m`: helpers for noise and jumps
  * `mainNoNoise.m`: empirical size & power without microstructure noise
  * `mainWithNoise.m`: empirical size & power with noise (with pre-averaging)

---

## 🛠️ Step-by-step: how `testLLNNY.m` constructs the statistic

1. **Choose a working barrier $c$.**
   In finite samples, set $c$ as $K$ times the standard deviation of tick-by-tick returns. Under microstructure noise, construct the sequence of pseudo-observations with selected pre-averaging windows with `wb_preaveraging`, then set $c$ from the pre-averaged returns:
   ```matlab
   X_pa = wb_preaveraging(X, round(0.5*sqrt(n)));
   c    = K * sqrt(var(diff(X_pa), 1));
   ```

3. **Compute $M_c$ and $\overline M_{c,\epsilon}$.**
   Use `ret_delta(X,c)` to extract $r^{(c)}$, then form the two empirical moments in `testLLNNY`.

4. **Invert $h_2$ and $h_{2,\epsilon}$.**
   With `H2_tab = [m,h_2(m)]` and `H2eps_tab = [m,h_{2,\epsilon}(m)]` (obtained via `finddata(eps, h_eps_vec, avar_r)` for the chosen $\epsilon$), use `invFunc` to get $\widehat m$ and $\widehat m_{\epsilon}$.

5. **Evaluate $V_{\epsilon}$ at $\widehat m_{\epsilon}$.**
   Use `linearInt` to look up $V_{\epsilon}(\widehat m_{\epsilon})$.

6. **Construct $T_{c,\epsilon}$ and test.**
   Under the null, $T_{c,\epsilon}$ is asymptotically standard normal. For Monte Carlo size-adjusted power, compare to the empirical 95th percentile from the null for each $K$.

---

## 🗄️ Reproducing the tables $h_{2}$, $h_{2,\epsilon}$, $V_{\epsilon}$

* **`h_simulate.m`** simulates a long standard-Gaussian random walk (length $10^9$) and sweeps an integer $m$-grid (e.g., $1$ to $90$). For each $m$ it computes
  $\mu_2(m)=\mathbb{E}[(r^{(m)})^2]$ and
  $\mu_{2,\epsilon}(m)=\mathbb{E}[\{(r^{(m)}\wedge m(1+\epsilon))^2\}]$ for $\epsilon \in [0.01,1.00]$,
  then tabulates $h_{2}(m)=\mu_2(m)/m^2$ and $h_{2,\epsilon}(m)=\mu_{2,\epsilon}(m)/m^2$.

* **Derivatives.** `h_first_derivative.m` estimates $h_2'(m)$ and $h_{2,\epsilon}'(m)$ by local-linear regression around each grid point of $m$, which feed the delta-method variance.

* **Variance function.** The script computes the variance and covariance terms for $(M_c, \overline M_{c,\epsilon})$ (denoted $v$, $v_{\epsilon}$, $c_{\epsilon}$ in code) and assembles $V_{\epsilon}(m)$. The file `avar_r.mat` stores $V_{\epsilon}(m)$ on the same $m$-grid.

---

## 📈 Monte Carlo simulations

* **Null:** $X = X^{c}$. Diffusion; optionally with noise in price observations.
* **Alternatives:** $X = X^{c} + X^{d}$. $X^{d}$ is a compound Poisson process with rate $\lambda$ and double-exponential (Laplace) jump sizes (location 0, scale $\beta$).
* **Barrier choice per path:** $c = K \sqrt{\mathrm{Var}(\Delta_i^n X)}$. With noise, use `wb_preaveraging` to build pseudo-observations and set $c$ from pre-averaged returns.
* **Power reporting:** condition on paths with at least one jump ($N>0$).
* **Size-adjusted power:** compare to the empirical 95th percentile from the null for each $K$.

---

## 📚 Reference

Li, Q., Li, Y., Nolte, I., Nolte, S., and Yu, S. (2025). Testing for jumps in a discretely observed price process with endogenous sampling times. Working Paper.
