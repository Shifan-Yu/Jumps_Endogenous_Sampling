# Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times

This repository contains the MATLAB code and simulation examples of “Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.” You can find the latest draft and supplemental materials on my website.

* **Barrier** $m$ (theoretical) and **working barrier** $c$ (chosen from data).
* **Barrier-hitting (first-passage) sampling** of returns $r^{(m)}$ / $r^{(c)}$.
* Scale-free maps $h_{2}(m)$ and $h_{2,\varepsilon}(m)$.
* **Implied barriers** $\widehat m, \widehat m_{\varepsilon}$ obtained by inverting these maps on an $m$–grid.
* **Ratio statistic** $Z_{\varepsilon}$ standardized with the paper’s $\mathrm{AVAR}_{r}(\cdot;\varepsilon)$, evaluated at $\widehat m_{\varepsilon}$.

> We deliberately avoid “RV/CRV” and “right-scaled” wording. All terms match the paper.

---

## Table of contents

* [Installation](#installation)
* [Quick start](#quick-start)
* [Paper setup (one-page recap)](#paper-setup-one-page-recap)
* [Folder structure](#folder-structure)
* [Reproducing the tables](#reproducing-the-tables)
* [Monte Carlo usage](#monte-carlo-usage)
* [Simulators](#simulators)
* [Helper functions](#helper-functions)
* [Common pitfalls](#common-pitfalls)
* [Citation](#citation)
* [License & contact](#license--contact)

---

## Installation

```matlab
cd <path-to-this-repo>
addpath(genpath(pwd));     % add all subfolders to the MATLAB path
savepath;                  % optional (persist across sessions)
```

Test your setup:

```matlab
which testLLNNY
which invFunc
which wb_preaveraging
```

You should see full paths pointing into this repository.

---

## Quick start

Minimal end-to-end example on a single path:

```matlab
% 0) Load precomputed tables
load numerical_h_hbar/h_vec.mat
load numerical_h_hbar/h_eps_vec.mat
load numerical_h_hbar/avar_r.mat

% 1) Choose ε and slice two-column tables (paper: h_{2,ε}, AVAR_r(·;ε))
eps = 0.05;
[H2eps_tab, AVAR_tab] = finddata(eps, h_eps_vec, avar_r);  % [m, h_{2,ε}(m)], [m, AVAR_r(m;ε)]
H2_tab = h_vec;                                            % [m, h_{2}(m)]

% 2) Path X (column vector); here a toy example
X = cumsum(randn(50_000,1));

% 3) Working barrier c = K * sqrt(Var(ΔX))  (practical finite-sample rule)
K = 6;
c = K * sqrt(var(diff(X), 1));    % population variance (1/N)

% 4) Ratio test (returns Z_ε, implied barriers, and #hits)
[Z, mhat_eps, mhat, N_c] = testLLNNY(X, c, eps, H2_tab, H2eps_tab, AVAR_tab);

% 5) Decision (5% right tail)
reject = (Z > 1.64485362695147);  % Φ^{-1}(0.95)
```

> **Tip.** Choose $K$ so the number of hits $N_c$ is comfortably large (hundreds+). If $N_c=0$, the function returns `NaN` (increase sample size or lower $c$).

---

## Paper setup (one-page recap)

**Barrier-hitting sampling.** Let $\tau_0=1$ and for $j\ge 1$,

$$
\tau_j=\inf\{t>\tau_{j-1}:|X_t-X_{\tau_{j-1}}|\ge m\},\qquad
r_j^{(m)}=X_{\tau_j}-X_{\tau_{j-1}}.
$$

**Scale-free maps.**

$$
h_{2}(m)=\frac{\mathbb{E}\!\left[(r^{(m)})^2\right]}{m^2},\qquad
h_{2,\varepsilon}(m)=\frac{\mathbb{E}\!\left[\min\!\big\{(r^{(m)})^2,(m(1+\varepsilon))^2\big\}\right]}{m^2}.
$$

**Empirical counterparts** at working barrier $c$ with $N_c$ hits:

$$
S_{2}=\frac{1}{N_c}\sum \frac{(r^{(c)})^2}{c^2},\qquad
S_{2,\varepsilon}=\frac{1}{N_c}\sum \frac{\min\!\{(r^{(c)})^2,(c(1+\varepsilon))^2\}}{c^2}.
$$

**Implied barriers** (inversions on the tabulated grid):

$$
\widehat m=h_{2}^{-1}(S_{2}),\qquad 
\widehat m_{\varepsilon}=h_{2,\varepsilon}^{-1}(S_{2,\varepsilon}).
$$

**Ratio statistic** (asymptotically standard normal under the null):

$$
Z_{\varepsilon}
=\frac{\widehat m_{\varepsilon}/\widehat m-1}
{\sqrt{\mathrm{AVAR}_{r}(\widehat m_{\varepsilon};\varepsilon)/n}},
$$

with $n$ the path length and $\mathrm{AVAR}_{r}(\cdot;\varepsilon)$ provided on the same $m$–grid.

---

## Folder structure

```
numerical_h_hbar/           % Precompute h_{2}, h_{2,ε}, derivatives, and AVAR_r on an m-grid
  h_simulate.m
  h_first_derivative.m
  ret_delta.m
  h_vec.mat                 % [m_grid, h_{2}(m)]
  h_eps_vec.mat             % [m_grid, h_{2,ε_k}(m)]  (wide across ε = 0.01:0.01:1)
  avar_r.mat                % [m_grid, AVAR_r(m; ε_k)] (wide across ε)

testToolbox/                % Ratio test and helpers (paper notation)
  testLLNNY.m               % builds Z_ε from implied barriers + AVAR_r
  invFunc.m                 % robust monotone inverse of a decreasing map
  linearInt.m               % clamped interp1 wrapper (table lookup)
  ret_delta.m               % barrier-hitting returns r^{(c)}
  finddata.m                % pick the ε column: wide → two-column tables
  wb_preaveraging.m         % pre-averaging; optional wild bootstrap (paper)
  pseudo.m                  % (compat) wrapper calling wb_preaveraging

MCExamples/                 % Reproducible Monte Carlo experiments
  simulatePrices/
    simPriceEfficient.m           % Heston + Laplace jumps (levels)
    simPriceNoise_autoGau_t.m     % Heston + jumps + auto Gau–t microstructure noise (levels)
    MA_noise.m, rounding.m, randLaplace.m (if used)
  mainNoNoise.m                   % size & power without microstructure noise
  mainWithNoise.m                 % with noise; optional pre-averaging for c

LLNNY_Jumps.pdf             % the paper
```

---

## Reproducing the tables

Run `numerical_h_hbar/h_simulate.m`:

1. Simulate a long Gaussian random walk and, for each grid point $m$, compute
   $h_{2}(m)$ and $h_{2,\varepsilon}(m)$ (censoring at $m(1+\varepsilon)$).
2. Estimate the derivatives $h_{2}'(m)$, $h_{2,\varepsilon}'(m)$ via local-linear slopes (`h_first_derivative.m`).
3. Combine variance components (delta method) and save $\mathrm{AVAR}_{r}(m;\varepsilon)$.
4. Save `h_vec.mat`, `h_eps_vec.mat`, `avar_r.mat`.

**Assumptions for inversion.** The maps are **monotone decreasing** and smooth on the grid. `invFunc.m` uses a tail quadratic (on $\mu_2$) and a local inverse quadratic (in $1/m$) with a monotone interpolation fallback.

---

## Monte Carlo usage

See `MCExamples/mainNoNoise.m` and `MCExamples/mainWithNoise.m`.

* **Null:** $X=X^{c}$ (diffusion; optionally with microstructure noise in levels).
* **Alternatives:** $X=X^{c}+X^{d}$ (compound Laplace jump level process).
* **Barrier choice:** per path $c_i = K\,\sqrt{\mathrm{Var}(\Delta X_i)}$.
  With noise, pre-average levels **only for barrier selection**:

  ```matlab
  kn = round(0.5*sqrt(length(X)));
  X_pa = wb_preaveraging(X, kn);     % pre-averaged levels (no flips/permutation)
  c_i  = K * sqrt(var(diff(X_pa),1));
  % Run the test on the intended series: original X (or WB pseudo if you adopt that variant).
  ```
* **Reporting power:** for jump alternatives, it is standard to **condition on paths with at least one jump** $N>0$.
* **Size-adjusted power:** compare alternatives to the empirical 95th percentile from the null **for each K**.

---

## Simulators

* `simPriceEfficient.m`
  Heston diffusion with correlated variance shocks + **compound Laplace jumps** (level process). Returns the total price $X$, continuous part $X^c$, jump level $X^d$, variance $v$, jump count $N$, and simple average $\int_0^1 v_t\,dt$.

* `simPriceNoise_autoGau_t.m`
  Same diffusion + jumps, plus **microstructure noise** added to **levels** (autocorrelated Gaussian + standardized Student-t mixture), with optional **rounding** to a price grid.

Jump arrivals are Bernoulli per grid step with probability $\lambda/\text{num}$ (≈ Poisson$(\lambda)$ over $[0,1]$); jump sizes are **Laplace** with scale $\beta$ (generated as $\text{sign}\times\mathrm{Exp}(\text{mean}=\beta)$).

---

## Helper functions

* `testLLNNY(X, c, eps, H2_tab, H2eps_tab, AVAR_tab)` → `[Z_eps, m_hat_eps, m_hat, N_c]`
  Builds $S_2$, $S_{2,\varepsilon}$; inverts to $\widehat m,\widehat m_{\varepsilon}$; looks up $\mathrm{AVAR}_{r}(\widehat m_{\varepsilon};\varepsilon)$; returns $Z_{\varepsilon}$.

* `finddata(eps, h_eps_vec, avar_r)`
  From wide matrices (across $\varepsilon$), extract two-column tables `[m, h_{2,\varepsilon}(m)]` and `[m, \mathrm{AVAR}_r(m;\varepsilon)]`.

* `invFunc(Y_trial, [m, Y(m)])`
  Robust monotone inverse of a **decreasing** tabulated map (tail quadratic + local inverse quadratic with monotone fallback).

* `linearInt(x, [xgrid, ygrid], method)`
  Clamped `interp1` lookup (`'linear'` default).

* `ret_delta(price, m)`
  Barrier-hitting returns $r^{(m)}$ (first-passage increments).

* `wb_preaveraging(x, kn, 'flip',bool,'perm',bool)`
  Pre-averaged path $x_{\text{pa}}$ and, if requested, a wild-bootstrapped pseudo path (sign-flip + permutation).
  `pseudo.m` is a compatibility wrapper that returns pre-averaged levels only.

---

## Common pitfalls

* **Mismatched $\varepsilon$.** Ensure the same $\varepsilon$ for `H2eps_tab` and `AVAR_tab` (use `finddata`).
* **Zero hits.** If $N_c=0$, the test returns `NaN`. Lower $c$ or lengthen the series.
* **Grid mismatch.** $m$–grids for $h$ and $\mathrm{AVAR}_{r}$ must match.
* **Noise design.** Pre-average for choosing $c$; run the test on the intended series (original levels, or the wild-bootstrapped pre-averaged series if you adopt that variant).
* **Reproducibility.** Set `rng(1)` before simulations; use `parpool` to speed up `parfor`.

---

## Citation

Please cite the accompanying paper `LLNNY_Jumps.pdf`. The implementation follows its barrier-hitting sampling, the maps $h_{2}$ and $h_{2,\varepsilon}$, the implied-barrier inversion, and the ratio statistic $Z_{\varepsilon}$ with $\mathrm{AVAR}_{r}$.

---

## License & contact

License: see `LICENSE` (if present); otherwise, all rights reserved by the authors.
Questions or suggestions: open an issue or contact the authors.
