# Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times — MATLAB Code

This repository contains the MATLAB code and simulation examples of **“Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.”** You can find the latest draft and supplemental materials on my website.

---

## What’s inside

* **Test implementation (`testLLNNY.m`)**
  Builds the test statistic from **implied barriers** $\widehat m$ and $\widehat m_{\varepsilon}$ by inverting the functions $h_{2}(m)$ and $h_{2,\varepsilon}(m)$.

* **Precomputed tables**:
  `h_vec.mat` $[\;m,\, h_{2}(m)\;]$, `h_eps_vec.mat` $[\;m,\, h_{2,\epsilon}(m)\;]$, `avar_r.mat` $[\;m,\, \mathrm{AVAR}_{r}(m;\varepsilon)\;]$.

* **Simulators**
  Heston model with jumps, with/without market microstructure noise.

* **Monte Carlo examples**
  Size and (size-adjusted) power across various barrier width $c$ and censoring parameter $\epsilon$.

---

## Folder structure

```
numerical_h_hbar/
  h_simulate.m                % build h2, h2,ε, derivatives, AVAR_r
  h_first_derivative.m
  ret_delta.m
  h_vec.mat
  h_eps_vec.mat
  avar_r.mat

testToolbox/
  testLLNNY.m                 % ratio test Z_ε
  invFunc.m                   % monotone inverse on [m, Y(m)]
  linearInt.m                 % safe interp1 wrapper
  ret_delta.m                 % barrier-hitting returns r^{(c)}
  finddata.m                  % pick ε column → two-col tables
  wb_preaveraging.m           % pre-averaging (+ optional wild bootstrap)
  pseudo.m                    % thin wrapper to wb_preaveraging (compat)

MCExamples/
  simulatePrices/
    simPriceEfficient.m
    simPriceNoise_autoGau_t.m
    MA_noise.m, rounding.m, randLaplace.m (if used)
  mainNoNoise.m
  mainWithNoise.m

LLNNY_Jumps.pdf               % paper (for notation and theory)
```

---

## Requirements

* MATLAB R2018b+ (tested on R2021a+).
* Optional: Parallel Computing Toolbox for `parfor`.

---

## Install

```matlab
cd <path-to-repo>
addpath(genpath(pwd));
savepath;    % optional
```

---

## Quick start

```matlab
% Load tables
load numerical_h_hbar/h_vec.mat
load numerical_h_hbar/h_eps_vec.mat
load numerical_h_hbar/avar_r.mat

% Pick ε and slice two-column tables
eps = 0.05;
[H2eps_tab, AVAR_tab] = finddata(eps, h_eps_vec, avar_r);
H2_tab = h_vec;

% A path X (column vector)
X = cumsum(randn(50_000,1));

% Working barrier: c = K * sqrt(Var(ΔX))  (finite-sample rule)
K = 6;
c = K * sqrt(var(diff(X), 1));   % population variance

% Ratio test
[Z, mhat_eps, mhat, N_c] = testLLNNY(X, c, eps, H2_tab, H2eps_tab, AVAR_tab);

% 5% right-tail decision
reject = (Z > 1.64485362695147);
```

**Tips**

* Choose $K$ so the number of hits $N_c$ is comfortably large (hundreds+).
* With microstructure noise, you may **pre-average** only for computing $c$:
  `X_pa = wb_preaveraging(X, round(0.5*sqrt(length(X))));  c = K*sqrt(var(diff(X_pa),1));`

---

## Reproducing the tables

Run `numerical_h_hbar/h_simulate.m` to compute $h_{2}$, $h_{2,\varepsilon}$, their local-linear derivatives, and $\mathrm{AVAR}_{r}$ on the $m$-grid, then save to the `.mat` files used by the test.

---

## Monte Carlo

See `MCExamples/mainNoNoise.m` and `mainWithNoise.m` for size and **size-adjusted** power. For jump alternatives, it’s standard to **condition on paths with at least one jump**.

---

## Citation

If you use this code, please cite the paper “Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.” (See `LLNNY_Jumps.pdf` and the latest draft on my website.)
