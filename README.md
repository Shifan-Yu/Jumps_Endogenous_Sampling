# Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times — MATLAB Code

This repository contains the MATLAB code and simulation examples of **“Testing for Jumps in a Discretely Observed Price Process with Endogenous Sampling Times.”** You can find the latest draft and supplemental materials on my website.

---

## What’s inside

* **Test implementation (`testLLNNY.m`)**:
  Construct the test statistic from $M_{c}$ and $M_{c,\epsilon}$ by inverting the functions $h_{2}(m)$ and $h_{2,\epsilon}(m)$.

* **Precomputed tables**:
  `h_vec.mat`: $h_{2}(m)$, `h_eps_vec.mat`: $h_{2,\epsilon}(m)$, `avar_r.mat`: $V_{\epsilon}(m)$.

* **Simulators**:
  Simulate a Heston model for the effcient price process and obtain its tick-level observations, to which we add jumps with different sizes.

* **Monte Carlo examples**:
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
  testLLNNY.m                 % test statistic
  invFunc.m                   % monotone inverse on [m, Y(m)]
  linearInt.m                 % safe interp1 wrapper
  ret_delta.m                 % price duration sampling
  finddata.m                  % select the columns (corresponding to select epsilon) in h_vec.mat, h_eps_vec.mat and avar_r.mat
  wb_preaveraging.m           % wild-bootstrapped preaveraging

MCExamples/
  simulatePrices/
    simPriceEfficient.m
    simPriceNoise_autoGau_t.m
    MA_noise.m, rounding.m, randLaplace.m
  mainNoNoise.m
  mainWithNoise.m
