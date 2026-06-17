# Bayesian Stochastic Process Model (Quadratic-Hazard) for Health & Actuarial Applications

A Bayesian reformulation of the Yashin/Arbeev **Stochastic Process Model (SPM)** /
**quadratic-hazard model** — a joint model of a longitudinal biomarker and a
time-to-event outcome — implemented in R. The model is fit by maximum likelihood
(RTMB) and, on the *same* objective, by full Bayesian NUTS sampling (tmbstan),
yielding **posterior credible bands on the risk-minimizing "optimal" biomarker
trajectory `f0(t)`** and posterior-predictive risk/survival with uncertainty.

Motivated by Arbeev et al. (2023), *Mech. Ageing Dev.* 211:111791
(`refs/Arbeev_etal_2023.pdf`), which applies the SPM to BMI trajectories and
Alzheimer's onset in the Health and Retirement Study.

## The model

**Latent biomarker (mean-reverting Ornstein–Uhlenbeck SDE):**

```
dY(t) = a(t) [ Y(t) - f1(t) ] dt + b(t) dW(t)
```

**Quadratic hazard for the event:**

```
mu(t, Y) = mu0(t) + Q(t) [ Y(t) - f0(t) ]^2
```

| component | meaning |
|-----------|---------|
| `f0(t)`   | **optimal / risk-minimizing biomarker level** ("physiological norm") — the headline estimand |
| `Q(t)`    | robustness / stress-resistance (curvature of the U-shaped risk) |
| `a(t)`    | adaptive capacity / resilience (mean-reversion strength) |
| `f1(t)`   | allostatic mean trajectory |
| `b(t)`    | diffusion (variability), here constant `sigma1` |
| `mu0(t)`  | baseline hazard (Gompertz-like) |

Estimation uses the **analytic marginalized likelihood** (Yashin 2007 Kalman-style
moment recursion): the conditional mean `m(t)` and variance `gamma(t)` of the latent
biomarker, given survival and observation history, evolve via ODEs that absorb the
survival-selection effect — so the latent path is never sampled.

## Repository layout

```
R/
  spm_nll.R            Shared machinery: data packing + the marginalized
                       likelihood / negative-log-posterior (one definition,
                       used by both the MLE and Bayesian fits).
  01_simulate_spm.R    Simulate an SPM cohort (irregular longitudinal obs + TTE).
  02_fit_rtmb.R        Maximum-likelihood fit (RTMB) + parameter recovery + SE bands.
  03_fit_tmbstan.R     Bayesian fit (NUTS via tmbstan) + credible bands on f0(t).
  04_plot_f0.R         Plot: Bayesian credible band vs ML vs truth for f0(t).
  05_risk_prediction.R Posterior-predictive risk/survival with credible bands.
  06_realdata_pbc.R    Real-data fit: SPM on PBC sequential data
                       (survival::pbcseq), log-bilirubin vs death (monotone).
  07_calibration.R     Frequentist coverage study over replicate cohorts.
  08_realdata_nafld.R  Real-data fit: SPM on NAFLD systolic blood pressure
                       (survival::nafld) vs death (U-shaped; interior optimum).
refs/                  Source paper.
output/                Generated artifacts (*.rds git-ignored); figures/ committed.
```

## Reproducing

Requires R (>= 4.5) with `RTMB`, `tmbstan`, `rstan`, `posterior`, `bayesplot`,
and a C++ toolchain.

```r
install.packages(c("RTMB","tmbstan","rstan","posterior","bayesplot"))
```

Run from the project root, in order:

```sh
Rscript R/01_simulate_spm.R     # -> output/sim_spm.rds
Rscript R/02_fit_rtmb.R         # -> output/fit_rtmb.rds   (MLE recovery)
Rscript R/03_fit_tmbstan.R      # -> output/fit_tmbstan.rds (Bayesian)
Rscript R/04_plot_f0.R          # -> output/figures/f0_credible_bands.png
Rscript R/05_risk_prediction.R  # -> output/figures/*.png
```

`03_fit_tmbstan.R` is configurable via environment variables
(`SPM_N`, `SPM_NSUB`, `SPM_ITER`, `SPM_WARM`, `SPM_CHAINS`) for subsample size,
ODE substeps, and sampler iterations.

## Status

- [x] SPM simulator (clean U-shaped risk, irregular/missing longitudinal obs)
- [x] Marginalized likelihood in RTMB; **all parameters recovered** on simulated data
- [x] Bayesian engine (tmbstan); **credible bands on `f0(t)` cover truth**, priors
      regularize weakly-identified baseline-hazard parameters
- [x] Posterior-predictive risk / survival with uncertainty
- [x] Real-data illustrations on two public cohorts: NAFLD systolic blood
      pressure (`survival::nafld`) shows a genuine U-shape with a firmly
      identified interior optimum (~138 mmHg [137, 139]); PBC bilirubin
      (`pbcseq`) is a monotone stress test where MLE collapses and MAP recovers
      a low-boundary optimum with wide bands
- [ ] Performance pass (vectorize likelihood) + full-scale production run
- [ ] Bayesian (NUTS) fit on real data; case study (CHS)
- [ ] Covariate stratification (e.g. genetic risk factor) of `f0/Q/a`

## Notes

- The classic SPM treats observations as exact; here a measurement-error SD `tau`
  is included as the general case (`tau -> 0` recovers the classic model).
- `stpm` (the original authors' CRAN package) was intended as an MLE cross-check
  but is archived and currently fails to build (links against an absent gfortran
  runtime); our own simulator + likelihood are the primary path.
