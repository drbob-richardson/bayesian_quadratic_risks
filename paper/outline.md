# Paper 1 — reframing for *Insurance: Mathematics and Economics* (IME)

**Working title:**
*A Bayesian quadratic-hazard model for biomarker-driven mortality and morbidity
risk: optimal health levels with uncertainty.*

Alternatives:
- *Risk-minimizing biological levels with credible bands: a Bayesian stochastic
  process model for longitudinal biomarkers and insured lives.*
- *Dynamic biomarker risk for life and health insurance: a Bayesian
  quadratic-hazard approach.*

**Keywords:** Longevity and mortality risk; Biomarkers; Joint models;
Bayesian inference; Uncertainty quantification; Underwriting.
**MSC:** 91G05 (actuarial mathematics); 62F15 (Bayesian inference); 62N02
(survival/event-history). **JEL:** G22; C11.

---

## The reframe in one paragraph
Same model and machinery we have built; the *story* moves from "aging biology"
to "insurance risk." The pitch to IME: insurers increasingly price and underwrite
on repeatedly measured biological markers, but standard actuarial hazard models
(proportional-hazards GLMs, aggregate mortality models) impose monotone or
multiplicative marker effects and cannot represent the **U-/J-shaped** risk that
is pervasive for BMI, blood pressure, cholesterol and glucose. We give a joint
model of the biomarker's irregular longitudinal dynamics and a **quadratic
hazard** whose minimum is an **age-varying risk-minimizing level `f0(t)`**, and a
**Bayesian** treatment that turns the model's latent risk quantities into
objects with **full posterior uncertainty** — directly usable where uncertainty
is the point (risk margins, capital, dynamic underwriting).

## Why IME, and what its readers will want
- Frame around **actuarial decisions**, not biology: underwriting, health/LTC/CI
  pricing, wellness-program targeting, and **uncertainty margins** (Solvency II
  risk margin, IFRS 17 risk adjustment) — credible bands map onto these.
- Mathematical rigor: clean statement of the model, a proposition for the
  marginalized likelihood (the Kalman-type moment recursion), explicit hazard.
- A reproducible numerical section (simulation + public data); the restricted
  CHS analysis is the companion applied paper (Paper 2), referenced not included.
- Position against actuarial precedent: the SPM/quadratic-hazard family has
  appeared in *North American Actuarial Journal* (Yashin et al. 2012, 2016) — all
  **frequentist**; the Bayesian contribution is new and aimed squarely at the
  uncertainty needs of modern risk management.

## Proposed abstract (draft)
> Insurers increasingly underwrite and price on repeatedly measured biological
> markers — body mass index, blood pressure, glucose, lipids — yet standard
> actuarial hazard models impose monotone or multiplicative marker effects that
> cannot capture the U- and J-shaped mortality and morbidity risks pervasive in
> epidemiology. We develop a Bayesian quadratic-hazard stochastic process model
> that jointly describes (i) the irregular longitudinal dynamics of a biomarker
> through a mean-reverting (Ornstein–Uhlenbeck) diffusion and (ii) a quadratic
> hazard whose minimum identifies an age-varying risk-minimizing level f0(t) — a
> "physiological optimum" — with curvature Q(t) and mean-reversion a(t)
> interpretable as stress-resistance and adaptive capacity. Estimation exploits
> the analytically marginalized likelihood, a Kalman-type moment recursion that
> integrates the quadratic hazard against the Gaussian latent state in closed
> form; we embed this in a Bayesian framework to obtain full posterior
> uncertainty. The treatment delivers (a) credible bands on the optimal
> trajectory f0(t) and the latent risk components, and (b) posterior-predictive,
> dynamically updated survival and hazard forecasts with coherent uncertainty —
> quantities of direct relevance to biomarker-based underwriting, wellness-program
> design, and the uncertainty margins required under modern solvency and reporting
> regimes. We show that priors are essential for identifiability rather than
> cosmetic: on monotone-risk markers the maximum-likelihood estimator collapses
> the quadratic term, while weakly-informative priors recover a stable,
> interpretable fit. The methods are implemented in R, validated by simulation,
> and illustrated on publicly available longitudinal survival data.

## Section plan (and what we already have)
1. **Introduction.** Biomarker underwriting & health insurance trends; failure of
   monotone/PH marker effects for U-shaped risk; longitudinal + irregular
   measurement; the case for joint modeling **with uncertainty**; contributions.
2. **The model.** SDE for the latent biomarker (OU, mean-reverting) + quadratic
   hazard `mu(t,Y)=mu0(t)+Q(t)(Y-f0(t))^2`; interpretation of f0/Q/a/f1 in
   *risk* terms (optimal level, sensitivity/robustness, recovery rate, drift).
   Component parameterizations.
3. **Likelihood (marginalized).** Proposition: conditional Gaussianity of the
   latent state under survival, the moment ODEs (m,γ) absorbing selection, and
   the effective hazard `mu0+Q[(m-f0)^2+γ]`. *(We have this implemented in
   `R/spm_nll.R`.)* Measurement error τ as the general case.
4. **Bayesian inference.** Priors and why they matter (identifiability;
   the MLE quadratic-collapse result); computation via the marginal likelihood
   in HMC/NUTS (RTMB→tmbstan); model comparison (WAIC/LOO) replacing AIC.
5. **Actuarial quantities.** (a) Credible bands on f0(t) — the "optimal metric"
   and its uncertainty; (b) posterior-predictive dynamic survival/hazard via
   landmarking — underwriting/reserving use; (c) tabulated survival/risk by
   deviation-from-optimal with credible intervals.
6. **Simulation study.** Parameter recovery (have it), calibration/coverage of
   credible bands, and the role of priors for identifiability.
7. **Illustration on public data.** PBC (`survival::pbcseq`): worked example on
   real, irregular, partially-missing data; MAP vs collapse; f0(t) on the
   marker scale. *(Have it.)*
8. **Discussion.** Risk margins/Solvency II/IFRS 17; wellness incentives;
   limitations (single biomarker, parametric components); extensions
   (multivariate biomarkers, covariates/underwriting factors, forecasting);
   pointer to the CHS case study (Paper 2).

## Existing assets → manuscript mapping
| Asset (this repo) | Goes into |
|---|---|
| `01_simulate_spm.R` + recovery (`02`) | §6 simulation |
| `spm_nll.R` marginalized likelihood | §3 proposition + §4 computation |
| `03_fit_tmbstan.R` credible bands on f0(t) | §5(a), Fig: `f0_credible_bands.png` |
| `05_risk_prediction.R` survival-by-deviation + dynamic | §5(b,c), Figs: `survival_by_deviation.png`, `dynamic_prediction.png` |
| `06_realdata_pbc.R` | §7, Fig: `pbc_f0_bilirubin.png` |

## To-write / to-check before submission
- Proposition + short proof/derivation of the moment recursion (cite
  Woodbury–Manton 1977; Yashin et al. 2007, 2012 for the ML version).
- Simulation-based **coverage** of the f0(t) credible bands (calibration table).
- Tighten differentiation from spline-based Bayesian joint models (Köhler et al.;
  `bamlss`) — emphasize the interpretable, structured f0/Q/a parametrization.
- Decide secondary illustration (PBC alone, or add a second public marker).
- Confirm IME formatting (elsarticle, structured abstract not required).
