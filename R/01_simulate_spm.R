# =====================================================================
# 01_simulate_spm.R
# Simulate data from the Yashin/Arbeev Stochastic Process Model (SPM) /
# quadratic-hazard model, matching the structure in Arbeev et al. (2023).
#
# Two coupled processes on a STANDARDIZED biomarker scale (as in the paper):
#   (1) Latent mean-reverting (Ornstein-Uhlenbeck) biomarker SDE:
#         dY(t) = a(t) * (Y(t) - f1(t)) dt + b(t) dW(t),     a(t) < 0
#   (2) Quadratic hazard for the time-to-event outcome:
#         mu(t, Y) = mu0(t) + Q(t) * (Y(t) - f0(t))^2
#
# Interpretable components (paper's terminology):
#   f1(t)  allostatic mean trajectory (where the body is pulled toward)
#   a(t)   adaptive capacity / resilience  (mean-reversion strength; |a| down with age)
#   b(t)   diffusion (variability), here constant = sigma1
#   f0(t)  OPTIMAL / risk-minimizing biomarker value  ("physiological norm")
#   Q(t)   robustness / stress-resistance (U-shape curvature; larger Q = steeper risk)
#   mu0(t) baseline hazard (Gompertz-like)
#   AL(t) = f0(t) - f1(t)   "allostatic load" gap
#
# Simulation strategy: fine-grid Euler-Maruyama for Y with a competing-risk
# event draw at each step (hazard -> per-step event prob). Longitudinal Y is
# then "observed" on an irregular schedule with measurement error, so missing /
# irregular observations are built in (mirrors CHS-style data).
# =====================================================================

set.seed(20260617)

# Generative simulator + component functions live in R/spm_simulate.R (shared
# with the calibration study so the data-generating process is defined once).
source("R/spm_simulate.R")

# ---- Run & summarize --------------------------------------------------------
p   <- make_spm_params()
sim <- simulate_spm(N = 2000, p = p)

cat("=== Simulated SPM cohort ===\n")
cat(sprintf("Individuals:           %d\n", nrow(sim$surv)))
cat(sprintf("Events:                %d (%.1f%%)\n",
            sum(sim$surv$status), 100 * mean(sim$surv$status)))
cat(sprintf("Longitudinal obs:      %d (median %.0f / person)\n",
            nrow(sim$long), median(table(sim$long$id))))
cat(sprintf("Entry age range:       %.1f - %.1f\n",
            min(sim$surv$t0), max(sim$surv$t0)))
cat(sprintf("Event/censor age range:%.1f - %.1f\n",
            min(sim$surv$time), max(sim$surv$time)))
cat("\nLongitudinal head:\n"); print(head(sim$long))
cat("\nSurvival head:\n");     print(head(sim$surv))

# Save for the estimation step
saveRDS(sim, "output/sim_spm.rds")
cat("\nSaved -> output/sim_spm.rds\n")

# Quick descriptive: empirical U-shape check (event rate by biomarker deviation)
last_obs <- do.call(rbind, lapply(split(sim$long, sim$long$id), function(d) d[nrow(d), ]))
merged   <- merge(last_obs, sim$surv, by = "id")
merged$dev <- merged$y - f0_fun(merged$age, p)
brks <- quantile(merged$dev, probs = seq(0, 1, 0.1), na.rm = TRUE)
merged$bin <- cut(merged$dev, breaks = brks, include.lowest = TRUE)
ushape <- aggregate(status ~ bin, data = merged, FUN = mean)
cat("\nEmpirical event rate by last-obs deviation-from-optimal decile (expect U/J):\n")
print(ushape)
