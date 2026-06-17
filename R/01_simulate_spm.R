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

# ---- Component functions (paper structure; tmin = 65 = Medicare eligibility) ----
make_spm_params <- function() {
  list(
    tmin   = 65,
    # adaptive capacity a(t) = aY + bY*(t-tmin), aY<0, bY>=0 -> |a| shrinks with age
    aY     = -0.15,
    bY     =  0.003,
    # diffusion (constant)
    sigma1 =  0.25,
    # initial condition Y(t0) ~ N(f1(t0), sigma0^2)
    sigma0 =  1.00,
    # allostatic mean f1(t) = af1 + bf1*(t-tmin)
    af1    =  0.00,
    bf1    =  0.020,
    # optimal/risk-minimizing f0(t) = af0 + bf0*(t-tmin)
    af0    =  0.00,
    bf0    =  0.000,   # optimal stays flat -> allostatic load AL = f1 - f0 grows with age
    # quadratic-hazard curvature Q(t) = aQ + bQ*(t-tmin) (kept > 0)
    aQ     =  0.020,
    bQ     =  0.0015,
    # baseline hazard mu0(t) = amu0 * exp(bmu0*(t-tmin)) (Gompertz-like)
    amu0   =  0.001,
    bmu0   =  0.060,
    # measurement error sd on observed biomarker
    tau    =  0.10
  )
}

a_fun   <- function(t, p) p$aY + p$bY * (t - p$tmin)
f1_fun  <- function(t, p) p$af1 + p$bf1 * (t - p$tmin)
f0_fun  <- function(t, p) p$af0 + p$bf0 * (t - p$tmin)
Q_fun   <- function(t, p) pmax(p$aQ + p$bQ * (t - p$tmin), 1e-6)   # robustness, >0
mu0_fun <- function(t, p) p$amu0 * exp(p$bmu0 * (t - p$tmin))
haz_fun <- function(t, y, p) mu0_fun(t, p) + Q_fun(t, p) * (y - f0_fun(t, p))^2

# ---- Simulate one cohort ----------------------------------------------------
simulate_spm <- function(N = 2000, p = make_spm_params(),
                         age_entry = c(65, 75),   # uniform entry age range
                         age_max   = 95,           # administrative censoring age
                         dt        = 1/12,         # Euler step (years)
                         obs_every = 2) {          # longitudinal obs spacing (years)
  long_list <- vector("list", N)
  surv_list <- vector("list", N)

  for (i in seq_len(N)) {
    t0 <- runif(1, age_entry[1], age_entry[2])
    y  <- rnorm(1, mean = f1_fun(t0, p), sd = p$sigma0)

    # scheduled observation ages for this individual
    obs_ages <- seq(t0, age_max, by = obs_every)

    t <- t0
    event_time <- NA_real_
    status <- 0L
    obs_age_rec <- numeric(0); obs_y_rec <- numeric(0)
    next_obs_idx <- 1L

    while (t < age_max) {
      # record a (noisy) longitudinal observation if we passed a scheduled age
      while (next_obs_idx <= length(obs_ages) && t >= obs_ages[next_obs_idx]) {
        obs_age_rec <- c(obs_age_rec, obs_ages[next_obs_idx])
        obs_y_rec   <- c(obs_y_rec,   y + rnorm(1, 0, p$tau))
        next_obs_idx <- next_obs_idx + 1L
      }
      # per-step event probability from current hazard
      h <- haz_fun(t, y, p)
      if (runif(1) < 1 - exp(-h * dt)) {
        event_time <- t + runif(1, 0, dt)   # event in this interval
        status <- 1L
        break
      }
      # Euler-Maruyama update of latent biomarker
      y <- y + a_fun(t, p) * (y - f1_fun(t, p)) * dt + p$sigma1 * sqrt(dt) * rnorm(1)
      t <- t + dt
    }
    if (status == 0L) event_time <- age_max   # administrative censoring

    # keep only longitudinal obs strictly before the event/censor time
    keep <- obs_age_rec <= event_time
    long_list[[i]] <- data.frame(id = i,
                                 age = obs_age_rec[keep],
                                 y   = obs_y_rec[keep])
    surv_list[[i]] <- data.frame(id = i, t0 = t0,
                                 time = event_time, status = status)
  }

  long <- do.call(rbind, long_list)
  surv <- do.call(rbind, surv_list)
  list(long = long, surv = surv, params = p)
}

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
