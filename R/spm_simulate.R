# =====================================================================
# spm_simulate.R -- shared generative simulator for the SPM / quadratic-hazard
# model (fine-grid Euler-Maruyama for the latent biomarker + competing-risk event
# draw; irregular noisy longitudinal observations). Sourced by 01 and by the
# calibration study so the data-generating process lives in one place.
# =====================================================================

make_spm_params <- function() {
  list(
    tmin   = 65,
    aY     = -0.15,  bY    =  0.003,   # adaptive capacity a(t)=aY+bY(t-tmin)
    sigma1 =  0.25,  sigma0 = 1.00,    # diffusion; initial-state SD
    af1    =  0.00,  bf1   =  0.020,   # allostatic mean f1(t)
    af0    =  0.00,  bf0   =  0.000,   # optimal f0(t) (flat -> AL grows with age)
    aQ     =  0.020, bQ    =  0.0015,  # curvature Q(t)
    amu0   =  0.001, bmu0  =  0.060,   # Gompertz baseline mu0(t)
    tau    =  0.10                     # measurement-error SD
  )
}

a_fun   <- function(t, p) p$aY + p$bY * (t - p$tmin)
f1_fun  <- function(t, p) p$af1 + p$bf1 * (t - p$tmin)
f0_fun  <- function(t, p) p$af0 + p$bf0 * (t - p$tmin)
Q_fun   <- function(t, p) pmax(p$aQ + p$bQ * (t - p$tmin), 1e-6)
mu0_fun <- function(t, p) p$amu0 * exp(p$bmu0 * (t - p$tmin))
haz_fun <- function(t, y, p) mu0_fun(t, p) + Q_fun(t, p) * (y - f0_fun(t, p))^2

simulate_spm <- function(N = 2000, p = make_spm_params(),
                         age_entry = c(65, 75), age_max = 95,
                         dt = 1/12, obs_every = 2) {
  long_list <- vector("list", N); surv_list <- vector("list", N)
  for (i in seq_len(N)) {
    t0 <- runif(1, age_entry[1], age_entry[2])
    y  <- rnorm(1, mean = f1_fun(t0, p), sd = p$sigma0)
    obs_ages <- seq(t0, age_max, by = obs_every)
    t <- t0; event_time <- NA_real_; status <- 0L
    obs_age_rec <- numeric(0); obs_y_rec <- numeric(0); next_obs_idx <- 1L
    while (t < age_max) {
      while (next_obs_idx <= length(obs_ages) && t >= obs_ages[next_obs_idx]) {
        obs_age_rec <- c(obs_age_rec, obs_ages[next_obs_idx])
        obs_y_rec   <- c(obs_y_rec,   y + rnorm(1, 0, p$tau))
        next_obs_idx <- next_obs_idx + 1L
      }
      h <- haz_fun(t, y, p)
      if (runif(1) < 1 - exp(-h * dt)) { event_time <- t + runif(1, 0, dt); status <- 1L; break }
      y <- y + a_fun(t, p) * (y - f1_fun(t, p)) * dt + p$sigma1 * sqrt(dt) * rnorm(1)
      t <- t + dt
    }
    if (status == 0L) event_time <- age_max
    keep <- obs_age_rec <= event_time
    long_list[[i]] <- data.frame(id = i, age = obs_age_rec[keep], y = obs_y_rec[keep])
    surv_list[[i]] <- data.frame(id = i, t0 = t0, time = event_time, status = status)
  }
  list(long = do.call(rbind, long_list), surv = do.call(rbind, surv_list), params = p)
}
