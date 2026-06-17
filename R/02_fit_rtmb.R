# =====================================================================
# 02_fit_rtmb.R
# Maximum-likelihood fit of the SPM / quadratic-hazard model via RTMB,
# using the analytic (marginalized) likelihood — the Yashin (2007)
# Kalman-style moment recursion. NO latent path sampling.
#
# For each individual the conditional moments of the latent biomarker
# given survival + observation history evolve (between observations) as:
#   dm/dt    = a(t)(m - f1(t)) - 2 Q(t) gamma (m - f0(t))
#   dgamma/dt= 2 a(t) gamma + sigma1^2 - 2 Q(t) gamma^2
# and the "effective" hazard given survival + obs is
#   mubar(t) = mu0(t) + Q(t)[ (m - f0(t))^2 + gamma ].
#
# At each observation y_j (with measurement error sd tau):
#   density contribution  N(y_j ; m^-, gamma^- + tau^2)
#   Kalman update         K = gamma^-/(gamma^- + tau^2),
#                         m^+ = m^- + K(y_j - m^-), gamma^+ = (1-K) gamma^-.
# Log-lik = sum log obs-density  - Lambda(t_end)  + status*log mubar(t_end).
#
# Targets: recover the TRUE simulation parameters, especially the optimal
# trajectory f0(t), and produce SE-based bands on f0(t) (the ML analog of
# the Bayesian credible bands we'll build next).
# =====================================================================

suppressPackageStartupMessages(library(RTMB))

sim  <- readRDS("output/sim_spm.rds")
ptru <- sim$params
tmin <- ptru$tmin

# ---- Pack data into per-individual lists -----------------------------------
ids   <- sort(unique(sim$surv$id))
nid   <- length(ids)
age_l <- vector("list", nid)
y_l   <- vector("list", nid)
for (k in seq_len(nid)) {
  d <- sim$long[sim$long$id == ids[k], ]
  age_l[[k]] <- d$age
  y_l[[k]]   <- d$y
}
tend   <- sim$surv$time[match(ids, sim$surv$id)]
status <- sim$surv$status[match(ids, sim$surv$id)]

dat <- list(nid = nid, age = age_l, y = y_l, tend = tend,
            status = status, tmin = tmin, nsub = 4L)  # nsub = ODE substeps/interval

# ---- Parameters (transformed for positivity where needed) -------------------
# a(t)   = -exp(laY) + bY*(t-tmin)     (intercept forced negative)
# f1(t)  = af1 + bf1*(t-tmin)
# f0(t)  = af0 + bf0*(t-tmin)          <-- the optimal trajectory (target)
# Q(t)   = exp(laQ) + bQ*(t-tmin)      (clamped > 0 inside)
# mu0(t) = exp(lamu0 + bmu0*(t-tmin))
# sigma1 = exp(lsig1), sigma0 = exp(lsig0), tau = exp(ltau)
parameters <- list(
  laY   = log(0.12),  bY = 0.0,
  af1   = 0.1,        bf1 = 0.01,
  af0   = 0.1,        bf0 = 0.0,
  laQ   = log(0.03),  bQ  = 0.0,
  lamu0 = log(0.002), bmu0 = 0.05,
  lsig1 = log(0.3),
  lsig0 = log(0.8),
  ltau  = log(0.15)
)

# ---- Negative log-likelihood ------------------------------------------------
nll_fun <- function(parms) {
  getAll(parms, dat)
  sig1 <- exp(lsig1); sig0 <- exp(lsig0); tau <- exp(ltau)
  sig1sq <- sig1^2; tausq <- tau^2

  a_f   <- function(t) -exp(laY) + bY * (t - tmin)
  f1_f  <- function(t) af1 + bf1 * (t - tmin)
  f0_f  <- function(t) af0 + bf0 * (t - tmin)
  Q_f   <- function(t) exp(laQ) + bQ * (t - tmin)   # intercept > 0 by construction
  mu0_f <- function(t) exp(lamu0 + bmu0 * (t - tmin))

  nll <- 0
  for (i in 1:nid) {
    ag <- age[[i]]; yy <- y[[i]]; nob <- length(ag)

    # initial prior at first observation age
    m <- f1_f(ag[1]); g <- sig0^2
    Lambda <- 0

    for (j in 1:nob) {
      if (j > 1) {
        # integrate moments + cumulative hazard from ag[j-1] to ag[j]
        t <- ag[j - 1]; dtv <- (ag[j] - ag[j - 1]) / nsub
        for (s in 1:nsub) {
          h0 <- mu0_f(t) + Q_f(t) * ((m - f0_f(t))^2 + g)
          # RK4 for (m,g)
          dm1 <- a_f(t)*(m-f1_f(t)) - 2*Q_f(t)*g*(m-f0_f(t))
          dg1 <- 2*a_f(t)*g + sig1sq - 2*Q_f(t)*g^2
          tm <- t + dtv/2; m2 <- m+dm1*dtv/2; g2 <- g+dg1*dtv/2
          dm2 <- a_f(tm)*(m2-f1_f(tm)) - 2*Q_f(tm)*g2*(m2-f0_f(tm))
          dg2 <- 2*a_f(tm)*g2 + sig1sq - 2*Q_f(tm)*g2^2
          m3 <- m+dm2*dtv/2; g3 <- g+dg2*dtv/2
          dm3 <- a_f(tm)*(m3-f1_f(tm)) - 2*Q_f(tm)*g3*(m3-f0_f(tm))
          dg3 <- 2*a_f(tm)*g3 + sig1sq - 2*Q_f(tm)*g3^2
          te <- t+dtv; m4 <- m+dm3*dtv; g4 <- g+dg3*dtv
          dm4 <- a_f(te)*(m4-f1_f(te)) - 2*Q_f(te)*g4*(m4-f0_f(te))
          dg4 <- 2*a_f(te)*g4 + sig1sq - 2*Q_f(te)*g4^2
          m <- m + (dm1+2*dm2+2*dm3+dm4)*dtv/6
          g <- g + (dg1+2*dg2+2*dg3+dg4)*dtv/6
          h1 <- mu0_f(te) + Q_f(te) * ((m - f0_f(te))^2 + g)
          Lambda <- Lambda + (h0 + h1)/2 * dtv     # trapezoidal cum-hazard
          t <- te
        }
      }
      # observation density + Kalman update
      nll <- nll - dnorm(yy[j], m, sqrt(g + tausq), log = TRUE)
      K <- g / (g + tausq)
      m <- m + K * (yy[j] - m)
      g <- (1 - K) * g
    }

    # tail: last obs age -> t_end, accumulate hazard, add event term
    t <- ag[nob]; te_end <- tend[i]
    if (te_end > t) {
      dtv <- (te_end - t) / nsub
      for (s in 1:nsub) {
        h0 <- mu0_f(t) + Q_f(t) * ((m - f0_f(t))^2 + g)
        dm1 <- a_f(t)*(m-f1_f(t)) - 2*Q_f(t)*g*(m-f0_f(t))
        dg1 <- 2*a_f(t)*g + sig1sq - 2*Q_f(t)*g^2
        tm <- t+dtv/2; m2 <- m+dm1*dtv/2; g2 <- g+dg1*dtv/2
        dm2 <- a_f(tm)*(m2-f1_f(tm)) - 2*Q_f(tm)*g2*(m2-f0_f(tm))
        dg2 <- 2*a_f(tm)*g2 + sig1sq - 2*Q_f(tm)*g2^2
        m3 <- m+dm2*dtv/2; g3 <- g+dg2*dtv/2
        dm3 <- a_f(tm)*(m3-f1_f(tm)) - 2*Q_f(tm)*g3*(m3-f0_f(tm))
        dg3 <- 2*a_f(tm)*g3 + sig1sq - 2*Q_f(tm)*g3^2
        te <- t+dtv; m4 <- m+dm3*dtv; g4 <- g+dg3*dtv
        dm4 <- a_f(te)*(m4-f1_f(te)) - 2*Q_f(te)*g4*(m4-f0_f(te))
        dg4 <- 2*a_f(te)*g4 + sig1sq - 2*Q_f(te)*g4^2
        m <- m + (dm1+2*dm2+2*dm3+dm4)*dtv/6
        g <- g + (dg1+2*dg2+2*dg3+dg4)*dtv/6
        h1 <- mu0_f(te) + Q_f(te) * ((m - f0_f(te))^2 + g)
        Lambda <- Lambda + (h0 + h1)/2 * dtv
        t <- te
      }
    }
    nll <- nll + Lambda
    if (status[i] == 1) {
      mubar <- mu0_f(te_end) + Q_f(te_end) * ((m - f0_f(te_end))^2 + g)
      nll <- nll - log(mubar)
    }
  }

  # ADREPORT f0(t) on an age grid -> SE-based bands (ML analog of credible bands)
  agrid <- seq(65, 95, by = 2.5)
  f0_grid <- af0 + bf0 * (agrid - tmin)
  ADREPORT(f0_grid)
  nll
}

# ---- Fit --------------------------------------------------------------------
cat("Building AD function (taping)...\n")
obj <- MakeADFun(nll_fun, parameters, silent = TRUE)
cat("Optimizing...\n")
t0 <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(iter.max = 500, eval.max = 500))
cat(sprintf("Done in %.1f s. Convergence code: %d (%s)\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            opt$convergence, opt$message))
cat(sprintf("Final nll: %.3f\n", opt$objective))

sdr <- sdreport(obj)

# ---- Compare estimates vs truth --------------------------------------------
est <- as.list(sdr, "Estimate")
se  <- as.list(sdr, "Std. Error")
trans <- function(nm, val) val
recover <- function(label, est_val, se_val, truth) {
  z <- (est_val - truth) / se_val
  cat(sprintf("  %-8s est=%9.4f  se=%8.4f  truth=%9.4f  z=%6.2f\n",
              label, est_val, se_val, truth, z))
}
cat("\n=== Parameter recovery (natural scale) ===\n")
recover("aY",    -exp(est$laY),  exp(est$laY)*se$laY,          ptru$aY)
recover("bY",     est$bY,        se$bY,                         ptru$bY)
recover("af1",    est$af1,       se$af1,                        ptru$af1)
recover("bf1",    est$bf1,       se$bf1,                        ptru$bf1)
recover("af0",    est$af0,       se$af0,                        ptru$af0)
recover("bf0",    est$bf0,       se$bf0,                        ptru$bf0)
recover("aQ",     exp(est$laQ),  exp(est$laQ)*se$laQ,           ptru$aQ)
recover("bQ",     est$bQ,        se$bQ,                         ptru$bQ)
recover("amu0",   exp(est$lamu0),exp(est$lamu0)*se$lamu0,       ptru$amu0)
recover("bmu0",   est$bmu0,      se$bmu0,                       ptru$bmu0)
recover("sigma1", exp(est$lsig1),exp(est$lsig1)*se$lsig1,       ptru$sigma1)
recover("sigma0", exp(est$lsig0),exp(est$lsig0)*se$lsig0,       ptru$sigma0)
recover("tau",    exp(est$ltau), exp(est$ltau)*se$ltau,         ptru$tau)

# ---- f0(t) band vs truth ----------------------------------------------------
rn <- rownames(summary(sdr, "report"))
f0_summ <- summary(sdr, "report")[grep("f0_grid", rn), ]
agrid <- seq(65, 95, by = 2.5)
f0_true <- ptru$af0 + ptru$bf0 * (agrid - tmin)
cat("\n=== Optimal trajectory f0(t): estimate +/- 1.96 SE vs truth ===\n")
for (k in seq_along(agrid)) {
  e <- f0_summ[k, 1]; s <- f0_summ[k, 2]
  cat(sprintf("  age %4.1f : %7.3f [%7.3f, %7.3f]  truth=%7.3f %s\n",
              agrid[k], e, e-1.96*s, e+1.96*s, f0_true[k],
              ifelse(f0_true[k] >= e-1.96*s & f0_true[k] <= e+1.96*s, "", "  <-- MISS")))
}

saveRDS(list(opt = opt, sdr = sdr, est = est, se = se), "output/fit_rtmb.rds")
cat("\nSaved -> output/fit_rtmb.rds\n")
