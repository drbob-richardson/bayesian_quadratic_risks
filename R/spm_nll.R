# =====================================================================
# spm_nll.R  -- shared SPM machinery used by both the MLE (RTMB) and the
# Bayesian (tmbstan) fits, so the likelihood is defined in ONE place.
#
#   pack_data(sim, nsub, use_priors)  -> data list for MakeADFun
#   make_nll(dat)                     -> nll/negative-log-posterior closure
#
# The nll is the analytic marginalized SPM likelihood (Yashin 2007 moment
# recursion). When dat$use_priors == 1 it ADDS negative-log-prior terms, so
# the objective becomes the negative log POSTERIOR -- exactly what tmbstan
# should sample. With use_priors == 0 it is the pure negative log-likelihood
# for maximum-likelihood estimation.
# =====================================================================

suppressPackageStartupMessages(library(RTMB))

# Pack a simulated (or real) dataset into the per-individual list MakeADFun needs.
# Expects sim$long (id, age, y) and sim$surv (id, time, status); tmin scalar.
pack_data <- function(sim, tmin = 65, nsub = 4L, use_priors = 0L,
                       ids = NULL, agrid = seq(65, 95, by = 2.5)) {
  if (is.null(ids)) ids <- sort(unique(sim$surv$id))
  nid <- length(ids)
  age_l <- vector("list", nid); y_l <- vector("list", nid)
  for (k in seq_len(nid)) {
    d <- sim$long[sim$long$id == ids[k], ]
    age_l[[k]] <- d$age; y_l[[k]] <- d$y
  }
  list(nid = nid, age = age_l, y = y_l,
       tend   = sim$surv$time[match(ids, sim$surv$id)],
       status = sim$surv$status[match(ids, sim$surv$id)],
       tmin   = tmin, nsub = as.integer(nsub),
       use_priors = as.integer(use_priors),
       agrid = agrid)
}

# Default starting / reference parameter list (transformed scale).
spm_init <- function() list(
  laY   = log(0.12),  bY = 0.0,
  af1   = 0.1,        bf1 = 0.01,
  af0   = 0.1,        bf0 = 0.0,
  laQ   = log(0.03),  bQ  = 0.0,
  lamu0 = log(0.002), bmu0 = 0.05,
  lsig1 = log(0.3),
  lsig0 = log(0.8),
  ltau  = log(0.15)
)

make_nll <- function(dat) {
  function(parms) {
    getAll(parms, dat)
    sig1 <- exp(lsig1); sig0 <- exp(lsig0); tau <- exp(ltau)
    sig1sq <- sig1^2; tausq <- tau^2

    a_f   <- function(t) -exp(laY) + bY * (t - tmin)
    f1_f  <- function(t) af1 + bf1 * (t - tmin)
    f0_f  <- function(t) af0 + bf0 * (t - tmin)
    Q_f   <- function(t) exp(laQ) + bQ * (t - tmin)   # intercept > 0
    mu0_f <- function(t) exp(lamu0 + bmu0 * (t - tmin))

    # one RK4 step of (m,g) + trapezoidal cum-hazard increment; returns list
    step <- function(m, g, t, dtv) {
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
      m_n <- m + (dm1+2*dm2+2*dm3+dm4)*dtv/6
      g_n <- g + (dg1+2*dg2+2*dg3+dg4)*dtv/6
      h1 <- mu0_f(te) + Q_f(te) * ((m_n - f0_f(te))^2 + g_n)
      list(m = m_n, g = g_n, dLam = (h0 + h1)/2 * dtv)
    }

    nll <- 0
    for (i in 1:nid) {
      ag <- age[[i]]; yy <- y[[i]]; nob <- length(ag)
      m <- f1_f(ag[1]); g <- sig0^2; Lambda <- 0
      for (j in 1:nob) {
        if (j > 1) {
          t <- ag[j-1]; dtv <- (ag[j]-ag[j-1])/nsub
          for (s in 1:nsub) { st <- step(m,g,t,dtv); m <- st$m; g <- st$g
            Lambda <- Lambda + st$dLam; t <- t + dtv }
        }
        nll <- nll - dnorm(yy[j], m, sqrt(g + tausq), log = TRUE)
        K <- g/(g + tausq); m <- m + K*(yy[j]-m); g <- (1-K)*g
      }
      t <- ag[nob]; te_end <- tend[i]
      if (te_end > t) {
        dtv <- (te_end - t)/nsub
        for (s in 1:nsub) { st <- step(m,g,t,dtv); m <- st$m; g <- st$g
          Lambda <- Lambda + st$dLam; t <- t + dtv }
      }
      nll <- nll + Lambda
      if (status[i] == 1) {
        mubar <- mu0_f(te_end) + Q_f(te_end) * ((m - f0_f(te_end))^2 + g)
        nll <- nll - log(mubar)
      }
    }

    # ---- weakly-informative priors (added only for Bayesian sampling) -------
    if (use_priors == 1) {
      lp <- 0
      lp <- lp + dnorm(laY,   log(0.10), 0.70, log=TRUE)   # mean-reversion scale
      lp <- lp + dnorm(bY,    0.0,       0.010, log=TRUE)
      lp <- lp + dnorm(af1,   0.0,       1.00, log=TRUE)
      lp <- lp + dnorm(bf1,   0.0,       0.050, log=TRUE)
      lp <- lp + dnorm(af0,   0.0,       1.00, log=TRUE)
      lp <- lp + dnorm(bf0,   0.0,       0.050, log=TRUE)
      lp <- lp + dnorm(laQ,   log(0.03), 1.00, log=TRUE)
      lp <- lp + dnorm(bQ,    0.0,       0.010, log=TRUE)
      lp <- lp + dnorm(lamu0, log(0.002),1.50, log=TRUE)   # baseline (weakly id.)
      lp <- lp + dnorm(bmu0,  0.06,      0.050, log=TRUE)   # Gompertz slope prior
      lp <- lp + dnorm(lsig1, log(0.30), 0.50, log=TRUE)
      lp <- lp + dnorm(lsig0, log(0.80), 0.50, log=TRUE)
      lp <- lp + dnorm(ltau,  log(0.15), 0.50, log=TRUE)
      nll <- nll - lp
    }

    # f0(t) on the age grid -> ADREPORT (SE bands for MLE; tracked for both)
    f0_grid <- af0 + bf0 * (agrid - tmin)
    ADREPORT(f0_grid)
    nll
  }
}
