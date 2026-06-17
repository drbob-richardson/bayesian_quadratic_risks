# =====================================================================
# 03_fit_tmbstan.R
# Bayesian fit of the SPM / quadratic-hazard model: NUTS (via tmbstan) on
# the SAME RTMB objective used for MLE, now with weakly-informative priors.
# Deliverable 1: POSTERIOR CREDIBLE BANDS on the optimal trajectory f0(t).
#
# Config via env vars (with sane defaults for a proof-of-concept run):
#   SPM_N      number of individuals to use (subsample)   [500]
#   SPM_NSUB   ODE substeps per interval                   [3]
#   SPM_ITER   total iterations per chain                  [600]
#   SPM_WARM   warmup iterations                           [300]
#   SPM_CHAINS chains                                      [2]
# =====================================================================

suppressPackageStartupMessages({
  library(RTMB); library(tmbstan); library(rstan)
})
source("R/spm_nll.R")

N      <- as.integer(Sys.getenv("SPM_N",      "500"))
NSUB   <- as.integer(Sys.getenv("SPM_NSUB",   "3"))
ITER   <- as.integer(Sys.getenv("SPM_ITER",   "600"))
WARM   <- as.integer(Sys.getenv("SPM_WARM",   "300"))
CHAINS <- as.integer(Sys.getenv("SPM_CHAINS", "2"))
options(mc.cores = min(CHAINS, parallel::detectCores()))

sim  <- readRDS("output/sim_spm.rds")
ptru <- sim$params; tmin <- ptru$tmin

set.seed(1)
all_ids <- sort(unique(sim$surv$id))
use_ids <- sort(sample(all_ids, min(N, length(all_ids))))
dat <- pack_data(sim, tmin = tmin, nsub = NSUB, use_priors = 1L, ids = use_ids)
cat(sprintf("Bayesian SPM: N=%d  events=%d  nsub=%d  %d chains x %d iter (warmup %d)\n",
            dat$nid, sum(dat$status), NSUB, CHAINS, ITER, WARM))

obj <- MakeADFun(make_nll(dat), spm_init(), silent = TRUE)

# warm start at the posterior mode (fast) so all chains begin in a good region
cat("Finding posterior mode for init...\n")
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(iter.max = 500, eval.max = 500))
cat(sprintf("  mode nlp=%.2f (code %d)\n", opt$objective, opt$convergence))

cat("Sampling (NUTS)...\n")
t0 <- Sys.time()
fit <- tmbstan(obj, chains = CHAINS, iter = ITER, warmup = WARM,
               init = "last.par.best", seed = 123,
               control = list(adapt_delta = 0.9, max_treedepth = 12))
cat(sprintf("Sampling done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Diagnostics ------------------------------------------------------------
post <- as.matrix(fit)
ndiv <- sum(sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                   function(x) sum(x[, "divergent__"])))
mon  <- rstan::monitor(fit, print = FALSE)
cat(sprintf("\nDivergences: %d | min Bulk-ESS: %.0f | max Rhat: %.3f\n",
            ndiv, min(mon$Bulk_ESS, na.rm = TRUE), max(mon$Rhat, na.rm = TRUE)))

# ---- Posterior on parameters (natural scale) vs truth -----------------------
draws <- list(
  aY     = -exp(post[,"laY"]),  bY = post[,"bY"],
  af1    = post[,"af1"],        bf1 = post[,"bf1"],
  af0    = post[,"af0"],        bf0 = post[,"bf0"],
  aQ     = exp(post[,"laQ"]),   bQ  = post[,"bQ"],
  amu0   = exp(post[,"lamu0"]), bmu0 = post[,"bmu0"],
  sigma1 = exp(post[,"lsig1"]), sigma0 = exp(post[,"lsig0"]),
  tau    = exp(post[,"ltau"]))
truth <- list(aY=ptru$aY,bY=ptru$bY,af1=ptru$af1,bf1=ptru$bf1,af0=ptru$af0,
              bf0=ptru$bf0,aQ=ptru$aQ,bQ=ptru$bQ,amu0=ptru$amu0,bmu0=ptru$bmu0,
              sigma1=ptru$sigma1,sigma0=ptru$sigma0,tau=ptru$tau)
cat("\n=== Posterior (median [95% CrI]) vs truth ===\n")
for (nm in names(draws)) {
  q <- quantile(draws[[nm]], c(.5,.025,.975))
  hit <- ifelse(truth[[nm]] >= q[2] & truth[[nm]] <= q[3], "", "  <-- MISS")
  cat(sprintf("  %-8s %9.4f [%9.4f, %9.4f]  truth=%9.4f%s\n",
              nm, q[1], q[2], q[3], truth[[nm]], hit))
}

# ---- CREDIBLE BANDS on f0(t) ------------------------------------------------
agrid <- dat$agrid
f0_draws <- sapply(agrid, function(t) post[,"af0"] + post[,"bf0"]*(t - tmin))
f0_q <- apply(f0_draws, 2, quantile, c(.5,.025,.975))
f0_true <- ptru$af0 + ptru$bf0*(agrid - tmin)
cat("\n=== f0(t) POSTERIOR median [95% credible band] vs truth ===\n")
for (k in seq_along(agrid)) {
  hit <- ifelse(f0_true[k] >= f0_q[2,k] & f0_true[k] <= f0_q[3,k], "", "  <-- MISS")
  cat(sprintf("  age %4.1f : %7.3f [%7.3f, %7.3f]  truth=%7.3f%s\n",
              agrid[k], f0_q[1,k], f0_q[2,k], f0_q[3,k], f0_true[k], hit))
}

saveRDS(list(fit = fit, draws = draws, f0_q = f0_q, agrid = agrid,
             use_ids = use_ids, config = list(N=N,NSUB=NSUB,ITER=ITER,
             WARM=WARM,CHAINS=CHAINS)), "output/fit_tmbstan.rds")
cat("\nSaved -> output/fit_tmbstan.rds\n")
