# =====================================================================
# 07_calibration.R -- Frequentist coverage study (the validation of record for
# the uncertainty bands). Over R replicate cohorts simulated from a fixed truth,
# fit by maximum likelihood and record whether the 95% Wald intervals contain the
# truth, for every parameter and for the optimal trajectory f0(t) on an age grid.
# Calibrated bands give coverage ~ 0.95.
#
# Config via env vars: SPM_R [200], SPM_N [500], SPM_NSUB [3].
# =====================================================================

suppressPackageStartupMessages(library(RTMB))
source("R/spm_simulate.R")
source("R/spm_nll.R")

R    <- as.integer(Sys.getenv("SPM_R",    "200"))
N    <- as.integer(Sys.getenv("SPM_N",    "500"))
NSUB <- as.integer(Sys.getenv("SPM_NSUB", "3"))
p    <- make_spm_params(); tmin <- p$tmin
agrid <- seq(65, 95, by = 2.5)
cat(sprintf("Coverage study: R=%d replicates, N=%d, nsub=%d\n", R, N, NSUB))

# natural-scale truth for the 13 parameters
truth <- c(aY=p$aY, bY=p$bY, af1=p$af1, bf1=p$bf1, af0=p$af0, bf0=p$bf0,
           aQ=p$aQ, bQ=p$bQ, amu0=p$amu0, bmu0=p$bmu0,
           sigma1=p$sigma1, sigma0=p$sigma0, tau=p$tau)
f0_truth <- p$af0 + p$bf0 * (agrid - tmin)

pnames <- names(truth)
cov_par <- setNames(numeric(length(pnames)), pnames)   # coverage counts
cov_f0  <- numeric(length(agrid))
nbad <- 0L; nok <- 0L
z <- 1.96

for (r in 1:R) {
  set.seed(1000 + r)
  fit_ok <- tryCatch({
    sim <- simulate_spm(N = N, p = p)
    dat <- pack_data(sim, tmin = tmin, nsub = NSUB, use_priors = 0L, agrid = agrid)
    obj <- MakeADFun(make_nll(dat), spm_init(), silent = TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr,
                  control = list(iter.max = 500, eval.max = 500))
    sdr <- sdreport(obj)
    e <- as.list(sdr, "Estimate"); s <- as.list(sdr, "Std. Error")
    # endpoint-transformed 95% Wald intervals on the natural scale
    ci <- list(
      aY     = sort(-exp(c(e$laY - z*s$laY, e$laY + z*s$laY))),
      bY     = e$bY + c(-1,1)*z*s$bY,
      af1    = e$af1 + c(-1,1)*z*s$af1,    bf1 = e$bf1 + c(-1,1)*z*s$bf1,
      af0    = e$af0 + c(-1,1)*z*s$af0,    bf0 = e$bf0 + c(-1,1)*z*s$bf0,
      aQ     = exp(c(e$laQ - z*s$laQ, e$laQ + z*s$laQ)),
      bQ     = e$bQ + c(-1,1)*z*s$bQ,
      amu0   = exp(c(e$lamu0 - z*s$lamu0, e$lamu0 + z*s$lamu0)),
      bmu0   = e$bmu0 + c(-1,1)*z*s$bmu0,
      sigma1 = exp(c(e$lsig1 - z*s$lsig1, e$lsig1 + z*s$lsig1)),
      sigma0 = exp(c(e$lsig0 - z*s$lsig0, e$lsig0 + z*s$lsig0)),
      tau    = exp(c(e$ltau  - z*s$ltau,  e$ltau  + z*s$ltau)))
    # f0(t) bands from ADREPORT
    sm <- summary(sdr, "report"); sm <- sm[grep("f0_grid", rownames(sm)), ]
    f0_lo <- sm[,1] - z*sm[,2]; f0_hi <- sm[,1] + z*sm[,2]
    if (any(!is.finite(unlist(ci))) || any(!is.finite(c(f0_lo,f0_hi))) ||
        opt$convergence != 0) stop("non-finite/non-converged")
    list(ci = ci, f0_lo = f0_lo, f0_hi = f0_hi)
  }, error = function(ee) NULL)

  if (is.null(fit_ok)) { nbad <- nbad + 1L; next }
  nok <- nok + 1L
  for (nm in pnames)
    cov_par[nm] <- cov_par[nm] + (truth[nm] >= fit_ok$ci[[nm]][1] &
                                  truth[nm] <= fit_ok$ci[[nm]][2])
  cov_f0 <- cov_f0 + (f0_truth >= fit_ok$f0_lo & f0_truth <= fit_ok$f0_hi)
  if (r %% 25 == 0) cat(sprintf("  ...%d/%d (ok=%d, failed=%d)\n", r, R, nok, nbad))
}

cat(sprintf("\nCompleted: %d valid fits, %d failed.\n", nok, nbad))
mc_se <- function(k, n) sqrt((k/n)*(1-k/n)/n)
cat("\n=== 95%% Wald interval coverage by parameter ===\n")
for (nm in pnames) {
  cv <- cov_par[nm]/nok
  cat(sprintf("  %-8s %.3f  (MC SE %.3f)\n", nm, cv, mc_se(cov_par[nm], nok)))
}
cat("\n=== 95%% band coverage for f0(t) by age ===\n")
for (k in seq_along(agrid))
  cat(sprintf("  age %4.1f : %.3f  (MC SE %.3f)\n",
              agrid[k], cov_f0[k]/nok, mc_se(cov_f0[k], nok)))
cat(sprintf("\nOverall f0(t) coverage (pooled): %.3f\n", sum(cov_f0)/(nok*length(agrid))))

saveRDS(list(R=R, N=N, NSUB=NSUB, nok=nok, nbad=nbad, agrid=agrid,
             cov_par=cov_par, cov_f0=cov_f0, pnames=pnames),
        "output/calibration.rds")
cat("Saved -> output/calibration.rds\n")
