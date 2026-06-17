# =====================================================================
# 08_realdata_nafld.R -- Headline real-data example with a GENUINE U/J-shaped
# marker: systolic blood pressure (SBP) vs mortality in the NAFLD cohort
# (survival::nafld1/nafld2). Unlike PBC (monotone bilirubin), SBP has an interior
# risk-minimizing optimum -- the actuarially compelling case for the model.
# Public and fully reproducible (ships with the survival package).
#
# Marker = standardized SBP; event = death; time axis = age.
# =====================================================================

suppressPackageStartupMessages({ library(survival); library(RTMB) })
source("R/spm_nll.R")

set.seed(2026)

# ---- Build SPM inputs from NAFLD -------------------------------------------
sbp <- nafld2[nafld2$test == "sbp", c("id","days","value")]
sbp <- sbp[is.finite(sbp$value) & sbp$value > 40 & sbp$value < 260, ]   # plausibility
base <- nafld1[, c("id","age","futime","status")]
sbp  <- merge(sbp, base, by = "id")
sbp$age_meas <- sbp$age + sbp$days / 365.25
sbp$tend     <- sbp$age + sbp$futime / 365.25

# standardize SBP (mmHg) -> z; record mean/sd to back-transform f0
mu_y <- mean(sbp$value); sd_y <- sd(sbp$value)
sbp$y <- (sbp$value - mu_y) / sd_y

long <- sbp[sbp$age_meas < sbp$tend, c("id","age_meas","y")]
names(long) <- c("id","age","y")
# SBP is noisy visit-to-visit; average within each integer-year window per
# subject to suppress measurement error and sharpen the curvature signal.
long$ageyr <- floor(long$age)
long <- aggregate(cbind(age, y) ~ id + ageyr, data = long, FUN = mean)
long <- long[, c("id","age","y")]
long <- long[order(long$id, long$age), ]
nob  <- table(long$id); keep_ids <- as.integer(names(nob)[nob >= 3])    # >=3 yearly SBPs
# subsample subjects for tractable taping; keep all their measurements
if (length(keep_ids) > 3000) keep_ids <- sort(sample(keep_ids, 3000))
long <- long[long$id %in% keep_ids, ]
bb   <- base[base$id %in% keep_ids, ]
surv <- data.frame(id = bb$id, time = bb$age + bb$futime/365.25,
                   status = as.integer(bb$status == 1))
sim_like <- list(long = long, surv = surv)

tmin <- floor(min(long$age)); arange <- range(c(long$age, surv$time))
agrid <- seq(floor(arange[1]), ceiling(arange[2]), by = 5)
cat("=== NAFLD / SBP SPM inputs ===\n")
cat(sprintf("Subjects: %d | deaths: %d (%.0f%%) | SBP obs: %d (median %.0f/subj)\n",
            nrow(surv), sum(surv$status), 100*mean(surv$status),
            nrow(long), median(table(long$id))))
cat(sprintf("Age range: %.1f - %.1f | tmin=%d | SBP mean=%.1f sd=%.1f mmHg\n",
            arange[1], arange[2], tmin, mu_y, sd_y))

dat <- pack_data(sim_like, tmin = tmin, nsub = 4L, use_priors = 1L,
                 ids = sort(keep_ids), agrid = agrid)

# ---- MAP fit ----------------------------------------------------------------
init <- spm_init()
init$af0 <- 0.0; init$bf0 <- 0          # interior optimum expected
init$laQ <- log(0.08); init$laY <- log(0.10)
init$lamu0 <- log(0.004); init$bmu0 <- 0.07
init$lsig1 <- log(0.40); init$lsig0 <- log(0.70); init$ltau <- log(0.45)
init$bQ <- 0; init$bf0 <- 0   # age-constant curvature AND age-constant optimum

# The age-trends of the curvature (bQ) and of the optimum (bf0) are confounded and
# not robustly identified from the SBP data: freeing them flips the sign of the
# optimum's slope. We therefore report the robust, confound-free specification with
# an age-constant curvature and an age-constant optimum (a single risk-minimizing
# SBP), fixing bQ = bf0 = 0. All remaining parameters are then identified.
obj <- MakeADFun(make_nll(dat), init,
                 map = list(bQ = factor(NA), bf0 = factor(NA)), silent = TRUE)
lo <- setNames(rep(-Inf, length(obj$par)), names(obj$par))
hi <- setNames(rep( Inf, length(obj$par)), names(obj$par))
lo["laY"]<-log(0.02); hi["laY"]<-log(3);  lo["laQ"]<-log(0.02); hi["laQ"]<-log(3)
lo["lamu0"]<-log(1e-4); hi["lamu0"]<-log(0.5); lo["bmu0"]<--0.05; hi["bmu0"]<-0.30
lo["lsig1"]<-log(0.05); hi["lsig1"]<-log(3); lo["lsig0"]<-log(0.05); hi["lsig0"]<-log(5)
lo["ltau"]<-log(0.02);  hi["ltau"]<-log(3)
cat("\nOptimizing (MAP)...\n")
t0 <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lo[names(obj$par)],
              upper = hi[names(obj$par)], control = list(iter.max=800, eval.max=800))
cat(sprintf("Done in %.1f s | code %d | npost=%.2f\n",
            as.numeric(difftime(Sys.time(), t0, units="secs")),
            opt$convergence, opt$objective))
sdr <- sdreport(obj); e <- as.list(sdr,"Estimate"); s <- as.list(sdr,"Std. Error")

show <- function(lab, est, se) cat(sprintf("  %-8s %9.4f  (SE %8.4f)\n", lab, est, se))
cat("\n=== MAP estimates (standardized SBP scale) ===\n")
show("a(tmin)", -exp(e$laY), exp(e$laY)*s$laY); show("bY", e$bY, s$bY)
show("af1", e$af1, s$af1); show("bf1", e$bf1, s$bf1)
show("af0", e$af0, s$af0); cat("  bf0      (fixed at 0; age-constant optimum)\n")
cat(sprintf("  --> optimal SBP = %.1f mmHg  [95%% CI %.1f, %.1f]\n",
            e$af0*sd_y+mu_y, (e$af0-1.96*s$af0)*sd_y+mu_y, (e$af0+1.96*s$af0)*sd_y+mu_y))
show("aQ", exp(e$laQ), exp(e$laQ)*s$laQ); cat("  bQ       (fixed at 0; age-constant curvature)\n")
show("amu0", exp(e$lamu0), exp(e$lamu0)*s$lamu0); show("bmu0", e$bmu0, s$bmu0)
show("sigma1", exp(e$lsig1), exp(e$lsig1)*s$lsig1); show("tau", exp(e$ltau), exp(e$ltau)*s$ltau)

# f0(t) -> SBP mmHg
rn <- rownames(summary(sdr,"report")); sm <- summary(sdr,"report")[grep("f0_grid",rn),]
to_sbp <- function(z) z*sd_y + mu_y
cat("\n=== Optimal (risk-minimizing) SBP f0(t) ===\n")
cat("  age   f0(std) [95% SE band]      -> SBP mmHg\n")
for (k in seq_along(agrid))
  cat(sprintf("  %4.0f  %6.3f [%6.3f,%6.3f]   %5.1f\n",
              agrid[k], sm[k,1], sm[k,1]-1.96*sm[k,2], sm[k,1]+1.96*sm[k,2],
              to_sbp(sm[k,1])))

# ---- Figure: U-shaped hazard + optimal SBP ----------------------------------
af0 <- e$af0; af0_lo <- e$af0 - 1.96*s$af0; af0_hi <- e$af0 + 1.96*s$af0
aQ  <- exp(e$laQ); tref <- 60
mu0r <- exp(e$lamu0 + e$bmu0*(tref - tmin))
sbp_grid <- seq(95, 200, length.out = 200)
zg <- (sbp_grid - mu_y)/sd_y
haz <- mu0r + aQ*(zg - af0)^2

png("output/figures/nafld_f0_sbp.png", width = 1900, height = 850, res = 200)
par(mfrow = c(1,2), mar = c(4.5, 4.5, 3.2, 1))
# (A) estimated U-shaped hazard as a function of SBP, with the optimum + CI
plot(sbp_grid, haz, type="l", lwd=3, col="firebrick",
     xlab="Systolic blood pressure (mmHg)", ylab=sprintf("Estimated hazard at age %d", tref),
     main="(A) Estimated U-shaped risk and optimal SBP")
abline(v = to_sbp(af0), col="grey30", lwd=2)
rect(to_sbp(af0_lo), -1, to_sbp(af0_hi), 1e9, col=adjustcolor("grey50",0.2), border=NA)
text(to_sbp(af0), max(haz)*0.92, sprintf("optimum %.0f mmHg\n[%.0f, %.0f]",
     to_sbp(af0), to_sbp(af0_lo), to_sbp(af0_hi)), pos=4, cex=0.8)
# (B) observed SBP over age with the (age-constant) optimal level band
plot(long$age, to_sbp(long$y), pch=16, col=adjustcolor("grey60",0.12), cex=0.4,
     xlab="Age (years)", ylab="Systolic blood pressure (mmHg)",
     main="(B) Observed SBP and risk-minimizing level")
rect(-1, to_sbp(af0_lo), 200, to_sbp(af0_hi), col=adjustcolor("firebrick",0.25), border=NA)
abline(h = to_sbp(af0), col="firebrick", lwd=3)
legend("topright", bty="n", legend=c("observed SBP","optimal SBP","95% band"),
       col=c("grey60","firebrick",adjustcolor("firebrick",0.5)), pch=c(16,NA,NA), lwd=c(NA,3,8))
dev.off()
cat("\nSaved -> output/figures/nafld_f0_sbp.png\n")

saveRDS(list(opt=opt, sdr=sdr, agrid=agrid, mu_y=mu_y, sd_y=sd_y,
             sim_like=sim_like), "output/fit_nafld_mle.rds")
cat("Saved -> output/fit_nafld_mle.rds\n")
