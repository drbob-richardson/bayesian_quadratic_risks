# =====================================================================
# 06_realdata_pbc.R -- First real-data shot: fit the SPM to the PBC sequential
# data (survival::pbcseq). Longitudinal marker = standardized log serum
# bilirubin; event = death; time axis = age. This is a MACHINERY test on real,
# irregular, partially-missing longitudinal data + survival -- NOT a U-shape
# story (bilirubin risk is monotonic, so f0 is expected at the low end).
#
# pbcseq: 312 subjects, sequential labs. status 0=censored,1=transplant,2=dead.
# We treat death (status==2) as the event and censor at transplant/last follow-up.
# =====================================================================

suppressPackageStartupMessages({ library(survival); library(RTMB) })
source("R/spm_nll.R")

# ---- Build SPM inputs from pbcseq ------------------------------------------
ps <- pbcseq
ps <- ps[!is.na(ps$bili) & ps$bili > 0, ]            # need the marker
ps$age_meas <- ps$age + ps$day / 365.25              # age at each measurement
ps$logbili  <- log(ps$bili)

# standardize the marker (zero mean, unit variance) as in the paper
mu_y <- mean(ps$logbili); sd_y <- sd(ps$logbili)
ps$y <- (ps$logbili - mu_y) / sd_y

# baseline (one row/subject) for the time-to-event part
base <- ps[!duplicated(ps$id), ]
base$tend   <- base$age + base$futime / 365.25       # age at event/censor
base$event  <- as.integer(base$status == 2)          # death = event

long <- data.frame(id = ps$id, age = ps$age_meas, y = ps$y)
long <- long[order(long$id, long$age), ]
# drop any measurement at/after the event age (numerical safety)
long <- merge(long, base[, c("id","tend")], by = "id")
long <- long[long$age < long$tend, c("id","age","y")]
# keep subjects with >= 2 longitudinal points
nob  <- table(long$id); keep_ids <- as.integer(names(nob)[nob >= 2])
long <- long[long$id %in% keep_ids, ]
surv <- data.frame(id = base$id, time = base$tend, status = base$event)
surv <- surv[surv$id %in% keep_ids, ]

sim_like <- list(long = long, surv = surv)
tmin <- floor(min(long$age))
arange <- range(c(long$age, surv$time))
agrid  <- seq(floor(arange[1]), ceiling(arange[2]), by = 5)

cat("=== PBC SPM inputs ===\n")
cat(sprintf("Subjects: %d | deaths: %d (%.0f%%) | long obs: %d (median %.0f/subj)\n",
            nrow(surv), sum(surv$status), 100*mean(surv$status),
            nrow(long), median(table(long$id))))
cat(sprintf("Age range: %.1f - %.1f  | tmin=%d\n", arange[1], arange[2], tmin))
cat(sprintf("log-bili standardization: mean=%.3f sd=%.3f\n", mu_y, sd_y))

# Use the weakly-informative priors (MAP fit): on real, monotonic-risk data the
# pure MLE collapses the quadratic term (aQ -> 0); the priors keep it identified.
dat <- pack_data(sim_like, tmin = tmin, nsub = 5L, use_priors = 1L,
                 ids = sort(keep_ids), agrid = agrid)

# ---- MAP fit (penalized likelihood) ----------------------------------------
# PBC-scale inits: standardized marker, optimum at the low-bilirubin end,
# non-degenerate quadratic curvature and mean reversion to start.
init <- spm_init()
init$af0 <- -0.6; init$bf0 <- 0
init$laQ <- log(0.10); init$laY <- log(0.10)
init$lamu0 <- log(0.003); init$bmu0 <- 0.06
init$lsig1 <- log(0.40); init$lsig0 <- log(0.60); init$ltau <- log(0.30)

obj <- MakeADFun(make_nll(dat), init, silent = TRUE)
# safety bounds (transformed scale) to keep the optimizer in a sane region
lo <- setNames(rep(-Inf, length(obj$par)), names(obj$par))
hi <- setNames(rep( Inf, length(obj$par)), names(obj$par))
lo["laY"]<-log(0.02); hi["laY"]<-log(3)
lo["laQ"]<-log(0.02); hi["laQ"]<-log(3)
lo["lamu0"]<-log(1e-4); hi["lamu0"]<-log(0.5)
lo["bmu0"]<--0.05;    hi["bmu0"]<-0.30
lo["lsig1"]<-log(0.05); hi["lsig1"]<-log(3)
lo["lsig0"]<-log(0.05); hi["lsig0"]<-log(5)
lo["ltau"]<-log(0.02);  hi["ltau"]<-log(3)
cat("\nOptimizing (MAP, with priors + bounds)...\n")
t0 <- Sys.time()
opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lo[names(obj$par)],
              upper = hi[names(obj$par)],
              control = list(iter.max = 800, eval.max = 800))
cat(sprintf("Done in %.1f s | code %d | npost=%.2f\n",
            as.numeric(difftime(Sys.time(), t0, units="secs")),
            opt$convergence, opt$objective))
sdr <- sdreport(obj)
est <- as.list(sdr, "Estimate"); se <- as.list(sdr, "Std. Error")

show <- function(lab, e, s) cat(sprintf("  %-8s %9.4f  (SE %8.4f)\n", lab, e, s))
cat("\n=== MAP estimates (standardized log-bili scale) ===\n")
show("a(tmin)", -exp(est$laY), exp(est$laY)*se$laY)
show("bY",       est$bY,        se$bY)
show("af1",      est$af1,       se$af1)
show("bf1",      est$bf1,       se$bf1)
show("af0",      est$af0,       se$af0)
show("bf0",      est$bf0,       se$bf0)
show("aQ",       exp(est$laQ),  exp(est$laQ)*se$laQ)
show("bQ",       est$bQ,        se$bQ)
show("amu0",     exp(est$lamu0),exp(est$lamu0)*se$lamu0)
show("bmu0",     est$bmu0,      se$bmu0)
show("sigma1",   exp(est$lsig1),exp(est$lsig1)*se$lsig1)
show("tau",      exp(est$ltau), exp(est$ltau)*se$ltau)

# f0(t) with SE bands, mapped back to the bilirubin (mg/dL) scale
rn <- rownames(summary(sdr,"report")); sm <- summary(sdr,"report")[grep("f0_grid",rn),]
f0_std <- sm[,1]; f0_se <- sm[,2]
to_bili <- function(z) exp(z * sd_y + mu_y)          # invert standardization+log
cat("\n=== Optimal (risk-minimizing) marker f0(t) ===\n")
cat("  age   f0(std)  [95% SE band]      -> bilirubin mg/dL (median)\n")
for (k in seq_along(agrid)) {
  cat(sprintf("  %4.0f  %7.3f [%7.3f,%7.3f]   %6.2f\n",
              agrid[k], f0_std[k], f0_std[k]-1.96*f0_se[k], f0_std[k]+1.96*f0_se[k],
              to_bili(f0_std[k])))
}

# ---- Figure: estimated optimal bilirubin f0(t) over the observed data --------
png("output/figures/pbc_f0_bilirubin.png", width = 1600, height = 1050, res = 200)
par(mar = c(4.5, 4.5, 3.5, 1))
obs_bili <- exp(long$y * sd_y + mu_y)                 # back to mg/dL
plot(long$age, obs_bili, log = "y", pch = 16,
     col = adjustcolor("grey60", 0.30), cex = 0.5,
     xlab = "Age (years)", ylab = "Serum bilirubin (mg/dL, log scale)",
     main = "PBC: observed bilirubin and estimated risk-minimizing f0(t)")
f0_bili <- to_bili(f0_std)
f0_lo   <- to_bili(f0_std - 1.96*f0_se)
f0_hi   <- to_bili(f0_std + 1.96*f0_se)
polygon(c(agrid, rev(agrid)), c(f0_lo, rev(f0_hi)),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(agrid, f0_bili, col = "steelblue", lwd = 3)
legend("topright", bty = "n",
       legend = c("observed bilirubin", "estimated optimal f0(t)", "95% SE band"),
       col = c("grey60", "steelblue", adjustcolor("steelblue",0.5)),
       pch = c(16, NA, NA), lwd = c(NA, 3, 8))
dev.off()
cat("Saved -> output/figures/pbc_f0_bilirubin.png\n")

saveRDS(list(opt=opt, sdr=sdr, dat=dat, agrid=agrid, mu_y=mu_y, sd_y=sd_y,
             sim_like=sim_like), "output/fit_pbc_mle.rds")
cat("Saved -> output/fit_pbc_mle.rds\n")
