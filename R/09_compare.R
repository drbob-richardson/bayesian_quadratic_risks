# =====================================================================
# 09_compare.R -- Benchmark the SPM against standard approaches on the simulated
# (U-shaped) cohort, as an interpretive comparison:
#   (1) Cox PH, marker entered linearly  -> a monotone effect MISSES the U-shape;
#   (2) Cox PH, penalized-spline marker  -> recovers a U-shape but only as a
#       static snapshot of the last observation, with no interpretable optimum,
#       no measurement-error/longitudinal-dynamics handling, and no uncertainty
#       quantification ON the optimum;
#   (3) SPM (this paper)                 -> recovers the optimum f0 with an
#       interval, using the whole trajectory.
# =====================================================================

suppressPackageStartupMessages(library(survival))

sim <- readRDS("output/sim_spm.rds")
mle <- readRDS("output/fit_rtmb.rds")       # SPM MLE fit (N=2000)
p   <- sim$params

# last-observed marker per subject (the usual naive time-fixed covariate)
last <- do.call(rbind, lapply(split(sim$long, sim$long$id),
                              function(d) d[nrow(d), ]))
df <- merge(last[, c("id","y")], sim$surv[, c("id","time","status")], by = "id")

# (1) linear Cox
cox_lin <- coxph(Surv(time, status) ~ y, data = df)
cl <- summary(cox_lin)$coefficients
cat("=== (1) Cox PH, linear marker ===\n")
cat(sprintf("  coef=%.4f  SE=%.4f  p=%.3f   <- monotone: imposes 'more is worse',\n",
            cl[1,"coef"], cl[1,"se(coef)"], cl[1,"Pr(>|z|)"]))
cat("      missing that LOW marker values are also high-risk (the lower arm).\n")

# (2) spline Cox
cox_sp <- coxph(Surv(time, status) ~ pspline(y, df = 4), data = df)
gp <- summary(cox_sp)$coefficients
cat("\n=== (2) Cox PH, penalized-spline marker ===\n")
print(round(anova(cox_sp), 4))
tp <- termplot(cox_sp, se = TRUE, plot = FALSE)$y
tp$y <- tp$y - min(tp$y)                       # center at its own minimum
spline_opt <- tp$x[which.min(tp$y)]            # spline-implied optimum
cat(sprintf("  spline-implied optimum: y = %.3f  (truth = %.3f)\n",
            spline_opt, p$af0))

# (3) SPM optimum (constant f0 here; estimate + 95% CI from sdreport)
spm_opt <- mle$est$af0; spm_se <- mle$se$af0
cat("\n=== (3) SPM optimum f0 ===\n")
cat(sprintf("  f0 = %.3f  [%.3f, %.3f]  (truth = %.3f)\n",
            spm_opt, spm_opt-1.96*spm_se, spm_opt+1.96*spm_se, p$af0))
cat("  ...with full uncertainty, using the whole trajectory (not just the last obs).\n")

# ---- Figure: estimated log relative hazard vs marker --------------------------
yg <- seq(min(df$y), max(df$y), length.out = 200)
lin_lh <- cl[1,"coef"] * yg; lin_lh <- lin_lh - min(lin_lh)

png("output/figures/method_comparison.png", width = 1600, height = 1080, res = 200)
par(mar = c(4.6, 4.6, 3.2, 1), cex.lab = 1.15, cex.axis = 1.05, cex.main = 1.15)
plot(tp$x, tp$y, type = "n", xlab = "Marker value (standardized)",
     ylab = "Estimated log relative hazard",
     main = "Recovering U-shaped biomarker risk: standard models vs. the SPM",
     xlim = range(yg), ylim = range(c(tp$y, lin_lh)))
# spline 95% band
polygon(c(tp$x, rev(tp$x)),
        c(tp$y - 1.96*tp$se, rev(tp$y + 1.96*tp$se)),
        col = adjustcolor("steelblue", 0.18), border = NA)
lines(tp$x, tp$y, col = "steelblue", lwd = 3)                 # Cox spline
lines(yg, lin_lh, col = "darkorange", lwd = 3, lty = 2)        # Cox linear
abline(v = p$af0, col = "black", lwd = 2, lty = 3)             # truth optimum
# SPM optimum with CI (drawn near the bottom)
yb <- par("usr")[3] + 0.04*diff(par("usr")[3:4])
segments(spm_opt-1.96*spm_se, yb, spm_opt+1.96*spm_se, yb, col = "firebrick", lwd = 4)
points(spm_opt, yb, pch = 18, col = "firebrick", cex = 1.6)
legend("topleft", bty = "n", cex = 1.0,
       legend = c("Cox PH, linear marker (monotone: misses the lower arm)",
                  "Cox PH, spline marker (shape only; no optimum interval, ignores dynamics)",
                  "SPM optimum f0 with 95% interval",
                  "true optimum"),
       col = c("darkorange", "steelblue", "firebrick", "black"),
       lwd = c(3, 3, 4, 2), lty = c(2, 1, 1, 3), pch = c(NA, NA, 18, NA))
dev.off()
cat("\nSaved -> output/figures/method_comparison.png\n")
