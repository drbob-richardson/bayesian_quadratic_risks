# =====================================================================
# 04_plot_f0.R -- Visualize the deliverable: posterior credible bands on the
# optimal trajectory f0(t), with the ML SE-bands and the simulation truth.
# Base-R graphics (no extra deps).
# =====================================================================

bay <- readRDS("output/fit_tmbstan.rds")
mle <- readRDS("output/fit_rtmb.rds")
sim <- readRDS("output/sim_spm.rds")
p   <- sim$params; tmin <- p$tmin

agrid   <- bay$agrid
f0_true <- p$af0 + p$bf0 * (agrid - tmin)
bq      <- bay$f0_q                      # 3 x G : median, lo, hi (Bayesian)

# ML SE-band for f0(t) from sdreport ADREPORT
sdr <- mle$sdr
sm  <- summary(sdr, "report"); sm <- sm[grep("f0_grid", rownames(sm)), ]
mle_e <- sm[,1]; mle_lo <- sm[,1]-1.96*sm[,2]; mle_hi <- sm[,1]+1.96*sm[,2]

png("output/figures/f0_credible_bands.png", width = 1600, height = 1100, res = 220)
par(mar = c(4.6, 4.8, 3.4, 1), cex.lab = 1.2, cex.axis = 1.05, cex.main = 1.2,
    mgp = c(2.8, 0.8, 0))
ylim <- range(c(bq, mle_lo, mle_hi, f0_true)) + c(-0.02, 0.02)
plot(agrid, bq[1,], type = "n", ylim = ylim,
     xlab = "Age (years)", ylab = expression("Optimal biomarker level  " * f[0](t)),
     main = expression("Optimal trajectory " * f[0](t) * ": Bayesian band vs. ML vs. truth"))
abline(h = 0, col = "grey85")
# Bayesian 95% credible band
polygon(c(agrid, rev(agrid)), c(bq[2,], rev(bq[3,])),
        col = adjustcolor("steelblue", 0.25), border = NA)
lines(agrid, bq[1,], col = "steelblue", lwd = 2.5)
# ML 95% band (dashed envelope)
lines(agrid, mle_lo, col = "darkorange", lwd = 1.5, lty = 2)
lines(agrid, mle_hi, col = "darkorange", lwd = 1.5, lty = 2)
lines(agrid, mle_e,  col = "darkorange", lwd = 2)
# truth
lines(agrid, f0_true, col = "black", lwd = 2.5, lty = 3)
legend("topleft", bty = "n", cex = 1.02, seg.len = 2.2,
       legend = c("Bayesian posterior median",
                  "Bayesian 95% credible band",
                  expression("ML estimate " %+-% " 1.96 SE"),
                  "Truth"),
       col = c("steelblue", adjustcolor("steelblue",0.5), "darkorange", "black"),
       lwd = c(3, 9, 2.5, 3), lty = c(1, 1, 2, 3))
dev.off()
cat("Saved -> output/figures/f0_credible_bands.png\n")
cat(sprintf("Note: Bayesian fit N=%d (subsample); ML fit N=%d (full). Bands not on identical data.\n",
            bay$config$N, length(unique(sim$surv$id))))
