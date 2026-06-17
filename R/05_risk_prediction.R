# =====================================================================
# 05_risk_prediction.R  -- Deliverable 2: posterior-predictive RISK / SURVIVAL
# with credible bands, computed from the existing tmbstan posterior draws
# (no re-fit). Two products:
#
#  (1) Population survival by deviation-from-optimal: S(t) for an individual
#      whose biomarker is maintained delta SD from the optimal f0(t), plus the
#      "natural" allostatic-mean path. Shows the survival value of optimal
#      health WITH uncertainty -- the actuarial pay-off of the quadratic hazard.
#
#  (2) Dynamic individual prediction: filter the latent state on a person's
#      observed history, then project conditional survival S(t | history) with
#      a credible band, vs their actual outcome.
# =====================================================================

bay <- readRDS("output/fit_tmbstan.rds")
sim <- readRDS("output/sim_spm.rds")
p   <- sim$params; tmin <- p$tmin
D   <- bay$draws                       # posterior draws, natural scale
nd  <- length(D$af0)
cat(sprintf("Using %d posterior draws.\n", nd))

# component functions evaluated as vectors over draws at a scalar age t
mu0_v <- function(t) D$amu0 * exp(D$bmu0 * (t - tmin))
Q_v   <- function(t) D$aQ   + D$bQ  * (t - tmin)
f0_v  <- function(t) D$af0  + D$bf0 * (t - tmin)
f1_v  <- function(t) D$af1  + D$bf1 * (t - tmin)
a_v   <- function(t) D$aY   + D$bY  * (t - tmin)

# ------------------------------------------------------------------ (1) -------
ages  <- seq(65, 95, by = 0.5); ng <- length(ages)
dt    <- diff(ages)[1]
# deviation scenarios (SD units from optimal) + the allostatic-mean path
devs  <- c(0, 0.5, 1.0, 1.5, 2.0)
scen_names <- c(sprintf("optimal (dev 0)"),
                sprintf("dev +0.5 SD"), sprintf("dev +1 SD"),
                sprintf("dev +1.5 SD"), sprintf("dev +2 SD"))
cols  <- c("forestgreen","gold3","darkorange","orangered","firebrick")

# survival matrices (ndraws x ngrid) per scenario via cumulative trapezoid
surv_scn <- lapply(devs, function(delta) {
  H <- matrix(0, nd, ng)                       # cumulative hazard
  hz_prev <- mu0_v(ages[1]) + Q_v(ages[1]) * delta^2
  for (k in 2:ng) {
    hz <- mu0_v(ages[k]) + Q_v(ages[k]) * delta^2
    H[, k] <- H[, k-1] + (hz_prev + hz)/2 * dt
    hz_prev <- hz
  }
  exp(-H)
})
# allostatic-mean path: deviation = f1(t) - f0(t) (the allostatic load), grows w/ age
H <- matrix(0, nd, ng)
d_prev <- f1_v(ages[1]) - f0_v(ages[1])
hz_prev <- mu0_v(ages[1]) + Q_v(ages[1]) * d_prev^2
for (k in 2:ng) {
  dd <- f1_v(ages[k]) - f0_v(ages[k])
  hz <- mu0_v(ages[k]) + Q_v(ages[k]) * dd^2
  H[, k] <- H[, k-1] + (hz_prev + hz)/2 * dt
  hz_prev <- hz
}
surv_mean <- exp(-H)

qband <- function(M) apply(M, 2, quantile, c(.5,.025,.975))

png("output/figures/survival_by_deviation.png", width = 1600, height = 1100, res = 200)
par(mar = c(4.5, 4.5, 3.5, 1))
plot(ages, qband(surv_scn[[1]])[1,], type="n", ylim=c(0,1),
     xlab="Age (years)", ylab="Survival  S(t | alive at 65)",
     main="Posterior-predictive survival by deviation from optimal biomarker")
for (i in seq_along(devs)) {
  qb <- qband(surv_scn[[i]])
  polygon(c(ages, rev(ages)), c(qb[2,], rev(qb[3,])),
          col=adjustcolor(cols[i],0.18), border=NA)
  lines(ages, qb[1,], col=cols[i], lwd=2.5)
}
qm <- qband(surv_mean)
lines(ages, qm[1,], col="grey30", lwd=2.5, lty=2)
legend("bottomleft", bty="n", lwd=2.5,
       col=c(cols,"grey30"), lty=c(rep(1,5),2),
       legend=c(scen_names, "allostatic-mean path"))
dev.off()
cat("Saved -> output/figures/survival_by_deviation.png\n")

# actuarial summary: P(survive 65 -> 85) by scenario, median [95% CrI]
i85 <- which.min(abs(ages - 85))
cat("\n=== P(survive to 85 | alive at 65), median [95% CrI] ===\n")
for (i in seq_along(devs)) {
  q <- quantile(surv_scn[[i]][, i85], c(.5,.025,.975))
  cat(sprintf("  %-16s %5.3f [%5.3f, %5.3f]\n", scen_names[i], q[1], q[2], q[3]))
}
qm85 <- quantile(surv_mean[, i85], c(.5,.025,.975))
cat(sprintf("  %-16s %5.3f [%5.3f, %5.3f]\n", "allostatic-mean", qm85[1], qm85[2], qm85[3]))

# ------------------------------------------------------------------ (2) -------
# plain-R RK4 moment step for ONE draw (scalars), mirroring spm_nll.R
step1 <- function(m, g, t, dtv, par) {
  af <- function(t) par$aY + par$bY*(t-tmin)
  f1f<- function(t) par$af1 + par$bf1*(t-tmin)
  f0f<- function(t) par$af0 + par$bf0*(t-tmin)
  Qf <- function(t) par$aQ + par$bQ*(t-tmin)
  mf <- function(t) par$amu0*exp(par$bmu0*(t-tmin))
  s1sq <- par$sigma1^2
  dm1 <- af(t)*(m-f1f(t)) - 2*Qf(t)*g*(m-f0f(t)); dg1 <- 2*af(t)*g + s1sq - 2*Qf(t)*g^2
  tm<-t+dtv/2; m2<-m+dm1*dtv/2; g2<-g+dg1*dtv/2
  dm2<-af(tm)*(m2-f1f(tm))-2*Qf(tm)*g2*(m2-f0f(tm)); dg2<-2*af(tm)*g2+s1sq-2*Qf(tm)*g2^2
  m3<-m+dm2*dtv/2; g3<-g+dg2*dtv/2
  dm3<-af(tm)*(m3-f1f(tm))-2*Qf(tm)*g3*(m3-f0f(tm)); dg3<-2*af(tm)*g3+s1sq-2*Qf(tm)*g3^2
  te<-t+dtv; m4<-m+dm3*dtv; g4<-g+dg3*dtv
  dm4<-af(te)*(m4-f1f(te))-2*Qf(te)*g4*(m4-f0f(te)); dg4<-2*af(te)*g4+s1sq-2*Qf(te)*g4^2
  mn<-m+(dm1+2*dm2+2*dm3+dm4)*dtv/6; gn<-g+(dg1+2*dg2+2*dg3+dg4)*dtv/6
  hz <- mf(te) + Qf(te)*((mn-f0f(te))^2 + gn)
  list(m=mn, g=gn, hz=hz)
}
par_d <- function(d) list(aY=D$aY[d],bY=D$bY[d],af1=D$af1[d],bf1=D$bf1[d],
  af0=D$af0[d],bf0=D$bf0[d],aQ=D$aQ[d],bQ=D$bQ[d],amu0=D$amu0[d],bmu0=D$bmu0[d],
  sigma1=D$sigma1[d],sigma0=D$sigma0[d],tau=D$tau[d])

# dynamic conditional survival S(t|history) for one individual, over draws.
# History is TRUNCATED at `landmark` age (landmarking) and survival is projected
# forward from the landmark -- so we predict the future given data-to-date and
# can overlay the actual (later) outcome.
predict_indiv <- function(ind_id, landmark=75, nsub=4, thin=2) {
  lon <- sim$long[sim$long$id==ind_id,]; srv <- sim$surv[sim$surv$id==ind_id,]
  keep <- lon$age <= landmark
  ag <- lon$age[keep]; yy <- lon$y[keep]; tlast <- landmark
  pg <- seq(tlast, 95, by=0.5); npg <- length(pg)
  use_d <- seq(1, nd, by=thin)
  Smat <- matrix(NA, length(use_d), npg)
  for (r in seq_along(use_d)) {
    par <- par_d(use_d[r]); tausq <- par$tau^2
    m <- (par$af1 + par$bf1*(ag[1]-tmin)); g <- par$sigma0^2  # prior at first obs
    for (j in seq_along(ag)) {                               # filter on history
      if (j>1) { t<-ag[j-1]; dtv<-(ag[j]-ag[j-1])/nsub
        for (s in 1:nsub){ st<-step1(m,g,t,dtv,par); m<-st$m; g<-st$g; t<-t+dtv } }
      K <- g/(g+tausq); m <- m + K*(yy[j]-m); g <- (1-K)*g
    }
    # advance filtered state from last observed age up to the landmark
    last_age <- ag[length(ag)]
    if (landmark > last_age) { t<-last_age; dtv<-(landmark-last_age)/nsub
      for (s in 1:nsub){ st<-step1(m,g,t,dtv,par); m<-st$m; g<-st$g; t<-t+dtv } }
    # project forward conditional survival from the landmark
    H <- 0; Smat[r,1] <- 1; mm<-m; gg<-g
    for (k in 2:npg) { t<-pg[k-1]; dtv<-(pg[k]-pg[k-1])/nsub
      for (s in 1:nsub){ st<-step1(mm,gg,t,dtv,par); H<-H+st$hz*dtv; mm<-st$m; gg<-st$g; t<-t+dtv }
      Smat[r,k] <- exp(-H)
    }
  }
  list(pg=pg, qb=apply(Smat,2,quantile,c(.5,.025,.975)),
       ag=ag, yy=yy, tlast=tlast, time=srv$time, status=srv$status)
}

# choose individuals observed through the landmark (>=3 obs by age 75) whose
# outcome occurs well AFTER the landmark, so the forecast has a real horizon.
LANDMARK <- 75
cand <- bay$use_ids
info <- do.call(rbind, lapply(cand, function(id){
  l<-sim$long[sim$long$id==id,]; s<-sim$surv[sim$surv$id==id,]
  nb_pre <- sum(l$age <= LANDMARK)
  data.frame(id=id, nobs_pre=nb_pre, status=s$status, time=s$time,
             meandev=mean(abs(l$y[l$age<=LANDMARK])))}))
elig <- info[info$nobs_pre>=3 & info$time>=LANDMARK+5, ]
ev <- elig[elig$status==1, ]; ce <- elig[elig$status==0, ]
ev_id <- ev$id[which.max(ev$meandev)]   # event case, more deviated history
ce_id <- ce$id[which.min(ce$meandev)]   # censored case, near-optimal history
cat(sprintf("\nDynamic prediction (landmark age %d): event id=%d, censored id=%d\n",
            LANDMARK, ev_id, ce_id))

pe <- predict_indiv(ev_id, landmark=LANDMARK)
pc <- predict_indiv(ce_id, landmark=LANDMARK)
png("output/figures/dynamic_prediction.png", width = 1700, height = 850, res = 200)
par(mfrow=c(1,2), mar=c(4.3,4.3,3,1))
draw_one <- function(pp, ttl) {
  plot(pp$pg, pp$qb[1,], type="n", ylim=c(0,1), xlim=c(min(pp$ag),95),
       xlab="Age (years)", ylab="Conditional survival  S(t | history)", main=ttl)
  polygon(c(pp$pg,rev(pp$pg)), c(pp$qb[2,],rev(pp$qb[3,])),
          col=adjustcolor("steelblue",0.25), border=NA)
  lines(pp$pg, pp$qb[1,], col="steelblue", lwd=2.5)
  abline(v=pp$tlast, col="grey60", lty=3)            # last observation
  rug(pp$ag, col="grey40")                            # observation ages
  if (pp$status==1) { abline(v=pp$time, col="firebrick", lwd=2)
    text(pp$time, 0.05, "event", col="firebrick", pos=4, cex=.8) }
  else { abline(v=pp$time, col="darkgreen", lwd=2, lty=2)
    text(pp$time, 0.05, "censored", col="darkgreen", pos=2, cex=.8) }
}
draw_one(pe, sprintf("Individual %d (had event)", ev_id))
draw_one(pc, sprintf("Individual %d (censored)", ce_id))
dev.off()
cat("Saved -> output/figures/dynamic_prediction.png\n")

saveRDS(list(ages=ages, surv_scn=surv_scn, surv_mean=surv_mean,
             devs=devs, pred_event=pe, pred_cens=pc),
        "output/risk_prediction.rds")
cat("Saved -> output/risk_prediction.rds\n")
