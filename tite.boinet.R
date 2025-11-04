#install.packages("boinet")
library(boinet)
?tite.boinet


eff.sc1<-c(0.05, 0.08, 0.24, 0.42, 0.44)
eff.sc2<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc3<-c(0.05, 0.08, 0.27, 0.28, 0.44)
eff.sc4<-c(0.05, 0.08, 0.27, 0.27, 0.27)
eff.sc5<-c(0.05, 0.08, 0.27, 0.42, 0.44)
eff.sc6<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc7<-c(0.05, 0.08, 0.42, 0.42, 0.37)

eff.scen<- list(eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6, eff.sc7)



for(i in 1:7){
print(boinet::tite.boinet(
  n.dose=5, start.dose=1, size.cohort=3, n.cohort=20,
  toxprob=c(0.01, 0.05, 0.10, 0.15, 0.20) , effprob=eff.scen[[i]],
  phi = 0.25, phi1 = 0.197, phi2 = 0.298,
  delta = 0.6, delta1 = 0.15,
  alpha.T1 = 0.5, alpha.E1 = 0.5, tau.T=28, tau.E=112,
  te.corr = 0.01, gen.event.time = "weibull",
  accrual = 10, gen.enroll.time = "uniform",
  stopping.npts = size.cohort*n.cohort,
  stopping.prob.T = 0.93, stopping.prob.E = 0.71,
  estpt.method = "obs.prob", obd.method = "utility.scoring",
  w1 = 0.33, w2 = 1.09,
  plow.ast = phi1, pupp.ast = phi2,
  qlow.ast = delta1/2, qupp.ast = delta,
  psi00 = 40, psi11 = 60,
  n.sim = 5000, seed.sim = 100))
}

for(m in 1:7){
  median_survival_month<-9
  n.dose=5
  start.dose=1
  size.cohort=3
  n.cohort=20
  toxprob=c(0.01, 0.05, 0.10, 0.15, 0.20)
  effprob=eff.scen[[i]]
  phi = 0.25
  phi1 = 0.197
  phi2 = 0.298
  delta = 0.6 
  delta1 = 0.15
  alpha.T1 = 0.5
  alpha.E1 = 0.5
  tau.T=28
  tau.E=112
  te.corr = 0.01
  gen.event.time = "weibull"
  accrual = 10
  gen.enroll.time = "uniform"
  stopping.npts = size.cohort*n.cohort
  stopping.prob.T = 0.93
  stopping.prob.E = 0.71
  estpt.method = "obs.prob"
  obd.method = "utility.scoring"
  w1 = 0.33
  w2 = 1.09
  plow.ast = phi1
  pupp.ast = phi2
  qlow.ast = delta1/2
  qupp.ast = delta
  psi00 = 40
  psi11 = 60
  n.sim = 5000
  seed.sim = 100
  
  dosen <- 1:n.dose
  dose <- paste("Dose", dosen, sep = "")
  toxp <- data.frame(t(toxprob))
  colnames(toxp) <- dose
  effp <- data.frame(t(effprob))
  colnames(effp) <- dose
  ncoh <- size.cohort
  nesc <- n.cohort
  nmax <- ncoh * nesc
  design.par <- gridoptim(phi = phi, phi1 = phi1, phi2 = phi2, 
                          delta = delta, delta1 = delta1)
  lambda1 <- design.par$lambda1
  lambda2 <- design.par$lambda2
  eta1 <- design.par$eta1
  pr.alpha <- 1
  pr.beta <- 1
  alpha.T1 <- alpha.T1
  alpha.T2 <- 0.5
  alpha.E1 <- alpha.E1
  alpha.E2 <- 0.5
  efftoxp <- list(toxp = toxp, effp = effp)
  ncop <- copula::normalCopula(te.corr, dim = 2, dispstr = "ex")
  mv.ncop <- NULL
  #if (gen.event.time == "weibull") {
  for (i in 1:n.dose) {
    psi.T <- efftoxp$toxp[i][[1]]
    zetta.T1 <- log(log(1 - psi.T)/log(1 - psi.T + 
                                         alpha.T1 * psi.T))/log(1/(1 - alpha.T2))
    zetta.T2 <- tau.T/(-log(1 - psi.T))^(1/zetta.T1)
    psi.E <- efftoxp$effp[i][[1]]
    zetta.E1 <- log(log(1 - psi.E)/log(1 - psi.E + 
                                         alpha.E1 * psi.E))/log(1/(1 - alpha.E2))
    zetta.E2 <- tau.E/(-log(1 - psi.E))^(1/zetta.E1)
    mv.ncop <- append(mv.ncop, copula::mvdc(copula = ncop, 
                                            margins = c("weibull", "weibull"), paramMargins = list(list(shape = zetta.T1, 
                                                                                                        scale = zetta.T2), list(shape = zetta.E1, 
                                                                                                                                scale = zetta.E2))))
  }
  #}else if (gen.event.time == "uniform") {
  #  for (i in 1:n.dose) {
  #    psi.T <- efftoxp$toxp[i][[1]]
  #    psi.E <- efftoxp$effp[i][[1]]
  #    mv.ncop <- append(mv.ncop, copula::mvdc(copula = ncop, 
  #                                            margins = c("unif", "unif"), paramMargins = list(list(min = 0, 
  #                                                                                                  max = tau.T * (1/psi.T)), list(min = 0, max = tau.E * 
  #                                                                                                                                   (1/psi.E)))))
  #  }
  #}
  data.obs.n <- array(0, dim = c(n.sim, n.dose))
  data.dur <- array(0, dim = c(n.sim))
  obd <- array(0, dim = c(n.sim))
  set.seed(seed.sim)
  for (ss in 1:n.sim) {
    obs.n <- numeric(n.dose)
    obs.tox <- numeric(n.dose)
    obs.tox.n <- numeric(n.dose)
    obs.eff <- numeric(n.dose)
    obs.eff.n <- numeric(n.dose)
    pe <- numeric(n.dose)
    pt <- numeric(n.dose)
    t.enter <- NULL
    t.decision <- 0
    tite.df <- NULL
    curdose <- start.dose
    early.stop <- 0
    for (i in 1:nesc) {
      dlab <- paste("Dose", curdose, sep = "")
      obs.n[curdose] <- obs.n[curdose] + ncoh
      for (j in 1:ncoh) {
        if (j == 1) {
          t.enter[j] <- t.decision
        }else {
          # if (gen.enroll.time == "uniform") {
          t.enter[j] <- t.enter[j - 1] + runif(1, 
                                               0, 2 * accrual)
          #  }else if (gen.enroll.time == "exponential") {
          #    t.enter[j] <- t.enter[j - 1] + rexp(1, 
          #                                        1/accrual)
          #  }
        }
      }
      if (i == nesc) {
        t.decision <- t.enter[length(t.enter)] + max(tau.T, 
                                                     tau.E)
      }else {
        #if (gen.enroll.time == "uniform") {
        t.decision <- t.enter[length(t.enter)] + 
          runif(1, 0, 2 * accrual)
        #}else if (gen.enroll.time == "exponential") {
        #  t.decision <- t.enter[length(t.enter)] + 
        #    rexp(1, 1/accrual)
        #}
      }
      time.te <- copula::rMvdc(ncoh, mv.ncop[[curdose]])
      DLT <- as.numeric(time.te[, 1] <= tau.T)
      ORR <- as.numeric(time.te[, 2] <= tau.E)
      new.tite.df<-data.frame(dose = curdose, 
                              enter = t.enter, endtox = t.enter + apply(as.matrix(1:ncoh), 
                                                                        1, function(x) {
                                                                          min(time.te[x, 1], tau.T)
                                                                        }), dlt = DLT, endeff = t.enter + apply(as.matrix(1:ncoh), 
                                                                                                                1, function(x) {
                                                                                                                  min(time.te[x, 2], tau.E)
                                                                                                                }), orr = ORR)
      lambda<- log(2)/((median_survival_month-1)*28)
      pat_cens<-rexp(nrow(new.tite.df), lambda)+new.tite.df$endtox
      new.tite.df$endeff <- pmin(new.tite.df$endeff, pat_cens)
      new.tite.df$orr <- ifelse(new.tite.df$endeff == pat_cens, NA, new.tite.df$orr)
      
      tite.df <- rbind(tite.df, new.tite.df)
      tite.df$endeff <- ifelse(tite.df$dlt == 1, tite.df$endtox, tite.df$endeff)
      tite.df$orr    <- ifelse(tite.df$dlt == 1, NA, tite.df$orr)
      size.cohort=3
      n.cohort=20
      
      tite.curdose <- tite.df[tite.df$dose == curdose, 
      ]
      gamma.T <- as.numeric(tite.curdose$endtox <= 
                              t.decision)
      gamma.E <- as.numeric(tite.curdose$endeff <= 
                              t.decision)
      while (mean(gamma.T * gamma.E) < 0.5) {
        nextdec.T <- tite.curdose[tite.curdose$endtox > 
                                    t.decision, "endtox"]
        nextdec.E <- tite.curdose[tite.curdose$endeff > 
                                    t.decision, "endeff"]
        t.decision <- min(nextdec.T, nextdec.E)
        gamma.T <- as.numeric(tite.curdose$endtox <= 
                                t.decision)
        gamma.E <- as.numeric(tite.curdose$endeff <= 
                                t.decision)
      }
      gamma.all.T <- as.numeric(tite.df$endtox <= t.decision)
      gamma.all.E <- as.numeric(tite.df$endeff <= t.decision)
      for (ds in 1:n.dose) {
        if (sum(tite.df$dose == ds) > 0) {
          compsub.T <- tite.df[(tite.df$endtox <= t.decision) & 
                                 (tite.df$dose == ds), ]
          pendsub.T <- tite.df[(tite.df$endtox > t.decision) & 
                                 (tite.df$dose == ds), ]
          compsub.E <- tite.df[(tite.df$endeff <= t.decision) & 
                                 (tite.df$dose == ds), ]
          pendsub.E <- tite.df[(tite.df$endeff > t.decision) & 
                                 (tite.df$dose == ds), ]
          x.DLT <- sum(compsub.T$dlt)
          n.DLT <- x.DLT + sum(1 - compsub.T$dlt) + 
            sum(t.decision - pendsub.T$enter)/tau.T
          x.ORR <- sum(compsub.E$orr, na.rm=TRUE)
          n.ORR <- x.ORR + sum(1 - compsub.E$orr, na.rm=TRUE) + 
            sum(t.decision - pendsub.E$enter)/tau.E
          obs.tox[ds] <- x.DLT
          obs.tox.n[ds] <- n.DLT
          pt[ds] <- x.DLT/n.DLT
          obs.eff[ds] <- x.ORR
          obs.eff.n[ds] <- n.ORR
          ifelse(is.na(x.ORR/n.ORR)==FALSE, pe[ds] <- x.ORR/n.ORR,pe[ds] <- 0.5)
        }
      }
      if ((pt[curdose] <= lambda1) & (pe[curdose] <= 
                                      eta1)) {
        nxtdose <- curdose + 1
      }else if ((pt[curdose] < lambda2) & (pe[curdose] > 
                                           eta1)) {
        nxtdose <- curdose
      }else if (pt[curdose] >= lambda2) {
        nxtdose <- curdose - 1
      }else if ((pt[curdose] > lambda1) & (pt[curdose] < 
                                           lambda2) & (pe[curdose] <= eta1)) {
        if (curdose == n.dose) {
          three <- c(curdose - 1, curdose)
          maxpe <- max(pe[three])
          maxpe.ds <- dosen[which((pe == maxpe) & (is.element(dosen, 
                                                              three)))]
          if (length(maxpe.ds) == 1) {
            nxtdose <- maxpe.ds
          }else {
            nxtdose <- sample(maxpe.ds, 1)
          }
        }else if (obs.n[curdose + 1] == 0) {
          nxtdose <- curdose + 1
        }else if (curdose == 1) {
          three <- c(curdose, curdose + 1)
          maxpe <- max(pe[three])
          maxpe.ds <- dosen[which((pe == maxpe) & (is.element(dosen, 
                                                              three)))]
          if (length(maxpe.ds) == 1) {
            nxtdose <- maxpe.ds
          }else {
            nxtdose <- sample(maxpe.ds, 1)
          }
        }else {
          three <- c(curdose - 1, curdose, curdose + 
                       1)
          maxpe <- max(pe[three])
          maxpe.ds <- dosen[which((pe == maxpe) & (is.element(dosen, 
                                                              three)))]
          if (length(maxpe.ds) == 1) {
            nxtdose <- maxpe.ds
          }else {
            nxtdose <- sample(maxpe.ds, 1)
          }
        }
      }
      po.shape1 <- pr.alpha + obs.tox
      po.shape2 <- pr.beta + (obs.n - obs.tox)
      tterm <- pbeta(phi, po.shape1, po.shape2)
      po.shape1 <- pr.alpha + obs.eff
      po.shape2 <- pr.beta + (obs.n - obs.eff)
      eterm <- 1 - pbeta(delta1, po.shape1, po.shape2)
      admflg <- !((eterm < (1 - stopping.prob.E)) | 
                    (tterm < (1 - stopping.prob.T)))
      admdose <- dosen[admflg]
      if (sum(admflg) == 0) {
        early.stop <- 1
        break
      }else if (sum(obs.n >= stopping.npts) > 0) {
        break
      }else {
        if (nxtdose == 0) {
          if (admflg[1]) {
            curdose <- 1
          }else {
            early.stop <- 1
            break
          }
        }else if (nxtdose == (n.dose + 1)) {
          curdose <- n.dose
        }else if (is.element(nxtdose, admdose)) {
          curdose <- nxtdose
        }else if (curdose < nxtdose) {
          if (sum(admdose >= nxtdose) != 0) {
            curdose <- min(admdose[admdose >= nxtdose])
          }
        }else if (curdose >= nxtdose) {
          if (sum(admdose <= nxtdose) != 0) {
            curdose <- max(admdose[admdose <= nxtdose])
          }else {
            early.stop <- 1
            break
          }
        }
      }
    }
    data.obs.n[ss, ] <- obs.n
    data.dur[ss] <- t.decision
    evadose <- dosen[obs.n != 0]
    obspt <- obs.tox[evadose]/obs.n[evadose]
    obspe <- obs.eff[evadose]/obs.n[evadose]
    tterm.obd <- numeric(n.dose)
    eterm.obd <- numeric(n.dose)
    for (i in evadose) {
      po.shape1 <- pr.alpha + obs.tox[i]
      po.shape2 <- pr.beta + (obs.n[i] - obs.tox[i])
      tterm.obd[i] <- pbeta(phi, po.shape1, po.shape2)
      po.shape1 <- pr.alpha + obs.eff[i]
      po.shape2 <- pr.beta + (obs.n[i] - obs.eff[i])
      eterm.obd[i] <- 1 - pbeta(delta1, po.shape1, 
                                po.shape2)
    }
    if (early.stop == 1) {
      obd[ss] <- 0
    }else if (length(evadose) == 1) {
      if ((tterm.obd[evadose] >= (1 - stopping.prob.T)) & 
          (eterm.obd[evadose] >= (1 - stopping.prob.E))) {
        obd[ss] <- evadose
      }
    }else if (sum((tterm.obd[evadose] >= (1 - stopping.prob.T)) & 
                  (eterm.obd[evadose] >= (1 - stopping.prob.E))) >= 
              1) {
      estpt <- Iso::pava(obspt)
      if (estpt.method == "multi.iso") {
        estpe <- multi.iso(obs = obs.eff[evadose], 
                           n = obs.n[evadose])
      }else if (estpt.method == "fp.logistic") {
        estpe <- fp.logit(obs = obs.eff[evadose], n = obs.n[evadose], 
                          dose = evadose)
      }else if (estpt.method == "obs.prob") {
        estpe <- obspe
      }
      obd[ss] <- obd.select(probt = estpt, probe = estpe, 
                            method = obd.method, phi = phi, phi1 = phi1, 
                            phi2 = phi2, delta = delta, delta1 = delta1, 
                            tterm = tterm.obd[evadose], eterm = eterm.obd[evadose], 
                            stopT = stopping.prob.T, stopE = stopping.prob.E, 
                            w1 = w1, w2 = w2, plow.ast = plow.ast, pupp.ast = pupp.ast, 
                            qlow.ast = qlow.ast, qupp.ast = qupp.ast, psi00 = psi00, 
                            psi11 = psi11)
    }
  }
  prop.select <- array(0, dim = c(n.dose))
  for (i in 1:n.dose) {
    prop.select[i] <- round(mean(obd == i) * 100, digits = 1)
  }
  names(prop.select) <- dose
  prop.stop <- round(mean(obd == 0) * 100, digits = 1)
  names(prop.stop) <- "No OBD %"
  n.patient <- round(apply(data.obs.n, 2, mean), digits = 1)
  names(n.patient) <- dose
  duration <- round(mean(data.dur), digits = 1)
  names(duration) <- "Trial duration (days)"
  names(toxprob) <- dose
  names(effprob) <- dose
  names(phi) <- "Target toxicity prob."
  names(delta) <- "Target efficacy prob."
  names(lambda1) <- "Lower toxicity boundary"
  names(lambda2) <- "Upper toxicity boundary"
  names(eta1) <- "Lower efficacy boundary"
  names(tau.T) <- "Tox. assessment window (days)"
  names(tau.E) <- "Eff. assessment window (days)"
  names(accrual) <- "Accrual rate (days)"
  names(estpt.method) <- "Efficacy prob. estimation"
  names(obd.method) <- "OBD selection"
  result <- list(toxprob = toxprob, effprob = effprob, 
                 phi = phi, delta = delta, lambda1 = lambda1, lambda2 = lambda2, 
                 eta1 = eta1, tau.T = tau.T, tau.E = tau.E, accrual = accrual, 
                 estpt.method = estpt.method, obd.method = obd.method, 
                 n.patient = n.patient, prop.select = prop.select, 
                 prop.stop = prop.stop, duration = duration)
  class(result) <- "tite.boinet"
  result
}