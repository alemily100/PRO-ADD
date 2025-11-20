#!/usr/bin/env Rscript


setwd("/home/ealger/pro_add")

args <- commandArgs(trailingOnly=TRUE)
sc<-args[1]



#setwd("M:/PhD/Trial Designs/mayo_clinic")

library(BOIN)
source("functions.R")


### Choosing values of a and b 
#shape <- read.csv("M:/PhD/Trial Designs/mayo_clinic/shape_param_inc.csv")[,-1]
#rate <- read.csv("M:/PhD/Trial Designs/mayo_clinic/rate_param_inc.csv")[1,2]



shape <- read.csv("/home/ealger/pro_add/shape_param_inc.csv")[,-1]
rate <- read.csv("/home/ealger/pro_add/rate_param_inc.csv")[1,2]

pro.sc1.2<- shape
pro.sc3<- shape[c(1:3,3,3),]
pro.sc4<- shape[c(1:3,5,5),]
pro.sc5.6<-shape[,c(1,4,5,6,7,8,9,9,9)]

pro.scen<- list(pro.sc1.2,pro.sc1.2, pro.sc3, pro.sc4, pro.sc5.6, pro.sc5.6, pro.sc1.2, pro.sc1.2)

eff.sc1<-c(0.05, 0.08, 0.24, 0.42, 0.44)
eff.sc2<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc3<-c(0.05, 0.08, 0.27, 0.28, 0.44)
eff.sc4<-c(0.05, 0.08, 0.27, 0.27, 0.27)
eff.sc5<-c(0.05, 0.08, 0.27, 0.42, 0.44)
eff.sc6<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc7<-c(0.05, 0.08, 0.42, 0.42, 0.37)
eff.sc8<-c(0.27, 0.29, 0.29, 0.29, 0.29)

eff.scen<- list(eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6, eff.sc7, eff.sc8)

dlt.5<- c(0.01, 0.05, 0.10, 0.15, 0.20)
dlt.3<- c(0.06, 0.13, 0.25, 0.40, 0.50)

mtd<- 3

for(m in sc:sc){
  scenario<-m


general_list<- list(n.patient.cohort=3, first.week.assessed=4, between.cohort.wk=4, n.doses=5,
                    n.timepoints=8,n.toxicity=78,n.grades=5, final.assessment.timepoint=16, mcmc.niter=1250,
                    mcmc.burnin.prop=0.20, n.cohorts.all=20, n.sim.final=600, cop.corr=0.9, correlation_dlt_response=0,
                    correlation_r1_r2=0.7)

boin_list<- list(cdlt_rates=eval(parse(text=paste0("dlt.",mtd))), esc_bound=0.197,des_bound=0.298,target=0.25, beta_a_safety=0.5,beta_b_safety=0.5 )

pro_list<- list(pro.schedule=2,beta_shape_sc=pro.scen[[scenario]], beta_rate_sc=rate,
                pro.between.cohort.timepoints=2, max.pro=0.6)

eff_list<- list(eff.schedule=c(8,16), eff_rates=eff.scen[[scenario]], beta_a=0.1, beta_b=0.9, min_eff=0.1, interim_complete_cohort1=6,
                interim_complete_cohort2=12)

n.sim<-5000
cl <- makeCluster(detectCores())
clusterSetRNGStream(cl,1915)
invisible(clusterEvalQ(cl,{
  library(rjags)
  library(tidyverse)
  library(clusterGeneration)
  library(parallel)
  library(mvtnorm)
  setwd("/home/ealger/pro_add")
  #setwd("M:/PhD/Trial Designs/mayo_clinic")
  source("functions.R")
}))
clusterExport(cl, c("general_list", "boin_list", "pro_list", "eff_list", "n.sim"))
val<-parLapply(cl, 1:n.sim, function (k) trial_design_hypothetical_all(general_list, boin_list, pro_list, eff_list))
stopCluster(cl)


final.rec<- rep(0, times =5)
efficacy.estimate<- matrix(NA, nrow=n.sim, ncol=5)
efficacy.estimate.bb<- matrix(NA, nrow=n.sim, ncol=5)
pro.estimate<- matrix(NA, nrow=n.sim, ncol=5)
n.patient.allocated<- matrix(NA, nrow=n.sim, ncol=5)
admiss<- matrix(NA, nrow=n.sim, ncol=5)
loss<- matrix(NA, nrow=n.sim, ncol=5)
#ncens<- matrix(NA, nrow=n.sim, ncol=9)
ndlt<- matrix(NA, nrow=n.sim, ncol=5)
loss.bb<- matrix(NA, nrow=n.sim, ncol=5)
futility1<-matrix(NA, nrow=n.sim, ncol=5)
futility2<-matrix(NA, nrow=n.sim, ncol=5)
futilityfinal<-matrix(NA, nrow=n.sim, ncol=5)
safetyfinal<-matrix(NA, nrow=n.sim, ncol=20)
sample_size<-matrix(NA, nrow=n.sim, ncol=2)

for(i in 1:n.sim){
    #final.rec[as.numeric(val[[i]][[1]])]<-final.rec[as.numeric(val[[i]][[1]])]+1
    ndlt[i,val[[i]][[2]]]<- val[[i]][[10]]
    efficacy.estimate[i,val[[i]][[2]]]<-val[[i]][[8]]
    pro.estimate[i,val[[i]][[2]]]<-val[[i]][[9]]
    loss[i,val[[i]][[2]]]<- val[[i]][[7]]
    n.patient.allocated[i,val[[i]][[2]]]<-val[[i]][[3]]
    efficacy.estimate.bb[i,val[[i]][[2]]]<-val[[i]][[12]]
    loss.bb[i,val[[i]][[2]]]<- val[[i]][[11]]
    admiss[i,val[[i]][[1]]]<- 1
    #ncens[i,]<-val[[i]][[10]]
    ndlt[i,val[[i]][[2]]]<-val[[i]][[10]]
    futility1[i,val[[i]][[2]]]<-(val[[i]][[13]])[val[[i]][[2]]]
    futility2[i,val[[i]][[2]]]<-(val[[i]][[14]])[val[[i]][[2]]]
    futilityfinal[i,val[[i]][[2]]]<-(val[[i]][[15]])[val[[i]][[2]]]
    safetyfinal[i,1:length(val[[i]][[16]])]<-val[[i]][[16]]
    sample_size[i,]<-c(val[[i]][[17]], val[[i]][[18]])
}
mtd.rec<- c()
mtd.prob<- matrix(ncol=5, nrow=n.sim)
for(i in 1:n.sim){
  npat<-sapply(1:5, function (k) sum(val[[i]][[4]][,2]==k))
  ntoxi<-sapply(1:5, function (k) sum(val[[i]][[4]][val[[i]][[4]][,2]==k,4]==1))
  boin<-select.mtd(0.25, npat, ntoxi, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05,
                   boundMTD=FALSE,p.tox=1.4*0.25)
  mtd.rec[i]<- boin$MTD
  mtd.prob[i,]<- as.numeric(boin$p_est[,2])
}
  #write.csv(final.rec/n.sim, paste0("/results/sens.sc",scenario,".", mtd,".60.r.csv"))
  

  write.csv(efficacy.estimate, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.e.csv"))
  write.csv(efficacy.estimate.bb, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.e.bb.csv"))
  write.csv(pro.estimate, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.p.csv"))
  write.csv(loss, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.l.csv"))
  write.csv(loss.bb, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.l.bb.csv"))
  write.csv(n.patient.allocated, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.n.csv"))
  write.csv(admiss, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.a.csv"))
  write.csv(mtd.rec, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.mr.csv"))
  write.csv(mtd.prob, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.mp.csv"))
  #write.csv(ncens, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.ncens.csv"))
  write.csv(ndlt, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.ndlt.csv"))
  write.csv(futility1, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.futility1.csv"))
  write.csv(futility2, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.futility2.csv"))
  write.csv(futilityfinal, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.futilityfinal.csv"))
  write.csv(safetyfinal, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.safetyfinal.csv"))
  write.csv(sample_size, paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.sc",scenario,".", mtd,".60.samplesize.csv"))
}

saveRDS(general_list, file = paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.general_list"))
saveRDS(boin_list, file = paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.boin_list"))
saveRDS(pro_list, file = paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.pro_list"))
saveRDS(eff_list, file = paste0("/home/ealger/pro_add/results/hypothetical.",general_list[[14]],".full/hypothetical.eff_list"))
