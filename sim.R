setwd("/annotated_code")
library(BOIN)
source("functions.R")

### Choosing values of a and b 
shape <- read.csv("shape_param_inc.csv")[,-1]
rate <- read.csv("rate_param_inc.csv")[1,2]


pro.sc1.2<- shape
pro.sc3<- shape[c(1:3,3,3),]
pro.sc4<- shape[c(1:3,5,5),]
pro.sc5.6<-shape[,c(1,4,5,6,7,8,9,9,9)]

pro.scen<- list(pro.sc1.2,pro.sc1.2, pro.sc3, pro.sc4, pro.sc5.6, pro.sc5.6)

eff.sc1<-c(0.05, 0.08, 0.24, 0.42, 0.44)
eff.sc2<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc3<-c(0.05, 0.08, 0.24, 0.28, 0.44)
eff.sc4<-c(0.05, 0.08, 0.24, 0.24, 0.24)
eff.sc5<-c(0.05, 0.08, 0.24, 0.42, 0.44)
eff.sc6<-c(0.05, 0.08, 0.42, 0.42, 0.42)

eff.scen<- list(eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6)

dlt.5<- c(0.01, 0.05, 0.10, 0.15, 0.20)
dlt.3<- c(0.06, 0.13, 0.25, 0.40, 0.50)

mtd<- 5

for(m in 1:6){
  scenario<-m


general_list<- list(n.patient.cohort=3, first.week.assessed=4, between.cohort.wk=4, n.doses=5,
                    n.timepoints=8, final.assessment.timepoint=16, mcmc.niter=1250,
                    mcmc.burnin.prop=0.20, n.cohorts.all=20,  cop.corr=0.9, median_survival_month=9)

boin_list<- list(cdlt_rates=eval(parse(text=paste0("dlt.",mtd))), esc_bound=0.197,des_bound=0.298,target=0.25, beta_a_safety=0.5,beta_b_safety=0.5 )

pro_list<- list(pro.schedule=2,beta_shape_sc=pro.scen[[scenario]], beta_rate_sc=rate,
                pro.between.cohort.timepoints=2)

eff_list<- list(eff.schedule=c(8,16), eff_rates=eff.scen[[scenario]], beta_a=0.1, beta_b=0.9, min_eff=0.1, interim_complete_cohort1=6,
                interim_complete_cohort2=12)
data.after.dlt<- TRUE
n.sim<-1000
cl <- makeCluster(detectCores())
clusterSetRNGStream(cl,1915)
invisible(clusterEvalQ(cl,{
  library(rjags)
  library(tidyverse)
  library(clusterGeneration)
  library(parallel)
  library(mvtnorm)
  setwd("/annotated_code")
  source("functions.R")
}))
clusterExport(cl, c("general_list", "boin_list", "pro_list", "eff_list", "n.sim", "data.after.dlt"))
val<-parLapply(cl, 1:n.sim, function (k) trial_design(general_list, boin_list, pro_list, eff_list, data.after.dlt))
stopCluster(cl)

final.sim<-sum(sapply(1:n.sim, function (k) sum(is.na(val[[k]][[1]]))==FALSE))
final.rec<- rep(0, times =5)
efficacy.estimate<- matrix(NA, nrow=final.sim, ncol=5)
pro.estimate<- matrix(NA, nrow=final.sim, ncol=5)
n.patient.allocated<- matrix(NA, nrow=final.sim, ncol=5)
admiss<- matrix(NA, nrow=final.sim, ncol=5)
loss<- matrix(NA, nrow=final.sim, ncol=5)
ncens<- matrix(NA, nrow=final.sim, ncol=5)
ndlt<- matrix(NA, nrow=final.sim, ncol=5)


j<-1
for(i in 1:n.sim){
  if(sum(is.na(val[[i]][[1]]))==FALSE){
    #final.rec[as.numeric(val[[i]][[1]])]<-final.rec[as.numeric(val[[i]][[1]])]+1
    ncens[j,val[[i]][[2]]]<- val[[i]][[10]]
    ndlt[j,val[[i]][[2]]]<- val[[i]][[11]]
    efficacy.estimate[j,val[[i]][[2]]]<-val[[i]][[8]]
    pro.estimate[j,val[[i]][[2]]]<-val[[i]][[9]]
    loss[j,val[[i]][[2]]]<- val[[i]][[7]]
    n.patient.allocated[j,val[[i]][[2]]]<-val[[i]][[3]]
    admiss[j,val[[i]][[1]]]<- 1
    j<-j+1
  }
}

j<-1
mtd.rec<- c()
mtd.prob<- matrix(ncol=5, nrow=final.sim)
for(i in 1:n.sim){
  if(sum(is.na(val[[i]][[1]]))==FALSE){
    npat<-sapply(1:5, function (k) sum(val[[i]][[4]][,2]==k))
    ntoxi<-sapply(1:5, function (k) sum(val[[i]][[4]][val[[i]][[4]][,2]==k,4]==1))
    boin<-select.mtd(0.25, npat, ntoxi, cutoff.eli=0.95, extrasafe=FALSE, offset=0.05,
            boundMTD=FALSE,p.tox=1.4*0.25)
    mtd.rec[j]<- boin$MTD
    mtd.prob[j,]<- as.numeric(boin$p_est[,2])
    j<-j+1
  }
}

if(data.after.dlt==FALSE){
  #write.csv(final.rec/n.sim, paste0("/results/sens.sc",scenario,".", mtd,".60.r.csv"))
  write.csv(efficacy.estimate, paste0("/sens.sc",scenario,".", mtd,".60.e.csv"))
  write.csv(pro.estimate, paste0("/sens.sc",scenario,".", mtd,".60.p.csv"))
  write.csv(loss, paste0("/sens.sc",scenario,".", mtd,".60.l.csv"))
  write.csv(n.patient.allocated, paste0("/sens.sc",scenario,".", mtd,".60.n.csv"))
  write.csv(admiss, paste0("/sens.sc",scenario,".", mtd,".60.a.csv"))
  write.csv(mtd.rec, paste0("/sens.sc",scenario,".", mtd,".60.mr.csv"))
  write.csv(mtd.prob, paste0("/sens.sc",scenario,".", mtd,".60.mp.csv"))
  write.csv(ncens, paste0("/sens.sc",scenario,".", mtd,".60.ncens.csv"))
  write.csv(ndlt, paste0("/sens.sc",scenario,".", mtd,".60.ndlt.csv"))
}else{
#write.csv(final.rec/n.sim, paste0("sc",scenario,".", mtd,".60.r.csv"))
write.csv(efficacy.estimate, paste0("/sc",scenario,".", mtd,".60.e.csv"))
write.csv(pro.estimate, paste0("/sc",scenario,".", mtd,".60.p.csv"))
write.csv(loss, paste0("/sc",scenario,".", mtd,".60.l.csv"))
write.csv(n.patient.allocated, paste0("/sc",scenario,".", mtd,".60.n.csv"))
write.csv(admiss, paste0("/sc",scenario,".", mtd,".60.a.csv"))
write.csv(mtd.rec, paste0("/sc",scenario,".", mtd,".60.mr.csv"))
write.csv(mtd.prob, paste0("/sc",scenario,".", mtd,".60.mp.csv"))
write.csv(ncens, paste0("/sc",scenario,".", mtd,".60.ncens.csv"))
write.csv(ndlt, paste0("/sc",scenario,".", mtd,".60.ndlt.csv"))
}
}
