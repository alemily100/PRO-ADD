######################load packages #################################
library(rjags)
library(tidyverse)
library(clusterGeneration)
library(parallel)
library(fitdistrplus)
library(VineCopula)
library(mvtnorm)

#####################################################

#Description: Simulate initial dose-escalation procedure using BOIN trial design

#Input:
#cdlt.rate: vector of true DLT rates for each dose
#n.cohort: number of patients at each cohort
#week.assessed: the week DLTs are assessed (week 4 in this paper)
#between.cohort.wk: the number of weeks between cohort recruitment (week 4 in this paper)
#n.dose : number of doses
#esc_bound: escalation threshold as per BOIN
#de_esc_bound: de-escalation threshold as per BOIN
#target: target DLT rate
#copula: vector of CDF values to impute correlation between DLT rate and PRO score
#alpha: shape parameter for prior as per stopping rule 
#beta: rate parameter for prior as per stopping rule 
#cohort_allot_interim1: cohort for first evaluation of futility stopping rule 

#Output:list; matrix detailing DLT observations for each patient and vector of dose recommendations  


boin<-function(cdlt.rate, n.cohort,week.assessed, between.cohort.wk, n.dose, esc_bound, de_esc_bound, target, copula, alpha, beta, cohort_allot_interim1){
  i<-current_dose<- 1
  inadmiss<-c()
  acc.dose<- 1:n.dose
  dosage<-c()
  dosage[i]<- current_dose
  M<- matrix(nrow=n.cohort, ncol=5)
  colnames(M)<-c("subj", "dose", "week", "c-dlt", "week_timeline")
  M[,1]<- 1:n.cohort
  M[,2]<- rep(current_dose, times=n.cohort)
  M[,3]<- rep(week.assessed, times=n.cohort)
  M[,4]<- sapply(1:n.cohort, function (k) qbinom(copula[1,k], 1, cdlt.rate[current_dose]))
  M[,5]<- rep(week.assessed, times=n.cohort)
  week<- week.assessed+between.cohort.wk
  stop<- 0
  while((sum(M[,2]==current_dose)<6 |(stop==0)) & (max(M[,1])/n.cohort)!=(cohort_allot_interim1-1)){
    if(1-pbeta(target,alpha+sum(M[M[,2]==current_dose,4]),sum(M[,2]==current_dose)-sum(M[M[,2]==current_dose,4]) + beta)>0.95){
      inadmiss<- unique(c(inadmiss, current_dose:n.dose))
    }
    if(length(inadmiss)==n.dose){
      dosage<- c(dosage,0)
      break
    }
    p.hat<-sum(M[M[,2]==current_dose,4])/sum(M[,2]==current_dose)
    if(p.hat<= esc_bound){
      ifelse(current_dose==n.dose, current_dose<- min(n.dose, min(inadmiss)-1), current_dose<- min(current_dose+1, min(inadmiss)-1))
    }
    if(p.hat>=de_esc_bound){
      ifelse(current_dose==1, current_dose<-1, current_dose<- min(current_dose-1, min(inadmiss)-1))
    }
    if(p.hat> esc_bound && p.hat<de_esc_bound){
      current_dose<- min(current_dose, min(inadmiss)-1)
    }
    if(current_dose==0){
      dosage<- c(dosage,0)
      break
    }
    i<- i+1
    dosage[i]<- current_dose
    new.subj<- ((week/between.cohort.wk-1)*n.cohort+1):((week/between.cohort.wk-1)*n.cohort+n.cohort)
    new.c.dlt<- sapply(new.subj, function (k) qbinom(copula[1,k], 1, cdlt.rate[current_dose]))
    M<-rbind(M, matrix(c(new.subj, rep(current_dose, times=n.cohort), rep(week.assessed, times=n.cohort), new.c.dlt, rep(week, times=n.cohort)), ncol=5))
    week<- week+between.cohort.wk
    stop_test<- sum(M[M[,2]==current_dose,4])/sum(M[,2]==current_dose)
    if(length(inadmiss)==0){
      if((((stop_test> esc_bound) && (stop_test<de_esc_bound)) == TRUE) | ((current_dose==5) && (stop_test<= esc_bound))==TRUE){
        stop<-1
      }
    }else{
      if((((stop_test> esc_bound) && (stop_test<de_esc_bound)) == TRUE) | ((current_dose==max((1:n.dose)[-inadmiss])) && (stop_test<= esc_bound))==TRUE){
        stop<-1
      }
    }
  }
  return(list(M,dosage))
}

#########################################################

#Description: Provide new dose recommendation as per boin following first futility assessment 

#Input:
#target: target DLT rate
#dose: vector of dose recommendations for each patient
#cdlt: vector of cdlt observations for each patient 
#esc_bound: escalation threshold as per BOIN
#n.dose : number of doses
#alpha: shape parameter for prior as per stopping rule 
#beta: rate parameter for prior as per stopping rule 

#Output: vector; set of admissible doses for safety for stage 2 of the design  

boin_decision<- function(target, dose, cdlt, esc_bound, n.dose, alpha, beta){
  max.d<-max(dose)
  admiss<- c()
  if(mean(cdlt[which(dose==max.d)])<=esc_bound & max.d<n.dose){
    ifelse(max.d+1> n.dose, admiss<- max.d, admiss<- max.d+1)
  }
  else{
    admiss<-which(sapply(1:max.d, function (k) 1-pbeta(target,alpha+sum(cdlt[which(dose==k)]),length(which(dose==k))-sum(cdlt[which(dose==k)])+ beta)<=0.95)==TRUE)
  }
  return(admiss)
}

###############################################

#Description: Futility stopping rule decision 

#Input:
#min_target: Minimum acceptable efficacy
#dose: vector of dose recommendations for each patient
#eff: vector of efficacy observations for each patient 
#alpha: shape parameter for prior as per stopping rule 
#beta: rate parameter for prior as per stopping rule 
#n.dose: number of doses

#Output: vector; set of admissible doses for efficacy for stage 2 of the design

futility_decision<- function(min_target, dose, eff, alpha, beta, n.dose){
  not.admiss<-sapply(1:n.dose, function (k) (pbeta(min_target,alpha+sum(eff[which(dose==k)]),length(which(dose==k))-sum(eff[which(dose==k)])+ beta)>0.7)&(sum(dose==k)>=6))
  admiss<- (1:n.dose)[!not.admiss]
  return(admiss)
}

#######################################

#Description: Updated data frame showing DLT observations for next cohort of patients as per the recommended dose 

#Input:
#clin_data: matrix of DLT obersvations to be updated
#n.pat.cohort: number of patients at each cohort
#between.cohort.wk:the number of weeks between cohort recruitment (week 4 in this paper)
#cdlt.rate: vector of true DLT rates for each dose
#rec_dose: dose next cohort of patients allocated to
#copula: vector of CDF values to impute correlation between DLT rate and PRO score

#Output: matrix; containing updated DLT information for all cohorts allocated so far

boin_next<- function(clin_data, n.pat.cohort, between.cohort.wk, cdlt.rate, rec_dose, copula){
  new.subj<-(max(clin_data[,1])+1):(max(clin_data[,1])+n.pat.cohort)
  rate<- sapply(new.subj, function (k) qbinom(copula[1,k], 1, cdlt.rate[rec_dose]))
  mat<-matrix(c(new.subj, rep(rec_dose, times=n.pat.cohort), rep(between.cohort.wk, times=n.pat.cohort),
                rate, rep(max(clin_data[,5])+between.cohort.wk, times=n.pat.cohort)), ncol=5)
  return(rbind(clin_data, mat))
}

############################################

#Description: Simulate PRO data for all patients 

#Input:
#n.cohort: number of patients at each cohort
#n.cohorts.assess: Number of cohorts included in trial
#n.timepoints: Number of timepoints (other than baseline) when PRO data is collected
#pro.schedule: Every number of weeks PRO-nAE data is collected (two weeks in this paper)
#c_dose: vector of recommended doses
#shape_mat: matrix of shape parameters for the doses at each timepoint for PRO-nAE data synthesis
#rate: fixed numerical rate parameter for PRO-nAE data synthesis 
#copula: vector of CDF values to impute correlation between DLT rate and PRO score
#between.cohort.wk: the number of weeks between cohort recruitment (week 4 in this paper)

#Output: matrix; matrix of PRO-nAE score for each patient at each timepoint

pro_sim<- function(n.cohort, n.cohorts.assess, n.timepoints, pro.schedule, c_dose, 
                   shape_mat, rate, copula, between.cohort.wk){
  M<- NULL
  for(i in 1:n.cohorts.assess){
    cohort<- i
    ((cohort-1)*n.cohort+1):(n.cohort*cohort)
    subj<-rep(((cohort-1)*n.cohort+1):(n.cohort*cohort), times=n.timepoints+1)
    currentd<- c_dose[cohort]
    dose<- c(rep(0, times=n.cohort), rep(currentd, times =n.cohort*(n.timepoints)))
    week<- rep(0:n.timepoints, each=n.cohort)
    cop_dim<- rep(2:9, times=3)
    score<-sapply(1:(length(subj)-n.cohort), function (k) qbeta(copula[cop_dim[k],subj[k]], shape1=shape_mat[currentd,week[-(1:n.cohort)][k]+1], shape2=rate))
    score<- c(rbeta(n.cohort, shape_mat[1,1], rate),score)
    week_time <- (week*2)+between.cohort.wk*(cohort-1)
    M<-rbind(M, cbind(subj, dose, rep((0:n.timepoints)*pro.schedule, each=n.cohort), score, week_time)) 
  }
  return(M)
}

################ 

#Description: RJAGS model to use to estimate PRO-nAE burden score 

model<- "model{
#Likelihood for Y
for (i in 1:length(Y)){
Y[i]~dbeta(alpha[i], beta[i])
alpha[i]<-m[i]*phi
beta[i]<- (1-m[i])*phi
logit(m[i])<-  (b0 + u[subj[i],1]) + b1*X[i]+ b2*Z[i] + b3*Z[i]^2
}
#Define random effect on each subject
for (j in n_patients){
    u[j,1] ~ dmnorm(0, sy)
  }
#Prior models
b1~dnorm(0, 1/100)
b2 ~ dnorm(0, 1/100)
b3~dnorm(0, 1/100)
b0 ~ dnorm(0,1/100)
#The variance  sigma^2_gamma is inverse gamma, and so the precision should be gamma
sy~dgamma(0.1,0.1)
phi~dgamma(0.1,0.1)
}"

model_all<- "model{
#Likelihood for Y
for (i in 1:length(Y)){
Y[i]~dbeta(alpha[i], beta[i])
alpha[i]<- m[i]*phi
beta[i]<- (1-m[i])*phi
logit(m[i])<-  b0 + u[subj[i],1] + b1[X[i]+1]+ b2*Z[i] + b3*Z[i]^2 + b4[X[i]+1]*Z[i]
}
#Define random effect on each subject
for (j in n_patients){
    u[j,1] ~ dmnorm(0, sy)
  }
#Prior models
b0 ~ dnorm(0,1/100)
b1[1]<-0
b1[2]~dnorm(0, 1/100)
b1[3]~dnorm(0, 1/100)
b1[4]~dnorm(0, 1/100)
b1[5]~dnorm(0, 1/100)
b1[6]~dnorm(0, 1/100)
b2 ~ dnorm(0, 1/100)
b3~dnorm(0, 1/100)
b4[1]<-0
b4[2]~dnorm(0, 1/100)
b4[3]~dnorm(0, 1/100)
b4[4]~dnorm(0, 1/100)
b4[5]~dnorm(0, 1/100)
b4[6]~dnorm(0, 1/100)
#The variance  sigma^2_gamma is inverse gamma, and so the precision should be gamma
sy~dgamma(0.1, 0.1)
phi~dgamma(0.1,0.1)
}"


#############################

#Description: estimate PRO-nAE score at final PRO assessment using rjags 

#Inputs:
#pro_data: matrix of PRO-nAE data as per `pro_sim` function 
#assessment.time.point: time point at which to estimate PRO-nAE burden score 
#n.iter: Number of iterations of rjags to run 
#runin.prop: Proportion of n.iter to be runin
#n.dose: number of doses

#Output: list; mean estimate for PRO-nAE burden score at each dose and mean phi parameter 

pro_estimate.all<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose){
  jags<- jags.model(textConnection(model_all),
                    data=list(Y=pro_data[,4], X=pro_data[,2], Z=pro_data[,3], n_patients=unique(pro_data[,1]), subj=pro_data[,1]),
                    inits=list(.RNG.name="base::Wichmann-Hill", .RNG.seed=10))
  val<- coda.samples(model=jags, variable.names = c("b0", "b1", "b2", "b3", "b4", "phi"), n.iter=n.iter)
  
  data<-data.frame(val[[1]])
  data<-data[-c(1:(runin.prop*n.iter)),]
  b0<- data$b0
  d1<- data$b1.2.
  d2<- data$b1.3.
  d3<- data$b1.4.
  d4<- data$b1.5.
  d5<- data$b1.6.
  b2<-data$b2
  b3<-data$b3
  intd1<-data$b4.2.
  intd2<-data$b4.3.
  intd3<-data$b4.4.
  intd4<-data$b4.5.
  intd5<-data$b4.6.
  phi<- data$phi
  est<-sapply(1:n.dose, function (k) b0+eval(parse(text=paste0("d",k))) + b2*(assessment.time.point)+b3*(assessment.time.point^2)+
                eval(parse(text=paste0("intd",k)))*assessment.time.point)
  est<- exp(est)/(1+exp(est))
  return(list(colMeans(est), mean(phi)))
}

#############################
#Description: function to synthesis efficacy according to some true value with non-treatment related death

#Input:
#x: value to solve for
#t1: time of first efficacy assessment 
#t2: time of second efficacy assessment 
#prob: true probability of efficacy 

f1 <- function(x,t1, t2, prob) {
  (1-prob)-(t1 + x*(1-t1))*(t2+x*(1-t2))
}

#############################################

#Description: Simulate efficacy observations for patients 

#Input: 
#n.cohort: number of patients at each cohort
#n.cohorts.assess: number of cohorts included in trial
#prob_efficacy: vector of true probabilities of efficacy 
#c_dose: vector of recommended doses
#between.cohort.wk: the number of weeks between cohort recruitment (week 4 in this paper)
#eff.schedule: vector of time (weeks) of efficacy assessments 
#med_survival_month: median time (weeks) of survival

#Output: matrix;  matrix of efficacy responses for each patient at each efficacy assessment 

eff_sim<- function(n.cohort, n.cohorts.assess, prob_efficacy, c_dose,between.cohort.wk, eff.schedule, med_survival_month){
  prob_efficacy_sim<-sapply(1:5, function (k) 1-uniroot(f1, t1=pexp(eff.schedule[1], rate=(log(2)/(med_survival_month*4))), 
                                                        t2= pexp(eff.schedule[2], rate=(log(2)/(med_survival_month*4))), prob=prob_efficacy[k], c(0,1))$root)
  
  M<- matrix(nrow=n.cohort*n.cohorts.assess, ncol=6)
  colnames(M)<-c("subj", "dose", "eff1", "eff2", "best.eff", "week_timeline")
  M[,1]<- 1:(n.cohort*n.cohorts.assess)
  M[,2]<- rep(c_dose, each=n.cohort)
  M[,3]<- sapply(1:n.cohorts.assess, function (k) rbinom(n.cohort, 1, prob_efficacy_sim[c_dose[k]]))
  M[,4]<- sapply(1:n.cohorts.assess, function (k) rbinom(n.cohort, 1, prob_efficacy_sim[c_dose[k]]))
  M[,5]<- pmax(M[,3],M[,4])
  M[,6]<- eff.schedule[2]+between.cohort.wk*(0:(n.cohorts.assess-1))
  return(M)
}  

###########################

#Description: Estimate probability of efficacy using iPipe method

#Input:
#eff_data: matrix of efficacy responses for each patient at each efficacy assessment  as per `pro_sim` function 
#e: value of epsilon
#a: shape parameter for beta prior
#b: rate parameter for beta prior
#n_t_grid: number of time points to evaluate over 

#Output: vector; estimated probability of efficacy for each dose as per iPipe

pipe_est<- function(eff_data,e, a, b, n_t_grid){
  ext<-sapply(1:5, function (k) na.omit(eff_data[eff_data[,2]==k,5]))
  sum <- sapply(ext, sum)
  n <- sapply(ext, length)
  Gamma<- t(matrix(c(0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1), nrow=5))
  t<- seq(from=0, to =1, length.out=n_t_grid)
  vec<-sapply(1:n_t_grid, function (i) which.max(sapply(1:6, function (j) prod(sapply(1:5, function (k) ((e*(1-pbeta(t[i], a+sum[k], n[k]-sum[k]+b)))^Gamma[j,k])*(
    ((1-e)*(pbeta(t[i], a+sum[k], n[k]-sum[k]+b)))^(1-Gamma[j,k])))))))
  change_points <- c(TRUE, diff(vec) != 0)
  new<-Gamma[vec[change_points],]
  first_row <- sapply(1:5, function (k) which(new[,k] != 1)[1])
  prob<-(t[change_points])[first_row]
  return(prob)
}

##########################

#Description: Samples from the posterior estimates of probability of efficacy as per iPipe design

#Input:
#eff_outcome: vector of patients responses for all patients at specified dose
#alpha: shape parameter for beta prior
#beta: rate parameter for beta prior
#n.sim.final: number of samples to collect from the posterior distribution


eff_estimate_sim<- function(eff_outcome, alpha, beta, n.sim.final){
  ep<- runif(n.sim.final)
  val<-sapply(1:n.sim.final, function(i) pipe_est(eff_outcome, ep[i], 0.5, 0.5, 100))
  return(t(val))
}

###################

#Description: calculate the expected loss for each dose using both efficacy and PRO-nAE burden score 

#Inputs:
#pro_data: matrix of PRO-nAE data as per `pro_sim` function 
#assessment.time.point: time point at which to estimate PRO-nAE burden score 
#n.iter: Number of iterations of rjags to run 
#runin.prop: Proportion of n.iter to be runin
#n.dose: number of doses
#eff_data: vector of patients responses for all patients at specified dose
#a: shape parameter for beta prior in iPipe
#b: rate parameter for beta prior in iPipe
#n.sim.final: number of samples to collect from the posterior distribution

#Output: list; estimated loss, PRO-nAE score and probability of response for each dose

loss.all<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, n.sim.final){
  pro<- pro_estimate.all(pro_data, assessment.time.point, n.iter, runin.prop, n.dose)
  pro_sample<- sapply(1:n.dose, function(k) rbeta(n.sim.final,pro[[1]][k]*pro[[2]],(1-pro[[1]][k])*pro[[2]]))
  eff<-eff_estimate_sim(eff_data, alpha, beta, n.sim.final)
  eff.est<-pipe_est(eff_data,0.5,alpha, beta, 100)
  loss.val<- matrix(nrow=n.sim.final, ncol=n.dose)
  for(i in 1:n.sim.final){
    loss.val[i,]<- ((pro_sample[i,]^2) + (eff[i,]-1)^2)^(1/2)
  }
  return(list(colMeans(loss.val), pro[[1]], eff.est))
}

#########################################

#Description: Recommend dose for stage 2 of design 

#Inputs:
#n.dose: number of doses
#target: target DLT rate
#dose: vector of dose recommendations for each patient
#cdlt: vector of cdlt observations for each patient 
#esc_bound: escalation threshold as per BOIN
#alpha: shape parameter for prior as per safety stopping rule 
#beta: rate parameter for prior as per safety stopping rule 
#eff_admiss: vector of admissible doses as per futility stopping rule
#n.admiss: number of admissible doses 

#Output: list; list of recommended dose and vector of admissible doses 

recommendation<- function(n.dose, target, dose, cdlt, esc_bound, alpha, beta, eff_admiss, n.admiss){
  admiss<- intersect(boin_decision(target, dose, cdlt, esc_bound, n.dose, alpha, beta), eff_admiss)
  #ifelse(length(admiss)==1, rec<- admiss, rec<- sample(admiss, size=1, prob=sqrt(2)-loss.est[[1]][admiss], replace=TRUE))
  if(length(admiss)==0){
    return(list(0,admiss))
  } else {
    ifelse(length(admiss)==1, rec<- admiss, rec<- sample(admiss, size=1, prob=(1/n.admiss)[admiss]))
    return(list(rec,admiss))
  }
}

################

#Description: calculate the final recommended dose using the expected loss for each dose using both efficacy and PRO-nAE burden score 

#Inputs:
#pro_data: matrix of PRO-nAE data as per `pro_sim` function 
#assessment.time.point: time point at which to estimate PRO-nAE burden score 
#n.iter: Number of iterations of rjags to run 
#runin.prop: Proportion of n.iter to be runin
#n.dose: number of doses
#eff_data: vector of patients responses for all patients at specified dose
#alpha: shape parameter for beta prior in iPipe
#beta: rate parameter for beta prior in iPipe
#target: target DLT rate
#dose: vector of dose recommendations for each patient
#cdlt: vector of cdlt observations for each patient 
#n.sim.final: number of samples to collect from the posterior distribution

#Output: list; recommended dose, estimated loss, PRO-nAE score and probability of response for each dose


final.recommendation.all<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, target, dose, cdlt,n.sim.final){
  loss.est<-loss.all(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, n.sim.final)
  admiss<- boin_admiss(target, dose, cdlt, alpha, beta)
  rec<- which(loss.est[[1]]==min(loss.est[[1]][admiss]))
  return(list(rec, loss.est[[1]], loss.est[[2]], loss.est[[3]]))
}

##################

#Description: find admissible doses as per safety stopping rule

#Inputs:
#target: target DLT rate
#dose: vector of dose recommendations for each patient
#cdlt: vector of cdlt observations for each patient 
#alpha: shape parameter for prior as per stopping rule 
#beta: rate parameter for prior as per stopping rule 

#Output:vector; admissible doses as per the safety stopping rule 

boin_admiss<- function(target, dose, cdlt, alpha, beta){
  max.d<-max(dose)
  admiss<- c()
  admiss<-which(sapply(1:max.d, function (k) 1-pbeta(target,1+sum(cdlt[which(dose==k)]),length(which(dose==k))-sum(cdlt[which(dose==k)])+ beta)<=0.95)==TRUE)
  return(admiss)
}

####################

#Description: Complete trial simulation 

#Inputs:
#general_ls: list containing all general parameters
    # n.patient.cohort: number of patients at each cohort
    # first.week.assessed: first week DLTs are assessed (week 4 in this paper)
    # between.cohort.wk: the number of weeks between cohort recruitment (week 4 in this paper)
    # n.doses: number of doses
    # n.timepoints: Number of timepoints (other than baseline) when PRO data is collected
    # final.assessment.timepoint: final time point at which to estimate PRO-nAE burden score 
    # mcmc.niter: Number of iterations of rjags to run 
    # mcmc.burnin.prop: Proportion of n.iter to be runin
    # n.cohorts.all:Number of cohorts included in trial
    # n.sim.final: number of samples to collect from the posterior distribution
    # cop.corr: correlation parameter for copula 


#boin_ls: list containing all dose-escalation parameters (Stage 1)
    # cdlt_rates: vector of true DLT rates for each dose
    # esc_bound:escalation threshold as per BOIN
    # des_bound: de-escalation threshold as per BOIN
    # target: target DLT rate
    # beta_a_safety: shape parameter for prior as per safety stopping rule 
    # beta_b_safety: rate parameter for prior as per safety stopping rule 

#pro_ls: list containing all parameters associated with PRO-nAE burden score
    # pro.schedule: Every number of weeks PRO-nAE data is collected (two weeks in this paper)
    # beta_shape_sc: marix of shape parameters for the doses at each timepoint for PRO-nAE data synthesis
    # beta_rate_sc: fixed numerical rate parameter for PRO-nAE data synthesis 
    # pro.between.cohort.timepoints: the number of weeks between cohort recruitment (week 4 in this paper)

#eff_ls: list containing all parameters associated with efficacy 
# eff.schedule:  vector of time (weeks) of efficacy assessments 
# eff_rates: vector of true probabilities of efficacy 
# beta_a: shape parameter for prior as per futility stopping rule 
# beta_b: rate parameter for prior as per futility stopping rule 
# min_eff: minimum acceptable probability of efficacy 
# interim_complete_cohort1:
# interim_complete_cohort2:

#data_after_dlt: TRUE/FALSE indicating if PRO/efficacy data is collected/available beyond a DLT observation

#Output: list; containing:
  #final.rec: vector of final recommendation 
  #boin.admiss: vector of final set of admissible doses
  #dose.explored: vector of doses explored within trial
  #pat.allocated: vector of number of patients allocated to each dose
  #cdlt: matrix of patients with DLT observations
  #eff: matrix of patients with efficacy observations 
  #pro: matrix of patients with pro observations
  #loss_est: vector of expected loss for each dose
  #eff_est: vector of expected efficacy for each dose
  #pro_est: vector of expected PRO-nAE score for each dose
  #n_cens: vector of number of patients experiencing non-treatment related death for each dose
  #n_dlt: vector of number of dlts at each dose  

 
trial_design<- function(general_ls, boin_ls, pro_ls, eff_ls, data_after_dlt){
  #general inputs
  n.patient.cohort<- general_ls[[1]]
  first.week.assessed<- general_ls[[2]]
  between.cohort.wk<-general_ls[[3]]
  n.doses<-general_ls[[4]]
  n.timepoints<-general_ls[[5]]
  final.assessment.timepoint<-general_ls[[6]]
  mcmc.niter<-general_ls[[7]]
  mcmc.burnin.prop<- general_ls[[8]]
  n.cohorts.all<- general_ls[[9]]
  cop.corr<- general_ls[[10]]
  n.sim.final<- mcmc.niter*(1-mcmc.burnin.prop)
  
  #dlt inputs 
  cdlt_rates<- boin_ls[[1]]
  esc_bound<- boin_ls[[2]]
  des_bound<- boin_ls[[3]]
  target<- boin_ls[[4]]
  beta_a_safety<- boin_ls[[5]]
  beta_b_safety<- boin_ls[[6]]
  
  #pro_inputs
  pro.schedule<-pro_ls[[1]]
  beta_shape_sc<- pro_ls[[2]]
  beta_rate_sc<- pro_ls[[3]]
  pro.between.cohort.timepoints<-pro_ls[[4]]
  
  #eff inputs
  eff.schedule<- eff_ls[[1]]
  eff_rates<-eff_ls[[2]]
  beta_a<- eff_ls[[3]]
  beta_b<- eff_ls[[4]]
  min_eff<- eff_ls[[5]]
  interim_complete_cohort1<- eff_ls[[6]]
  interim_complete_cohort2<- eff_ls[[7]]
  
  ###stage 1 - create data up until cohort 6
  mat<-diag(2, nrow=n.timepoints+1, ncol=n.timepoints+1)
  mat<-matrix(rep(0.9, times=(n.timepoints+1)^2), nrow=n.timepoints+1)
  diag(mat)<- 1
  t<-rmvnorm(n=n.cohorts.all*n.patient.cohort,mean = rep(0, times=n.timepoints+1), sigma = mat)
  copula<-sapply(1:(n.cohorts.all*n.patient.cohort), function (k) pnorm(t[k,], mean = rep(0), sd = sqrt(1)))
  
  
  #compute time of random non-treatment related death
  median_survival_month<- (general_ls[[11]]-1)
  lambda<- log(2)/(median_survival_month*4)
  pat_cens<-rexp(n.patient.cohort*n.cohorts.all, lambda)+first.week.assessed
  wk_interim1<- eff.schedule[2]+(interim_complete_cohort1-1)*between.cohort.wk
  wk_interim2<- eff.schedule[2]+(interim_complete_cohort2-1)*between.cohort.wk
  cohort_allot_interim1<-(wk_interim1/between.cohort.wk)+1 
  cohort_allot_interim2<-(wk_interim2/between.cohort.wk)+1 
  clin_set<-boin(cdlt_rates, n.patient.cohort, first.week.assessed, 
                 between.cohort.wk, n.doses, esc_bound, des_bound, target, copula, beta_a_safety, beta_b_safety, cohort_allot_interim1)
  clin_mat<- clin_set[[1]]
  initial_rec<- clin_set[[2]]
  week<- ((nrow(clin_mat)/n.patient.cohort+1)*between.cohort.wk)-between.cohort.wk
  if(initial_rec[length(initial_rec)]==0){
    return(list(final.rec=NA, boin.admiss=NA, dose.explored=NA, pat.allocated=NA,
                cdlt = NA, eff=NA, pro=NA))
  }
  #make next dose decision using BOIN data  
  eff_admiss<- 1:n.doses
  next.recommendation<-recommendation(n.doses,target, clin_mat[,2], clin_mat[,4], esc_bound, beta_a_safety, beta_b_safety, eff_admiss, as.numeric(table(clin_mat[,2])))
  next.dose<-next.recommendation[[1]] 
  boin.admiss<-next.recommendation[[2]]
  while(nrow(clin_mat)<n.cohorts.all*n.patient.cohort){
    #update tables for next dose recommendation
    clin_mat<-boin_next(clin_mat, n.patient.cohort, between.cohort.wk, cdlt_rates,next.dose, copula)
    next.recommendation<-recommendation(n.doses,target, clin_mat[,2], clin_mat[,4], esc_bound, beta_a_safety, beta_b_safety, eff_admiss, as.numeric(table(clin_mat[,2])))
    next.dose<-next.recommendation[[1]] 
    if(next.dose[length(next.dose)]==0){
      return(NA)
    }
    if(max(clin_mat[,1])/n.patient.cohort==cohort_allot_interim1){
      eff<- eff_sim(n.patient.cohort, interim_complete_cohort1, eff_rates, clin_mat[(1:interim_complete_cohort1)*n.patient.cohort,2], between.cohort.wk, eff.schedule, general_ls[[11]]-1)
      week_trunc<-pat_cens[1:(n.patient.cohort*interim_complete_cohort1)]
      for(k in 1:max(eff[,1])){
        if(week_trunc[k]<eff.schedule[2]){
          if(week_trunc[k]>=eff.schedule[1]+1 & week_trunc[k]<eff.schedule[2]){
            eff[k, 5]<-eff[k,3] 
            eff[k, 4]<-0
          }else{
            eff[k, 5]<-0
            eff[k, 4]<-0
            eff[k, 3]<-0
          }
        }
      }
      if(data_after_dlt==FALSE && sum(clin_mat[,4])>0){
        eff_dlt<-which(clin_mat[1:(n.patient.cohort*interim_complete_cohort1),4]==1)
        pat_rm<- which(eff[,1] %in% eff_dlt)
        if(length(pat_rm)>0){
          eff[pat_rm,c(3,4,5)]<- rep(0, times=3)
        }
      }
      eff_admiss<-futility_decision(min_eff, eff[,2], eff[,5], beta_a, beta_b, n.doses)
    }
    if(max(clin_mat[,1])/n.patient.cohort==cohort_allot_interim2){
      eff<- rbind(eff, eff_sim(n.patient.cohort, cohort_allot_interim2-cohort_allot_interim1, eff_rates, clin_mat[((interim_complete_cohort1+1):(interim_complete_cohort2))*n.patient.cohort,2],
                               between.cohort.wk, eff.schedule, general_ls[[11]]-1))
      eff[,1]<- 1:(interim_complete_cohort2*n.patient.cohort)
      week_trunc<-pat_cens[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2)]
      for(k in (n.patient.cohort*interim_complete_cohort1+1):max(eff[,1])){
        i<- k-n.patient.cohort*interim_complete_cohort1
        if(week_trunc[i]<eff.schedule[2]){
          if(week_trunc[i]>=eff.schedule[1]+1 & week_trunc[i]<eff.schedule[2]){
            eff[k, 5]<-eff[k,3] 
            eff[k, 4]<-0
          }else{
            eff[k, 5]<-0
            eff[k, 4]<-0
            eff[k, 3]<-0
          }
        }
      }
      if(data_after_dlt==FALSE && sum(clin_mat[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2),4])>0){
        eff_dlt<-which(clin_mat[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2),4]==1)
        pat_rm<- which(eff[,1] %in% (eff_dlt+(n.patient.cohort*interim_complete_cohort1)))
        if(length(pat_rm)>0){
          eff[pat_rm,c(3,4,5)]<- rep(0, times=3)
        }
      }
      eff_admiss<-futility_decision(min_eff, eff[,2], eff[,5], beta_a, beta_b, n.doses)
    }
    boin.admiss<-intersect(next.recommendation[[2]], eff_admiss)
    next.recommendation<-recommendation(n.doses,target, clin_mat[,2], clin_mat[,4], esc_bound, beta_a_safety, beta_b_safety, eff_admiss, as.numeric(table(clin_mat[,2])))
    next.dose<-next.recommendation[[1]] 
    if(next.dose[length(next.dose)]==0){
      return(NA)
    }
    print(next.dose)
  }
  #stage 2 - evaluate right at the end of the trial 
  c_dose<- clin_mat[(1:n.cohorts.all)*n.patient.cohort,2]
  #create pro matrix
  pro<-pro_sim(n.patient.cohort, n.cohorts.all, n.timepoints, pro.schedule, c_dose, 
               beta_shape_sc, beta_rate_sc, copula, between.cohort.wk)
  
  eff<- rbind(eff, eff_sim(n.patient.cohort, length((interim_complete_cohort2+1):n.cohorts.all), eff_rates, c_dose[(interim_complete_cohort2+1):n.cohorts.all],
                           between.cohort.wk, eff.schedule, general_ls[[11]]-1))
  eff[,1]<- 1:(n.cohorts.all*n.patient.cohort)
  week_trunc<-pat_cens[(n.patient.cohort*interim_complete_cohort2+1):(n.patient.cohort*n.cohorts.all)]
  for(k in (n.patient.cohort*interim_complete_cohort2+1):max(eff[,1])){
    i<- k-n.patient.cohort*interim_complete_cohort2
    if(week_trunc[i]<eff.schedule[2]){
      if(week_trunc[i]>=eff.schedule[1]+1 & week_trunc[i]<eff.schedule[2]){
        eff[k, 5]<-eff[k,3] 
        eff[k, 4]<-0
      }else{
        eff[k, 5]<-0
        eff[k, 4]<-0
        eff[k, 3]<-0
      }
    }
  }
  row.trunc<-c(sapply(1:max(pro[,1]), function (k) pro[pro[,1]==k,3]<=pat_cens[k]))
  pro<- pro[order(pro[,1]),]
  pro<-pro[row.trunc,]
  if(data_after_dlt==FALSE && sum(clin_mat[,4])>0){
    dlt<-which(clin_mat[,4]==1)
    pat_rm<- intersect(which(pro[,1] %in% dlt), which(pro[,3]>4))
    if(length(pat_rm)>0){
      pro<-pro[-pat_rm,]
    }
    eff_dlt<-which(clin_mat[(n.patient.cohort*interim_complete_cohort2+1):(n.patient.cohort*n.cohorts.all),4]==1)
    pat_rm<- which(eff[,1] %in% (eff_dlt+n.patient.cohort*interim_complete_cohort2))
    if(length(pat_rm)>0){
      eff[pat_rm,c(3,4,5)]<- rep(0, times=3)
    }
  }
  #find final recommendation
  final.run<- final.recommendation.all(pro, final.assessment.timepoint, mcmc.niter, mcmc.burnin.prop, n.doses, eff,
                                       beta_a, beta_b, target, clin_mat[,2], clin_mat[,4], n.sim.final)
  final<- final.run[[1]]
  final.admiss<-intersect(boin_admiss(target,clin_mat[,2],clin_mat[,4], beta_a, beta_b), futility_decision(min_eff,eff[,2],eff[,5], beta_a, beta_b, n.doses))
  explored<- unique(clin_mat[,2])
  allocated<- sapply(1:n.doses, function (k) nrow(clin_mat[clin_mat[,2]==k,]))
  cens_val<-cbind(as.numeric(pat_cens<=eff.schedule[2]),clin_mat[,2] )
  ncens<-sapply(1:n.doses, function (k) sum(cens_val[cens_val[,2]==k,1]))[explored]
  dlt_val<- cbind(clin_mat[,4],clin_mat[,2])
  ndlt<-sapply(1:n.doses, function (k) sum(dlt_val[dlt_val[,2]==k,1]))[explored]

  # return(list(final.rec=final, boin.admiss=final.admiss, dose.explored=explored, pat.allocated=allocated, pro.est=pro.estimate, 
  #             eff.est=efficacy.estimate, loss=loss.val, cdlt = clin_mat, eff=cdlt_eff_mat, pro=cdlt_pro_mat))
  return(list(final.rec=final, 
    boin.admiss=final.admiss, dose.explored=explored, pat.allocated=allocated[explored],
    cdlt = clin_mat,eff=eff, pro=pro, loss_est=final.run[[2]][explored], eff_est=final.run[[4]][explored], pro_est=final.run[[3]][explored], 
    n_cens=ncens,n_dlt= ndlt))
}
