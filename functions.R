library(rjags)
library(tidyverse)
library(clusterGeneration)
library(parallel)
library(fitdistrplus)
library(VineCopula)
library(mvtnorm)
######################inverse gamma posterior################################
posterior<- function(mu, sigma_2, mu_0j, lambda_j, zeta_j, zeta_dash_j){
  pos<- (sqrt(lambda_j)/sqrt(2*sigma_2*pi))*((zeta_dash_j^zeta_j)/gamma(zeta_j))*(1/sigma_2)^(zeta_j+1)*
    exp(-((2*zeta_dash_j)+lambda_j*(mu-mu_0j)^2)/(2*sigma_2))
  return(pos)
}

posterior_mean<- function(mu_0j, lambda_j, x_bar_j, n_j){
  pos_m<- (lambda_j*mu_0j + (n_j*x_bar_j))/(lambda_j + n_j)
  return(pos_m)
}

posterior_zeta_j<- function(n_j, zeta_j){
  pos_z<- (n_j/2)+zeta_j
  return(pos_z)
}

posterior_lambda_j<-function(lambda_j, n_j){
  pos_l<- lambda_j+n_j
  return(pos_l)
}

posterior_zeta_j_dash <- function(mu_0j, lambda_j, x_bar_j, n_j, zeta_j_dash, x_j){
  pos_z_dash<- zeta_j_dash + (1/2)*sum((x_j-x_bar_j)^2) + ((n_j*lambda_j)/(lambda_j + n_j))*((x_bar_j-mu_0j)^2)/2 
  return(pos_z_dash)
}


###################################################################
############################## SIMULATING PRO-CTCAE ##########################################


#for a vector of probabilities of toxicity for dose 1, induce a dose effect 
dose_effect<- function(a,b, vec_initial_probs, modified_a){
  ngrades<- length(vec_initial_probs)
  C<- sapply(1:(ngrades-1), function(j) qbeta(sum(vec_initial_probs[1:j]), a, b))
  vec<- c()
  vec[1]<-pbeta(C[1], modified_a, b)
  for(i in 2:(ngrades-1)){
    vec[i]<- pbeta(C[i], modified_a, b)-pbeta(C[i-1],modified_a, b)
  }
  vec[ngrades]<- 1-pbeta(C[ngrades-1], modified_a,b) 
  return(vec)
}

generate_matrix_baseline<- function(ngrades, a, b){
  u<-sort(runif(ngrades-1))
  baseline<- c(pbeta(u[1],a,b), pbeta(u[2],a,b)-pbeta(u[1],a,b), pbeta(u[3],a,b)-pbeta(u[2],a,b),pbeta(u[4],a,b)-pbeta(u[3],a,b),pbeta(u[4],a,b, lower.tail = FALSE))
  return(baseline)
}

generate_matrix_beta<- function(ndose, ngrades, vec_a, a, b){
  u<-sort(runif(ngrades-1))
  dose1<- c(pbeta(u[1],a,b), pbeta(u[2],a,b)-pbeta(u[1],a,b), pbeta(u[3],a,b)-pbeta(u[2],a,b),pbeta(u[4],a,b)-pbeta(u[3],a,b),pbeta(u[4],a,b, lower.tail = FALSE))
  M<- rbind(dose1, t(sapply(1:(ngrades-1), function(k) dose_effect(a,b, dose1, vec_a[k]))))
  return(M)
}

mod_a<- seq(from=1, to =1.5, length.out=5)[-1]
n.toxicity<- 78
ls<- vector(mode='list', length=n.toxicity)
for(i in 1:n.toxicity){
  ls[[i]]<- generate_matrix_beta(5,5,mod_a, 1, 4)
}

#for a matrix of probabilities of toxicity for multiple doses, induce a time effect 
time_effect<- function(matrix_tox, b, vec_modified_a, a, modified_b){
  vec<-c()
  ngrades=ncol(matrix_tox)
  dose<-0
  M<- matrix(nrow=nrow(matrix_tox), ncol=ncol(matrix_tox))
  while(dose< nrow(matrix_tox)){
    dose<- dose+1
    modified_a<-c(a, vec_modified_a)[dose] 
    C<- sapply(1:(ngrades-1), function(j) qbeta(sum(matrix_tox[dose,1:j]), modified_a, b))
    vec[1]<-pbeta(C[1], modified_a, modified_b)
    for(i in 2:(ngrades-1)){
      vec[i]<- pbeta(C[i], modified_a, modified_b)-pbeta(C[i-1], modified_a, modified_b)
    }
    vec[ngrades]<- 1-pbeta(C[ngrades-1], modified_a, modified_b) 
    M[dose,]<- vec
  }
  return(M)
}


################# Simulate AE burden score totally ##########################

#generate 78 matrices showing probability of each severity grade for n.doses
# n.toxicity<- 78
# ls<- vector(mode='list', length=n.toxicity)
# for(i in 1:n.toxicity){
#   ls[[i]]<- generate_matrix(5, 5, c(0.1,0.2,0.3,0.4))
# }

# For each PRO-CTCAE item, create a new matrix showing the toxicity probability across each cycle. 
ae.score.cycle<- function(toxicity_matrix, a, vec_a,b, vec_modified_cycle_b){
  cyc<-lapply(c(b,vec_modified_cycle_b), function (k) time_effect(toxicity_matrix,b, vec_a,a,k))
  return(cyc)
}

ae_score_baseline_wish<- function(n.dose, n.toxicity, toxicity_mat, sigma){
  z<- mvrnorm(mu=rep(0, times=n.toxicity), Sigma=sigma)
  u<-sapply(1:n.toxicity, function (k) pnorm(z[k], mean =0, sd=sqrt(diag(sigma)[k])))
  scores<- sapply(1:n.toxicity, function (j) min(which(sapply(1:n.dose, function (k) sum(toxicity_mat[[j]][1:k]))>u[j])))
  return(scores)
}

ae_score_patient_wish<- function(cycle, dose, n.toxicity, toxicity_dose_grade_cycle_mat, sigma){
  z<- mvrnorm(mu=rep(0, times=n.toxicity), Sigma=sigma)
  u<-sapply(1:n.toxicity, function (k) pnorm(z[k], mean =0, sd=sqrt(diag(sigma)[k])))
  scores<- sapply(1:n.toxicity, function (j) min(which(sapply(1:5, function(k) sum(toxicity_dose_grade_cycle_mat[,j][[cycle]][dose,1:k]))>u[j])))
  return(scores)
}

#generate matrix of shapes and rates to generate the PRO data 
beta_matrix<- function(n.toxicity, shape_dose_vec, rate_cycle_vec, baseline_a, baseline_b, n.sample, n.grade, 
                       n.dose, n.timepoints){
  M<- matrix(nrow=n.toxicity, ncol=5)
  ls<- vector(mode='list', length=n.toxicity)
  for(i in 1:n.toxicity){
    ls[[i]]<- generate_matrix_beta(n.dose, n.grade, shape_dose_vec[-1],shape_dose_vec[1],rate_cycle_vec[1])
  }
  
  ls_baseline<- vector(mode='list', length=n.toxicity)
  for(i in 1:n.toxicity){
    ls_baseline[[i]]<-generate_matrix_baseline(n.grade, baseline_a,baseline_b)
  }
  cov<-lapply(rep(n.toxicity, times=n.sample), function(k) genPositiveDefMat(dim = k,covMethod = "eigen")$Sigma)
  prob_t<-sapply(1:n.toxicity, function (k) ae.score.cycle(ls[[k]], shape_dose_vec[1], shape_dose_vec[-1], rate_cycle_vec[1], rate_cycle_vec[-1]))
  baseline_reading<-sapply(1:n.sample, function (k) ae_score_baseline_wish(n.grade, n.toxicity, ls_baseline, cov[[k]]))
  baseline_score<- colSums(baseline_reading-1)/(n.toxicity*(n.grade-1))
  if(min(baseline_score)==0){
    baseline_score[which(baseline_score==0)]<- 0.001
  }
  baseline_beta<-fitdist(baseline_score, "beta")
  fix.sh2<-as.numeric(baseline_beta$estimate[2])
  b<- matrix(nrow=n.dose, ncol=n.timepoints+1)
  b[,1]<- rep(baseline_beta$estimate[1], times=nrow(b))
  for(i in 1:n.dose){
    for(j in 1:n.timepoints){
      currentd<-i
      cycle<-j
      M_score<-sapply(1:n.sample, function (k) ae_score_patient_wish(cycle, currentd, n.toxicity,prob_t, cov[[k]]))
      score<- colSums(M_score-1)/(n.toxicity*(n.grade-1))
      if(min(score)==0){
        score[which(score==0)]<- 0.001
      }
      b[i,j+1]<-fitdist(score, "beta", fix.arg = list(shape2=fix.sh2))$estimate[1]
    }
  }
  return(list(shape_matrix=b, fixed_rate=fix.sh2))
}

######## BOIN functions 
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

futility_decision<- function(min_target, dose, eff, alpha, beta, n.dose){
  not.admiss<-sapply(1:n.dose, function (k) (pbeta(min_target,alpha+sum(eff[which(dose==k)]),length(which(dose==k))-sum(eff[which(dose==k)])+ beta)>0.7)&(sum(dose==k)>=6))
  admiss<- (1:n.dose)[!not.admiss]
  return(admiss)
}

boin_next<- function(clin_data, n.pat.cohort, between.cohort.wk, cdlt.rate, rec_dose, copula){
  new.subj<-(max(clin_data[,1])+1):(max(clin_data[,1])+n.pat.cohort)
  rate<- sapply(new.subj, function (k) qbinom(copula[1,k], 1, cdlt.rate[rec_dose]))
  mat<-matrix(c(new.subj, rep(rec_dose, times=n.pat.cohort), rep(between.cohort.wk, times=n.pat.cohort),
                rate, rep(max(clin_data[,5])+between.cohort.wk, times=n.pat.cohort)), ncol=5)
  return(rbind(clin_data, mat))
}

## PRO data generation functions 
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


pro_estimate<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose){
  jags<- jags.model(textConnection(model),
                    data=list(Y=pro_data[,4], X=pro_data[,2], Z=pro_data[,3], n_patients=unique(pro_data[,1]), subj=pro_data[,1]),
                    inits=list(.RNG.name="base::Wichmann-Hill", .RNG.seed=10))
  val<- coda.samples(model=jags, variable.names = c("b0", "b1", "b2", "b3", "phi"), n.iter=n.iter)
  data<-data.frame(val[[1]])
  data<-data[-c(1:(runin.prop*n.iter)),]
  b0<- data$b0
  b1<- data$b1
  b2<- data$b2
  b3<- data$b3
  phi<- data$phi
  est<-sapply(1:n.dose, function (k) b0+(b1*k)+(b2*assessment.time.point)+b3*(assessment.time.point^2))
  est<- exp(est)/(1+exp(est))
  return(list(colMeans(est), mean(phi)))
}

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

##efficacy functions 
f1 <- function(x,t1, t2, prob) {
  (1-prob)-(t1 + x*(1-t1))*(t2+x*(1-t2))
}

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


eff_estimate_sim<- function(eff_outcome, alpha, beta, n.sim.final){
  ep<- runif(n.sim.final)
  val<-sapply(1:n.sim.final, function(i) pipe_est(eff_outcome, ep[i], 0.5, 0.5, 100))
  return(t(val))
}

beta_binom_eff<- function(eff_data, alpha, beta, n.sim.final){
  ext<-sapply(1:5, function (k) na.omit(eff_data[eff_data[,2]==k,5]))
  sum <- sapply(ext, sum)
  n <- sapply(ext, length)
  alpha_new<-sapply(1:5, function (k) sum[k]+alpha)
  beta_new<- sapply(1:5, function (k) n[k] - sum[k]+ beta)
  beta.sim<-sapply(1:5, function (k) rbeta(n.sim.final, alpha_new[k], beta_new[k]))
  mean.eff<-sapply(1:5, function (k) alpha_new[k]/(alpha_new[k]+beta_new[k]))
  return(list(beta.sim,mean.eff))
}



################ RJAGS 

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

loss.all<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, n.sim.final){
  pro<- pro_estimate.all(pro_data, assessment.time.point, n.iter, runin.prop, n.dose)
  pro_sample<- sapply(1:n.dose, function(k) rbeta(n.sim.final,pro[[1]][k]*pro[[2]],(1-pro[[1]][k])*pro[[2]]))
  eff<-eff_estimate_sim(eff_data, alpha, beta, n.sim.final)
  eff.est<-colMeans(eff)
  loss.val<- matrix(nrow=n.sim.final, ncol=n.dose)
  for(i in 1:n.sim.final){
    loss.val[i,]<- ((pro_sample[i,]^2) + (eff[i,]-1)^2)^(1/2)
  }
  return(list(colMeans(loss.val), pro[[1]], eff.est))
}

loss.all.bb<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, n.sim.final){
  pro<- pro_estimate.all(pro_data, assessment.time.point, n.iter, runin.prop, n.dose)
  pro_sample<- sapply(1:n.dose, function(k) rbeta(n.sim.final,pro[[1]][k]*pro[[2]],(1-pro[[1]][k])*pro[[2]]))
  efficacy<- beta_binom_eff(eff_data, alpha, beta, n.sim.final)
  eff<-efficacy[[1]]
  eff.est<-efficacy[[2]]
  loss.val<- matrix(nrow=n.sim.final, ncol=n.dose)
  for(i in 1:n.sim.final){
    loss.val[i,]<- ((pro_sample[i,]^2) + (eff[i,]-1)^2)^(1/2)
  }
  return(list(colMeans(loss.val), pro[[1]], eff.est))
}

# final.recommendation<- function(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta, target, dose, cdlt,n.sim.final){
#   loss.est<-loss(pro_data, assessment.time.point, n.iter, runin.prop, n.dose, eff_data, alpha, beta,n.sim.final)
#   admiss<- boin_admiss(target, dose, cdlt)
#   rec<- which(loss.est==min(loss.est[admiss]))
#   return(rec)
# }

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

#final admiss set
boin_admiss<- function(target, dose, cdlt, alpha, beta){
  max.d<-max(dose)
  admiss<- c()
  admiss<-which(sapply(1:max.d, function (k) 1-pbeta(target,1+sum(cdlt[which(dose==k)]),length(which(dose==k))-sum(cdlt[which(dose==k)])+ beta)<=0.95)==TRUE)
  return(admiss)
}


#TRIAL DESIGN 
trial_design_treatment_policy<- function(general_ls, boin_ls, pro_ls, eff_ls, data_after_dlt){
  #general inputs
  n.patient.cohort<- general_ls[[1]]
  first.week.assessed<- general_ls[[2]]
  between.cohort.wk<-general_ls[[3]]
  n.doses<-general_ls[[4]]
  n.timepoints<-general_ls[[5]]
  n.toxicity<- general_ls[[6]]
  n.grades<- general_ls[[7]]
  final.assessment.timepoint<-general_ls[[8]]
  mcmc.niter<-general_ls[[9]]
  mcmc.burnin.prop<- general_ls[[10]]
  n.cohorts.all<- general_ls[[11]]
  n.sim.final<- general_ls[[12]]
  cop.corr<- general_ls[[13]]
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
  max.pro<- pro_ls[[5]]
  
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
  median_survival_month<- (general_ls[[14]]-1)
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
      eff<- eff_sim(n.patient.cohort, interim_complete_cohort1, eff_rates, clin_mat[(1:interim_complete_cohort1)*n.patient.cohort,2], between.cohort.wk, eff.schedule, general_ls[[14]]-1)
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
                               between.cohort.wk, eff.schedule, general_ls[[14]]-1))
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
                           between.cohort.wk, eff.schedule, general_ls[[14]]-1))
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
  
  loss.est.pipe<-loss.all(pro, final.assessment.timepoint, mcmc.niter, mcmc.burnin.prop, n.doses, eff,
                     0.5, 0.5, n.sim.final)
  loss.est.bb<-loss.all.bb(pro, final.assessment.timepoint, mcmc.niter, mcmc.burnin.prop, n.doses, eff,
                     0.5, 0.5, n.sim.final)
  final.admiss<-intersect(boin_admiss(target,clin_mat[,2],clin_mat[,4], beta_a, beta_b), futility_decision(min_eff,eff[,2],eff[,5], beta_a, beta_b, n.doses))
  explored<- unique(clin_mat[,2])
  allocated<- sapply(1:n.doses, function (k) nrow(clin_mat[clin_mat[,2]==k,]))
  cens_val<-cbind(as.numeric(pat_cens<=eff.schedule[2]),clin_mat[,2] )
  ncens<-sapply(1:n.doses, function (k) sum(cens_val[cens_val[,2]==k,1]))[explored]
  dlt_val<- cbind(clin_mat[,4],clin_mat[,2])
  ndlt<-sapply(1:n.doses, function (k) sum(dlt_val[dlt_val[,2]==k,1]))[explored]

  return(list(#final.rec=final, 
    boin.admiss=final.admiss, dose.explored=explored, pat.allocated=allocated[explored],
    cdlt = clin_mat,eff=eff, pro=pro, loss_est=loss.est.pipe[[1]][explored], eff_est=loss.est.pipe[[3]][explored], pro_est=loss.est.pipe[[2]][explored], 
    n_cens=ncens,n_dlt= ndlt, loss_est_bb=loss.est.bb[[1]][explored], eff_est_bb= loss.est.bb[[3]][explored]))
}



#TRIAL DESIGN 
trial_design_hypothetical<- function(general_ls, boin_ls, pro_ls, eff_ls, data_after_dlt){
  #general inputs
  n.patient.cohort<- general_ls[[1]]
  first.week.assessed<- general_ls[[2]]
  between.cohort.wk<-general_ls[[3]]
  n.doses<-general_ls[[4]]
  n.timepoints<-general_ls[[5]]
  n.toxicity<- general_ls[[6]]
  n.grades<- general_ls[[7]]
  final.assessment.timepoint<-general_ls[[8]]
  mcmc.niter<-general_ls[[9]]
  mcmc.burnin.prop<- general_ls[[10]]
  n.cohorts.all<- general_ls[[11]]
  n.sim.final<- general_ls[[12]]
  cop.corr<- general_ls[[13]]
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
  max.pro<- pro_ls[[5]]
  
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
  median_survival_month<- (general_ls[[14]]-1)
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
      eff<- eff_sim(n.patient.cohort, interim_complete_cohort1, eff_rates, clin_mat[(1:interim_complete_cohort1)*n.patient.cohort,2], between.cohort.wk, eff.schedule, general_ls[[14]]-1)
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
          #eff[pat_rm,c(3,4,5)]<- rep(0, times=3)
          eff<- eff[-pat_rm,]
        }
      }
      eff_admiss<-futility_decision(min_eff, eff[,2], eff[,5], beta_a, beta_b, n.doses)
    }
    if(max(clin_mat[,1])/n.patient.cohort==cohort_allot_interim2){
      new<- eff_sim(n.patient.cohort, cohort_allot_interim2-cohort_allot_interim1, eff_rates, clin_mat[((interim_complete_cohort1+1):(interim_complete_cohort2))*n.patient.cohort,2],
                    between.cohort.wk, eff.schedule, general_ls[[14]]-1)
      relabel<- c(eff[,1], seq(from=(max(eff[,1])+1), by=1, length.out=nrow(new)))
      eff<- rbind(eff,new)
      eff[,1]<- relabel
      week_trunc<-pat_cens[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2)]
      for(k in (n.patient.cohort*interim_complete_cohort1+1):max(eff[,1])){
        i<- k-n.patient.cohort*interim_complete_cohort1
        if(week_trunc[i]<eff.schedule[2]){
          if(week_trunc[i]>=eff.schedule[1]+1 & week_trunc[i]<eff.schedule[2]){
            eff[which(eff[,1]==k), 5]<-eff[which(eff[,1]==k),3] 
            eff[which(eff[,1]==k), 4]<-0
          }else{
            eff[which(eff[,1]==k), 5]<-0
            eff[which(eff[,1]==k), 4]<-0
            eff[which(eff[,1]==k), 3]<-0
          }
        }
      }
      if(data_after_dlt==FALSE && sum(clin_mat[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2),4])>0){
        eff_dlt<-which(clin_mat[(n.patient.cohort*interim_complete_cohort1+1):(n.patient.cohort*interim_complete_cohort2),4]==1)
        pat_rm<- which(eff[,1] %in% (eff_dlt+(n.patient.cohort*interim_complete_cohort1)))
        if(length(pat_rm)>0){
          eff<- eff[-pat_rm,]
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
  new<- eff_sim(n.patient.cohort, length((interim_complete_cohort2+1):n.cohorts.all), eff_rates, c_dose[(interim_complete_cohort2+1):n.cohorts.all],
                between.cohort.wk, eff.schedule, general_ls[[14]]-1)
  relabel<- c(eff[,1], seq(from=(max(eff[,1])+1), by=1, length.out=nrow(new)))
  eff<- rbind(eff,new)
  eff[,1]<- relabel
  week_trunc<-pat_cens[(n.patient.cohort*interim_complete_cohort2+1):(n.patient.cohort*n.cohorts.all)]
  for(k in (n.patient.cohort*interim_complete_cohort2+1):max(eff[,1])){
    i<- k-n.patient.cohort*interim_complete_cohort2
    if(week_trunc[i]<eff.schedule[2]){
      if(week_trunc[i]>=eff.schedule[1]+1 & week_trunc[i]<eff.schedule[2]){
        eff[which(eff[,1]==k), 5]<-eff[which(eff[,1]==k),3] 
        eff[which(eff[,1]==k), 4]<-0
      }else{
        eff[which(eff[,1]==k), 5]<-0
        eff[which(eff[,1]==k), 4]<-0
        eff[which(eff[,1]==k), 3]<-0
        print(k)
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
      eff<- eff[-pat_rm,]
    }
  }
  #find final recommendation
  loss.est.pipe<-loss.all(pro, final.assessment.timepoint, mcmc.niter, mcmc.burnin.prop, n.doses, eff,
                          0.5, 0.5, n.sim.final)
  loss.est.bb<-loss.all.bb(pro, final.assessment.timepoint, mcmc.niter, mcmc.burnin.prop, n.doses, eff,
                           0.5, 0.5, n.sim.final)
  final.admiss<-intersect(boin_admiss(target,clin_mat[,2],clin_mat[,4], beta_a, beta_b), futility_decision(min_eff,eff[,2],eff[,5], beta_a, beta_b, n.doses))
  explored<- unique(clin_mat[,2])
  allocated<- sapply(1:n.doses, function (k) nrow(clin_mat[clin_mat[,2]==k,]))
  cens_val<-cbind(as.numeric(pat_cens<=eff.schedule[2]),clin_mat[,2] )
  ncens<-sapply(1:n.doses, function (k) sum(cens_val[cens_val[,2]==k,1]))[explored]
  dlt_val<- cbind(clin_mat[,4],clin_mat[,2])
  ndlt<-sapply(1:n.doses, function (k) sum(dlt_val[dlt_val[,2]==k,1]))[explored]
  
  return(list(#final.rec=final, 
    boin.admiss=final.admiss, dose.explored=explored, pat.allocated=allocated[explored],
    cdlt = clin_mat,eff=eff, pro=pro, loss_est=loss.est.pipe[[1]][explored], eff_est=loss.est.pipe[[3]][explored], pro_est=loss.est.pipe[[2]][explored], 
    n_cens=ncens,n_dlt= ndlt, loss_est_bb=loss.est.bb[[1]][explored], eff_est_bb= loss.est.bb[[3]][explored]))
}
