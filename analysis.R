setwd("C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/PRO_ADD_resubmission/PRO-ADD")
#loss for each scenario

shape <- read.csv("shape_param_inc.csv")[,-1]
rate <- read.csv("rate_param_inc.csv")[1,2]

matrix<- matrix(nrow=5, ncol=9)
for(i in 1:5){
  matrix[i,]<-sapply(1:9, function (k) shape[i,k]/(shape[i,k]+rate))
}

pro.sc1.2<- shape
pro.sc3<- shape[c(1:3,3,3),]
pro.sc4<- shape[c(1:3,5,5),]
pro.sc5.6<-shape[,c(1,4,5,6,7,8,9,9,9)]

pro.scen<- list(pro.sc1.2,pro.sc1.2, pro.sc3, pro.sc4, pro.sc5.6, pro.sc5.6)


eff.sc1<-c(0.05, 0.08, 0.24, 0.42, 0.44)
eff.sc2<-c(0.05, 0.08, 0.42, 0.42, 0.42)
eff.sc3<-c(0.05, 0.08, 0.27, 0.28, 0.44)
eff.sc4<-c(0.05, 0.08, 0.27, 0.27, 0.27)
eff.sc5<-c(0.05, 0.08, 0.27, 0.42, 0.44)
eff.sc6<-c(0.05, 0.08, 0.42, 0.42, 0.42)

eff.scen<- list(eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6)


matrix.list<- list()
for(j in 1:6){
  matrix<- matrix(nrow=5, ncol=9)
  for(i in 1:5){
    matrix[i,]<-sapply(1:9, function (k) pro.scen[[j]][i,k]/(pro.scen[[j]][i,k]+rate))
  }
  matrix.list[[j]]<- matrix
}

#Results and recommendations 
loss.target<-sqrt(0^2 + (1-0.1)^2)

loss.mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  loss.mat[i,]<-sqrt(matrix.list[[i]][,9]^2 + (1-eff.scen[[i]])^2)
}

round(loss.mat,2)

mat<- matrix(nrow=6, ncol=5)

################
### MAIN MANUSCRIPT
#Table 3 and 4 - Data available after dlt
mtd<- 5
target<- 0.9

#optimal
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}



for(i in 1:6){
eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
)))
  mat[i,]<-vec/5000
}
round(mat,2)


#acceptable
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('results/hypothetical/hypothetical.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}
round(mat,2)


################
#sensitivity - moderate correlation
mtd<- 5
target<- 0.9

#optimal
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}



for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}
round(mat,2)


#acceptable
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('results/hypothetical.dltresp0.5/hypothetical.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}
round(mat,2)


## looked up until here and updated with the new analysis
##### n allocated to each dose 
n.mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".n[is.na(sc",i,".n)] <- 0")))
  n.mat[i,]<-eval(parse(text=paste0("colMeans(sc",i,".n)")))
}

colMeans(n.mat)



#Table 3 and 4 - Data not available after dlt (sensitivity analysis 1)
mtd<- 3
target<- 0.9

#optimal
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}


for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#acceptable
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}
round(mat,2)


################
### Supplementary materials
#Table S4 and S5 - Data available after dlt

#optimal
mtd<- 3
target<- 0.9

for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.l.bb.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#acceptable
mtd<-3
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)


#Table S4 and S5 - Data not available after dlt (sensitivity analysis 1)
#optimal
mtd<- 3
target<- 0.9

for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.l.bb.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#acceptable
mtd<-3
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.l.bb.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.patrm.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)


#Table S6 and S7 - Data not available after dlt (sensitivity analysis 2)

#iPipe
#optimal
mtd<- 3
target<- 0.9
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}


for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#acceptable
mtd<- 3
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#beta-binomial 
#optimal
mtd<- 3
target<- 0.9
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.l.bb.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}


for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,]==min(sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]))
  if(length(val)>0){
    if(sc",i,".loss[j,val]<=target){
      vec[val]<- vec[val]+1
    }
  }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

#acceptable
mtd<- 3
mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.l.bb.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sens.sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}
for(i in 1:6){
  eval(parse(text=paste0("
vec<-rep(0, times=5)
for(j in 1:length(sc",i,".mtd)){
  val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
  if(length(val)>0){
      vec[val]<- vec[val]+1
    }
}"
  )))
  mat[i,]<-vec/5000
}

round(mat,2)

################# bias for efficacy

eff_mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.e.csv')[,-1]")))
  eff_mat[i,]<-colMeans(e1, na.rm=TRUE)
}
bias_ipipe<-eff_mat-rbind(NULL, eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6)

eff_matbb<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.e.bb.csv')[,-1]")))
  eff_matbb[i,]<-colMeans(e1, na.rm=TRUE)
}
bias_bb<-eff_matbb-rbind(NULL, eff.sc1, eff.sc2, eff.sc3, eff.sc4, eff.sc5, eff.sc6)
round(bias_bb,3)

######## mse for efficacy

eff_mat_pipe_mse<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.e.csv')[,-1]")))
  eval(parse(text=paste0("scen<- eff.sc",i)))
  val1<-na.omit(e1[,1])-scen[1]
  val2<-na.omit(e1[,2])-scen[2]
  val3<-na.omit(e1[,3])-scen[3]
  val4<-na.omit(e1[,4])-scen[4]
  val5<-na.omit(e1[,5])-scen[5]
  eff_mat_pipe_mse[i,]<-sapply(1:5, function (k) eval(parse(text=paste0("(1/length(val",k,"))*sum(val",k,"^2)"))))
}
round(eff_mat_pipe_mse,3)

eff_mat_bb_mse<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.e.bb.csv')[,-1]")))
  eval(parse(text=paste0("scen<- eff.sc",i)))
  val1<-na.omit(e1[,1])-scen[1]
  val2<-na.omit(e1[,2])-scen[2]
  val3<-na.omit(e1[,3])-scen[3]
  val4<-na.omit(e1[,4])-scen[4]
  val5<-na.omit(e1[,5])-scen[5]
  eff_mat_bb_mse[i,]<-sapply(1:5, function (k) eval(parse(text=paste0("(1/length(val",k,"))*sum(val",k,"^2)"))))
}
round(eff_mat_bb_mse,3)

####### bias for PRO-nAE burden score
pro_scen<- sapply(1:6, function (i) matrix.list[[i]][,9])
round(pro_scen,2)
pro_mat<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.p.csv')[,-1]")))
  pro_mat[i,]<-colMeans(e1, na.rm=TRUE)
}
bias_pro<-pro_mat-rbind(NULL, pro_scen[,1], pro_scen[,2], pro_scen[,3], pro_scen[,4], pro_scen[,5], pro_scen[,6])
round(bias_pro,3)

#mse for PRO-nAE burden score
pro_mat_mse<- matrix(nrow=6, ncol=5)
for(i in 1:6){
  eval(parse(text=paste0("e1<- read.csv('M:/PhD/Trial Designs/mayo_clinic/results/sc",i,".5.60.p.csv')[,-1]")))
  eval(parse(text=paste0("scen<- pro_scen[,",i,"]")))
  val1<-na.omit(e1[,1])-scen[1]
  val2<-na.omit(e1[,2])-scen[2]
  val3<-na.omit(e1[,3])-scen[3]
  val4<-na.omit(e1[,4])-scen[4]
  val5<-na.omit(e1[,5])-scen[5]
  pro_mat_mse[i,]<-sapply(1:5, function (k) eval(parse(text=paste0("(1/length(val",k,"))*sum(val",k,"^2)"))))
}
round(eff_mat_pipe_mse,3)

###
#stop early for futility saftey stopping rules
mtd<- 5
target<- 0.9
for(i in 1:6){
  eval(parse(text=paste0("sc",i,".admiss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.a.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".loss <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.l.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".mtd <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.mr.csv')[,-1]")))
  eval(parse(text=paste0("sc",i,".n <- read.csv('C:/Users/ealger/OneDrive - The Institute of Cancer Research/M/PhD/Trial Designs/mayo_clinic/results/sc",i,".",mtd,".60.n.csv')[,-1]")))
  eval(parse(text=paste0("for(i in 1:length(sc",i,".mtd)){if(5-sc",i,".mtd[i]>0){sc",i,".admiss[i,(sc",i,".mtd[i]+1):5]<-NA}}")))
}

early_stop<-c()
for(i in 1:6){
  early_stop[i]<-eval(parse(text=paste0("(5000-nrow(sc",i,".admiss)+length(which(rowSums(sc",i,".admiss, na.rm=TRUE)==0)))/5000")))
}

#stop early because of no loss
loss_stop<-c()
for(i in 1:6){
  none<-0
  eval(parse(text=paste0("
  for(j in 1:length(sc",i,".mtd)){
    if(sum(sc",i,".admiss[j,], na.rm=TRUE)>0){
      val<-which(sc",i,".loss[j,] %in% sc",i,".loss[j,!is.na(sc",i,".admiss[j,])][sc",i,".loss[j,!is.na(sc",i,".admiss[j,])]<=target])
      if(length(val)==0){
        none<- none+1
      }
    }
  }
  loss_stop[i]<-none/5000")))
}

