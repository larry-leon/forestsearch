
rm(list=ls())
# Local github
codepath<-c("/media/larryleon/My Projects/GitHub/Forest-Search/R/")
source(paste0(codepath,"source_forestsearch_v0.R"))
source_fs_functions(file_loc=codepath)

library(kableExtra)
library(knitr)
library(ggplot2)
library(gridExtra)
library(cubature)
library(aVirtualTwins)
library(randomForest)
library(survival)
library(survminer)
library(grf)
library(policytree)
library(data.table)
library(plyr)
library(dplyr)
library(glmnet)
library(corrplot)

maxFollow<-84
cens.type<-"weibull"

# m1 -censoring adjustment
muC.adj<-log(1.5)
k.z3<-1.0
k.treat<-0.9

z1_frac<-0.25 # Default model index 'm1' (The 1st quartile of z1=er)
pH_super<-0.125 # non-NULL re-defines z1_frac

if(is.null(pH_super)){
#pH_check<-with(gbsg,mean(pgr<=quantile(pgr,c(z3_frac),1,0) & er<=quantile(er,z1_frac)))
pH_check<-with(gbsg,mean(meno==0 & er<=quantile(er,z1_frac)))
cat("Underlying pH_super",c(pH_check),"\n")
}
# pH_super specified
# If pH_super then override  z1_frac and find z1_frac to yield pH_super

if(!is.null(pH_super)){
  # Approximate Z1 quantile to yield pH proportion  
  z1_q<-uniroot(propH.obj4,c(0,1),tol=0.0001,pH.target=pH_super)$root
  #pH_check<-with(gbsg,mean(pgr<=quantile(pgr,c(z3_frac),1,0) & er<=quantile(er,z1_q)))
  pH_check<-with(gbsg,mean(meno==0 & er<=quantile(er,z1_q)))
  cat("pH",c(pH_check),"\n")
  rel_error<-(pH_super-pH_check)/pH_super
  if(abs(rel_error)>=0.1) stop("pH_super approximation relative error exceeds 10%")
  z1_frac<-z1_q
  cat("Underlying pH_super",c(pH_check),"\n")
  }

# Bootstrap on log(hr) scale converted to HR (est.loghr=TRUE & est.scale="hr")
t.start.all<-proc.time()[3]

#########################
# Forest search criteria
#########################
hr.threshold<-1.25   # Initital candidates 
hr.consistency<-1.0  # Candidates for many splits
pconsistency.threshold <- 0.9
stop.threshold <- 0.95
maxk<-2
nmin.fs<-60
pstop_futile<-0.5
# Limit timing for forestsearch
max.minutes<-3.0
m1.threshold<-Inf # Turning this off (Default)
#pconsistency.threshold<-0.70 # Minimum threshold (will choose max among subgroups satisfying)
fs.splits<-400 # How many times to split for consistency
# vi is % factor is selected in cross-validation --> higher more important
vi.grf.min<-0.2
# Null, turns off grf screening
d.min<-10 # Min number of events for both arms (d0.min=d1.min=d.min)
# default=5
##########################
# Virtual twins analysis
##########################
# Counter-factual difference (C-E) >= vt.threshold
# Large values in favor of C (control)
vt.threshold<-0.225  # For VT delta
treat.threshold<-0.0

maxdepth<-2
n.min<-60
ntree<-1000
# GRF criteria
dmin.grf<-12.0 # For GRF delta
# Note: For CRT this represents dmin.grf/2 RMS for control (-dmin.grf/2 for treatment)
frac.tau<-0.60

outcome.name<-c("y.sim")
event.name<-c("event.sim")
id.name<-c("id")
treat.name<-c("treat")

cox.formula.sim<-as.formula(paste("Surv(y.sim,event.sim)~treat"))
cox.formula.adj.sim<-as.formula(paste("Surv(y.sim,event.sim)~treat+v1+v2+v3+v4+v5"))

mod.harm<-"alt"
hrH.target <- 2.5
# out.loc = NULL turns off file creation
N <- 700
this.dgm<-get.dgm4(mod.harm=mod.harm,N=N,k.treat=k.treat,
hrH.target=hrH.target,cens.type=cens.type,out.loc=NULL,details=TRUE,parms_torand=FALSE)

n_add_noise <- 3

dgm<-this.dgm$dgm

Nsims <-100

SimsToLook <- matrix(NA, nrow=Nsims, ncol=5)

for(ss in 1:Nsims){
sim <- ss
x<-sim_aftm4_gbsg(dgm=dgm,n=N,maxFollow=maxFollow,muC.adj=muC.adj,simid=sim)

kmfit <- survfit(Surv(y.sim,event) ~ treat, data=x)

ggsurvplot(kmfit, data=x, main="K-M curves for simulated data",
legend="top", legend.title="Treatment",
legend.labs=c("Control","Experimental"),
palette="grey", risk.table=TRUE, risk.table.col="strata")

  if(n_add_noise==0){
  confounders.name <- c("z1","z2","z3","z4","z5","size","grade")
  }
    
  if(n_add_noise==5){
  set.seed(8316951+1000*sim)
  # Add 5 noise 
  x$noise1 <- rnorm(N)
  x$noise2 <- rnorm(N)
  x$noise3 <- rnorm(N)
  x$noise4 <- rnorm(N)
  x$noise5 <- rnorm(N)
  confounders.name <- c("z1","z2","z3","z4","z5","size","grade","noise1","noise2","noise3","noise4","noise5")
  }
  
  if(n_add_noise==3){
    set.seed(8316951+1000*sim)
    # Add 3 noise 
    x$noise1 <- rnorm(N, sd=1)
    x$noise2 <- rnorm(N, sd=1)
    x$noise3 <- rnorm(N, sd=1)
    confounders.name <- c("z1","z2","z3","z4","z5","size","grade","noise1","noise2","noise3")
  }

# More challenging
#which_replace <- which(confounders.name=="z1")
#confounders.name[which_replace] <- "er"

cox.formula.check<-as.formula(paste("Surv(y.sim,event.sim)~treat+v1+v2+v3+v4+v5+size+grade+noise1+noise2+noise3"))
coxph(cox.formula.check,data=x)

Zm <- cor(as.matrix(x[,c(confounders.name)]))
corrplot(Zm)

# Options
# Allconfounders.name is list of confounders
# within analysis dataset
# (1) use_lasso=TRUE & use_grf=FALSE
# Lasso used to possibly reduce dimension
# Any continuous factors are cut at medians
# (2) use_lasso=TRUE & use_grf=TRUE
# Lasso used to reduce dimension
# Continuous covariates are cut at medians
# However, if GRF selects a covariate cut
# then only that cut is used:
# For example if "age <= median(ag)" is 
# called for per Lasso but GRF includes
# "age <= 54", then only the latter is used
# (3) use_grf_only = TRUE (overrides use_lasso and use_grf)
# Only factors selected via GRF are used
# (4) use_lasso = F & use_grf =T
# Median cuts (unless selected via GRF) as  in (2)
# However no possible dimension reduction via lasso
# All categorical factors included

use_lasso <- TRUE
use_grf <- TRUE
use_grf_only <- FALSE

fs.est <- forestsearch(df.analysis=x, Allconfounders.name=confounders.name,
details=TRUE,use_lasso=use_lasso, use_grf=use_grf, use_grf_only=use_grf_only,
dmin.grf=12, frac.tau=0.6,
conf_force=NULL, outcome.name=outcome.name,treat.name=treat.name,
event.name=event.name,id.name=id.name,n.min=nmin.fs,hr.threshold=hr.threshold,
hr.consistency=hr.consistency,fs.splits=fs.splits,d0.min=d.min,d1.min=d.min,
pstop_futile=pstop_futile,pconsistency.threshold=pconsistency.threshold, 
stop.threshold=stop.threshold,max.minutes=max.minutes,maxk=maxk,by.risk=12,
plot.sg=TRUE,vi.grf.min=vi.grf.min)

if(!is.null(fs.est$sg.harm)){
dfH <- subset(fs.est$df.est, flag.harm==1)
dfH_est <- subset(fs.est$df.est, treat.recommend==0)
dfSens <- subset(fs.est$df.est, flag.harm==1 & treat.recommend==0)
sensH <- nrow(dfSens)/nrow(dfH)
ppvH <- nrow(dfSens)/nrow(dfH_est)

dfHc <- subset(fs.est$df.est, flag.harm==0)
dfHc_est <- subset(fs.est$df.est, treat.recommend==1)
dfSensc <- subset(fs.est$df.est, flag.harm==0 & treat.recommend==1)
sensHc <- nrow(dfSensc)/nrow(dfHc)
ppvHc <- nrow(dfSensc)/nrow(dfHc_est)

SimsToLook[ss,] <- c(ss, sensH, ppvH, sensHc, ppvHc)
}

cat("Simulation sensitivity=",c(ss,sensH,ppvH,sensHc,ppvHc),"\n")
}

sims_look <- as.data.table(SimsToLook)
names(sims_look) <- c("sim","sensH","ppvH","sensHc","ppvHc")

summary(sims_look)

sims_look <- na.omit(sims_look)

tolook <- subset(sims_look, ppvHc >= 0.7 & ppvHc <1)

print(tolook)


#sim     ppv_H
#1:   8 0.7115385
#2:  41 0.7631579
#3:  46 0.7037037
#4:  50 0.7699115
#5:  52 0.7281553
#6:  58 0.7142857
#7:  64 0.8160000
#8:  71 0.8080808
#9:  75 0.8100000
#10:  79 0.7596154
#11:  80 0.7757009
#12:  86 0.7818182
#13:  93 0.8017241
#14:  98 0.7818182


#sim      sensH       ppvH    sensHc     ppvHc
#1:  11 0.50000000 0.57317073 0.9422442 0.9239482
#2:  14 0.68918919 0.56666667 0.9376997 0.9622951
#3:  17 0.32967033 0.28301887 0.8752053 0.8973064
#4:  18 0.44897959 0.49438202 0.9252492 0.9116203
#5:  20 0.80459770 0.43478261 0.8515498 0.9684601
#6:  23 0.81012658 0.52032520 0.9049919 0.9740035
#7:  30 0.78021978 0.51824818 0.8916256 0.9644760
#8:  32 0.65934066 0.58823529 0.9310345 0.9481605
#9:  36 0.32142857 0.38571429 0.9301948 0.9095238
#10:  39 0.56179775 0.53191489 0.9279869 0.9356436
#11:  61 0.39759036 0.41772152 0.9254457 0.9194847
#12:  63 0.51351351 0.50666667 0.9408946 0.9424000
#13:  75 0.54320988 0.53012048 0.9369952 0.9400324
#14:  81 0.47872340 0.53571429 0.9356436 0.9204545
#15:  87 0.61111111 0.44444444 0.9124204 0.9534110
#16:  91 0.51724138 0.63380282 0.9575856 0.9332273
#17:  97 0.51923077 0.60674157 0.9412752 0.9181669
#18:  98 0.06976744 0.08450704 0.8941368 0.8728140
#19:  99 0.24358974 0.15702479 0.8360129 0.8981002
#20: 100 0.48101266 0.47500000 0.9323671 0.9338710




