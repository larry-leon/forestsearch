
rm(list=ls())

# Local github popOS
codepath<-c("/media/larryleon/My Projects/GitHub/Forest-Search/R/")

# On MAC
#codepath <- c("/Users/larryleon/Documents/GitHub/Forest-Search/R/")

source(paste0(codepath,"source_forestsearch_v0.R"))
source_fs_functions(file_loc=codepath)

library(kableExtra)
library(knitr)
library(ggplot2)
library(gridExtra)
library(randomForest)
library(survival)
library(grf)
library(policytree)
library(data.table)
library(plyr)
library(dplyr)
library(glmnet)
library(cli)
library(gtsummary)

load(c("applications/simulated_ex1/output/sim-FS4_N1k_Noise=3_hrH=2_sim99.Rdata"))

df.analysis <- x

confounders.name <- c(fs.est$Allconfounders.name)

outcome.name<-c("y.sim")
event.name<-c("event.sim")
id.name<-c("id")
treat.name<-c("treat")

itt_tab<-SGtab(df=df.analysis,SG_flag="ITT",outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)
itt_tab$res_out

cat("True H and Hc marginal hazard ratios",c(round(c(dgm$hr.H.true,dgm$hr.Hc.true),2)),"\n")

# Limit timing for forestsearch

max.minutes<-1.0

set.seed(8316952)

dfa <- fs.est$df.predict[,c(confounders.name,outcome.name,event.name,id.name,treat.name)]

dfa$treat.recommend.original <- fs.est$df.predict[,c("treat.recommend")]

SG_tab<-SGtab(df=dfa,SG_flag="treat.recommend.original",sg1_name="Recommended",sg0_name="Not Recommended",
              outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)
SG_tab$res_out

# Compare to GRF predictions
X<-as.matrix(dfa[,confounders.name])
# Convert to numeric
X<-apply(X,2,as.numeric)
Y<-dfa[,outcome.name]
W<-dfa[,treat.name]
D<-dfa[,event.name]
tau.rmst<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

cs.forest <- causal_survival_forest(X,Y,W,D,horizon=0.6*tau.rmst,seed=8316951)
tau.hat <- predict(cs.forest)$predictions

dfnew <- as.data.frame(dfa)
dfnew$tauhat.grf <- c(tau.hat)

#Kfolds <- round(2*nrow(dfnew)/3,0)

Kfolds <- 10

set.seed(8316953)

df_scrambled <- dfnew[sample(nrow(dfnew)),]
folds <- cut(seq(1,nrow(dfnew)),breaks=Kfolds,labels=FALSE)

cat("Range of unique left-out fold sample sizes",c(range(unique(table(folds)))),"\n")

t.start<-proc.time()[3]

#temp <- cv_forparallel(cv_index=5)

library(doRNG)
library(doFuture)

registerDoFuture()
registerDoRNG()
plan("multisession", workers=124)

resCV <- foreach(
  cv_index = seq_len(Kfolds),
  .options.future=list(seed=TRUE),
  .combine="rbind",
  .errorhandling="pass"
) %dofuture% {
  ans <- cv_forparallel(cv_index)
}

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60

cat("Minutes for Cross-validation",c(Kfolds,t.min),"\n")
cat("Projection per 100",c(t.min*(100/Kfolds)),"\n")

# Extract sg1 and sg2

sg1 <- sg2 <- rep(NA,Kfolds)
for(kk in 1:Kfolds){
  # First element since all identical for same cvindex (=kk)  
  sg1[kk] <- subset(resCV, cvindex==kk)$sg1[1] 
  sg2[kk] <- subset(resCV, cvindex==kk)$sg2[1] 
}
SGs_found <- cbind(sg1,sg2)
temp <- CV_sgs(sg1=sg1, sg2=sg2, confs=fs.est$Allconfounders.name, sg_analysis=fs.est$sg.harm)

# Note: Current code does not handle this scenario {z3 <= 0} vs {z3};
# Treats them as different!
# Needs revision ...

cat("Any found",c(mean(temp$any_found)),"\n")
cat("Exact match",c(mean(temp$exact_match)),"\n")
cat("At least 1 match",c(mean(temp$one_match)),"\n")
cat("Cov 1 any",c(mean(temp$cov1_any)),"\n")
cat("Cov 2 any",c(mean(temp$cov2_any)),"\n")

CV_summary <- temp

# Propn agreement in H and H^c
# between sample estimate (original) and cross-validation
tabit <- with(resCV,table(treat.recommend,treat.recommend.original))

sens_H <- tabit[1,1]/sum(tabit[,1])
sens_Hc <- tabit[2,2]/sum(tabit[,2])
ppv_H <- tabit[1,1]/sum(tabit[1,])
ppv_Hc <- tabit[2,2]/sum(tabit[2,])

cat("Agreement (sens, ppv) in H and Hc:",c(sens_H,sens_Hc,ppv_H,ppv_Hc),"\n")

df.test.out <- resCV

# GRF
df.test.out$treat.agree_grf <- ifelse(df.test.out$treat.recommend==0 & df.test.out$tauhat.grf <=0,
                                      2,ifelse(df.test.out$treat.recommend==1 & df.test.out$tauhat.grf >0, 1, 3))

df.test.out$treat_grf <- ifelse(df.test.out$tauhat.grf>0,1,0)
df.test.out$Notreat_grf <- 1-df.test.out$treat_grf

tbl_est <- 
  df.test.out %>% 
  select(treat.recommend,treat.recommend.original,tauhat.grf,treat.agree_grf,treat_grf,Notreat_grf) %>%
  tbl_summary(
    by= treat.recommend,
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c("{mean} ({p25}, {p75})", "{min}, {max}"),
    missing="no"
  )
tbl_est

tabit <- with(df.test.out,table(treat.recommend,treat_grf))

sens_H <- tabit[1,1]/sum(tabit[,1])
sens_Hc <- tabit[2,2]/sum(tabit[,2])
ppv_H <- tabit[1,1]/sum(tabit[1,])
ppv_Hc <- tabit[2,2]/sum(tabit[2,])
cat("Agreement (sens, ppv) in H and Hc:",c(sens_H,sens_Hc,ppv_H,ppv_Hc),"\n")

cat("Range of unique left-out fold sample sizes",c(range(unique(table(df.test.out$cvindex)))),"\n")

cat("K-fold aggregated sample size = n?",c(length(unique(df.test.out$id)),nrow(dfnew)),"\n")

# Compare with raw-unadjusted originals

SG_tab<-SGtab(df=as.data.frame(df.test.out),SG_flag="treat.recommend.original",sg1_name="Recommended",sg0_name="Not Recommended",
              outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)

print(SG_tab$res_out)

SG_tab_Kfold<-SGtab(df=as.data.frame(df.test.out),SG_flag="treat.recommend",sg1_name="Recommended",sg0_name="Not Recommended",
                    outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,draws=0)

print(SG_tab_Kfold$res_out)

save(SG_tab_Kfold,SGs_found,CV_summary,Kfolds,df.test.out,H_est,Hc_est,file="applications/simlated_ex1/output/sim-FS4_Sim=99_crossvalidation_Kfolds=twothirdsN_v0.Rdata")

