
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

library(cli)

library(corrplot)

library(table1)

load(c("output/sim-FS4_N1k_Noise=3_hrH=2_sim99.Rdata"))

# k-fold splitting
#Randomly shuffle the data
#yourdata<-yourdata[sample(nrow(yourdata)),]
#Create 10 equally size folds

# Limit timing for forestsearch
max.minutes<-3.0

outcome.name<-c("y.sim")
event.name<-c("event.sim")
id.name<-c("id")
treat.name<-c("treat")

conf_force <- NULL

set.seed(8316951)

xx <- x[sample(nrow(x)),]
folds <- cut(seq(1,nrow(xx)),breaks=10,labels=FALSE)

X<-as.matrix(xx[,c("z1","z2","z3","z4","z5","size")])
# Convert to numeric
X<-apply(X,2,as.numeric)
Y<-xx[,outcome.name]
W<-xx[,treat.name]
D<-xx[,event.name]
tau.rmst<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

cs.forest <- causal_survival_forest(X,Y,W,D,horizon=0.8*tau.rmst,seed=8316951)

sg1<-rep(NA,10)
sg2<-rep(NA,10)
df.test.out <- NULL

for(ii in 1:10){
testIndexes <- which(folds==ii,arr.ind=TRUE)
x.test <- xx[testIndexes, ]
x.train <- xx[-testIndexes, ]

# GRF 
X.test<-as.matrix(x.test[,c("z1","z2","z3","z4","z5","size")])
X.test<-apply(X.test,2,as.numeric)

tau.hat <- predict(cs.forest, X.test)$predictions
x.test$tauhat.grf <- tau.hat

fs.train <- forestsearch(df.analysis=x.train, df.predict=x.test,
Allconfounders.name=fs.est$Allconfounders.name,
details=TRUE,use_lasso=fs.est$use_lasso, use_grf=fs.est$use_grf, use_grf_only=fs.est$use_grf_only,
dmin.grf=fs.est$dmin.grf, frac.tau=fs.est$frac.tau,
conf_force=conf_force, outcome.name=outcome.name,treat.name=treat.name,
event.name=event.name,id.name=id.name,
n.min=fs.est$n.min,hr.threshold=fs.est$hr.threshold,hr.consistency=fs.est$hr.consistency,
fs.splits=fs.est$fs.splits,d0.min=fs.est$d0.min,d1.min=fs.est$d1.min,
pstop_futile=fs.est$pstop_futile,pconsistency.threshold=fs.est$pconsistency.threshold, 
stop.threshold=fs.est$stop.threshold,
max.minutes=max.minutes,
maxk=fs.est$maxk,by.risk=12,
plot.sg=TRUE)

sg1[ii] <- fs.train$sg.harm[1]
sg2[ii] <- fs.train$sg.harm[2]

df.test <- fs.train$df.predict

temp <- df.test[,c("id","z1","z2","z3","z4","z5","size","hlin.ratio","treat.recommend","tauhat.grf")]

df.test.out <- rbind(df.test.out,temp)

rm("temp")

mean(df.test$hlin.ratio)
with(df.test,mean(h1.potential/h0.potential))

# Estimated harm group
with(df.test,plot(age,hlin.ratio,ylim=c(-0.5,3)))
mean(subset(df.test,treat.recommend==0)$hlin.ratio)
mean(subset(df.test,treat.recommend==1)$hlin.ratio)
}
cbind(sg1,sg2)

df.test.out$treat.agree <- ifelse(df.test.out$treat.recommend==0 & df.test.out$tauhat.grf <=0,
2,ifelse(df.test.out$treat.recommend==1 & df.test.out$tauhat.grf >0, 1, 3))

plot(df.test.out$tauhat.grf,df.test.out$id,col=df.test.out$treat.agree)

#plot(df.test.out$hlin.ratio,df.test.out$id,col=df.test.out$treat.recommend+1)
#legend("none",c("control","treat"),col=c(1,2), pch=19, bty="n")

mean(subset(df.test.out,treat.recommend==0)$hlin.ratio)
mean(subset(df.test.out,treat.recommend==0)$tauhat.grf)

mean(subset(df.test.out,treat.recommend==1)$hlin.ratio)
mean(subset(df.test.out,treat.recommend==1)$tauhat.grf)



