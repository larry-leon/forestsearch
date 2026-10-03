title(sub="(b) Estrogen > 0")
with(rotterdam,table(hormon,chemo))
gbsg_validate <-within(rotterdam,{
rfstime <- ifelse(recur==1,rtime,dtime)
#rfstime <- pmin(rtime,dtime)
# censor if beyond 84 months
t_months<-rfstime/30.4375
time_months <- pmin(t_months,84)
#time_months <- t_months
#status <- ifelse((recur==1 | death==1) & t_months <= 84,1,0)
status <- pmax(recur,death)
ignore <- (recur==0 & death==1 & rtime < dtime)
status2 <- ifelse(recur==1 | ignore, recur, death)
rfstime2 <- ifelse(recur==1 | ignore, rtime,dtime)
time_months2 <- rfstime2/30.4375
grade3 <- ifelse(grade=="3",1,0)
treat <- hormon
id<-as.numeric(c(1:nrow(rotterdam)))
SG0 <- ifelse(er<=0,0,1)
})
# Check subject 40 who was censored at recurrence but death afterwards
# Event at 6.6 years
subset(gbsg_validate,pid==40)
tte.name<-c("time_months")
event.name<-c("status")
treat.name<-c("treat")
# nodes >=1 to match gbsg?
df.analysis <- subset(gbsg_validate, nodes >= 1)
# PS weighting
glm_covs <- glm(treat ~ age+meno+size+grade3+nodes+pgr+er+chemo, data=df.analysis, family="binomial")
ps_covs <- glm_covsfitted.values
# nodes >=1 to match gbsg?
df.analysis <- subset(gbsg_validate, nodes >= 1)
# PS weighting
glm_covs <- glm(treat ~ age+meno+size+grade3+nodes+pgr+er+chemo, data=df.analysis, family="binomial")
ps_covs <- glm_covs$fitted.values
glm_null <- glm(treat ~ 1, data=df.analysis, family="binomial")
ps_null <- glm_null$fitted.values
# Stabilized weight
df.analysis$weight <- ifelse(df.analysis$treat==1, ps_null/ps_covs, (1-ps_null)/(1-ps_covs))
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.85
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-1,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.035,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.035,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.85
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-1,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.035,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.035,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.10, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.15, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.75
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-1,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.75
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-2,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3), time.zero = -1,
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3), time.zero = -2,
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.65
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-2,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
gbsg_validate <-within(rotterdam,{
rfstime <- ifelse(recur==1,rtime,dtime)
#rfstime <- pmin(rtime,dtime)
# censor if beyond 84 months
t_months<-rfstime/30.4375
time_months <- pmin(t_months,84)
#time_months <- t_months
#status <- ifelse((recur==1 | death==1) & t_months <= 84,1,0)
status <- pmax(recur,death)
# Second version (conservative per Therneau)
ignore <- (recur==0 & death==1 & rtime < dtime)
status2 <- ifelse(recur==1 | ignore, recur, death)
rfstime2 <- ifelse(recur==1 | ignore, rtime,dtime)
time_months2 <- rfstime2/30.4375
grade3 <- ifelse(grade=="3",1,0)
treat <- hormon
id<-as.numeric(c(1:nrow(rotterdam)))
SG0 <- ifelse(er<=0,0,1)
})
# Check subject 40 who was censored at recurrence but death afterwards
# Event at 6.6 years
#subset(gbsg_validate,pid==40)
tte.name<-c("time_months")
event.name<-c("status")
treat.name<-c("treat")
# nodes >=1 to match gbsg?
df.analysis <- subset(gbsg_validate, nodes >= 1)
# PS weighting (stabilized)
glm_covs <- glm(treat ~ age+meno+size+grade3+nodes+pgr+er+chemo, data=df.analysis, family="binomial")
ps_covs <- glm_covs$fitted.values
glm_null <- glm(treat ~ 1, data=df.analysis, family="binomial")
ps_null <- glm_null$fitted.values
# Stabilized weight
df.analysis$weight <- ifelse(df.analysis$treat==1, ps_null/ps_covs, (1-ps_null)/(1-ps_covs))
par(mfrow=c(1,2))
byrisk <- 12
risk.cex<-0.65
legend.cex<-0.90
# Add all events
evs<-sort(unique(df.analysis$time_months[which(df.analysis$status==1)]))
tpoints.add<-c(-1,evs,max(df.analysis$time_months))
dfp <- subset(df.analysis, er<=0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
stop.onerror=TRUE,check.KM=TRUE,
Xlab="Months",Ylab="Recurrence-free Survival",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(a) Estrogen = 0")
dfp <- subset(df.analysis, er>0)
kmH.fit<-KM.plot.2sample.weighted(Y=dfp[,c(tte.name)],E=dfp[,c(event.name)],Treat=dfp[,c(treat.name)],
Weight=dfp[,c("weight")],
risk.set=TRUE,by.risk=byrisk,tpoints.add=tpoints.add,risk.cex=risk.cex,
show.Y.axis=FALSE,
stop.onerror=TRUE, check.KM=TRUE,
Xlab="Months",Ylab="",details=FALSE,
show.ticks=TRUE,
col.1="black", col.2="darkgrey",
ltys=c(1,2),lwds=c(3,3),
plotit=TRUE, risk_offset=0.125, risk_delta=0.05,
show.logrank=FALSE,show.med=FALSE,show.cox=TRUE)
title(sub="(b) Estrogen > 0")
library(speff2trial)
help(actg175)
help(ACTG175)
with(ACTG175,quantile(preanti,0.75))
df <- with(ACTG175, arms %in% c(0,1))
dim(df)
df <- subset(ACTG175, arms %in% c(0,1))
dim(df)
with(df,quantile(preanti,0.75))
summary(df)
df <- subset(ACTG175, arms %in% c(3,1))
dim(df)
summary(df$preanti)
summary(df$age)
1/0.8
df_check <- subset(ACTG175,preanti <= 744.5)
library(survival)
df_check <- subset(ACTG175,preanti <= 744.5)
library(speff2trial)
df_check <- subset(ACTG175,preanti <= 744.5)
summary(df_check$z30)
df_check <- subset(ACTG175,preanti >0 & preanti <= 744.5)
summary(df_check$z30)
table(df_check$oprior)
help(ACTG175)
30/737
nrow(ACTG175)
2*191/1083
2*191
cm <- citation("glmnet")
cm
cm <- citation("policytree")
cm
print(cm, bibtex=TRUE)
vt.tree
libary(aVirtualtwins)
libary(aVirtualTwins)
library(aVirtualTwins)
vt.tree
help(vt.tree)
help(rpart)
rm(list=ls())
library(doFuture)
library(doRNG)
registerDoFuture()
registerDoRNG()
# workers=48 works well on popOS
plan("multisession")
library(randomForestSRC)
dim(shf)
data(peakVO2, package = "randomForestSRC")
df.analysis <- peakVO2
confounders.name <- names(df.analysis)
loc_remove <- which(confounders.name %in% c("died","ttodead","aspirin"))
confounders.name <- confounders.name[-c(loc_remove)]
dim(df.analysis)
length(confounders.name)
help(peakVO2)
help(cite)
cite("R, randomForestSRC",refs,textual=TRUE)
revs <- NULL
refs <- NULL
cite("R, randomForestSRC",refs,textual=TRUE)
refs
refs <- "1"
cite("R, randomForestSRC",refs,textual=TRUE)
cm <- citation(package = "randomForestSRC")
print(cm, bibentry=TRUE)
print(cm, bibtex=TRUE)
cm <- citation(package="speff2trial")
print(cm, bibtex=TRUE)
install.packages("~/Desktop/Link to GitHub/htesim_0.0.1.tar.gz", repos = NULL, type = "source")
install.packages("tram")
install.packages("~/Desktop/Link to GitHub/htesim_0.0.1.tar.gz", repos = NULL, type = "source")
library(htesim)
help(htesim)
dgp
sanitize_fct
sanitize
load("results-input//pAnyH_approx_nsg_nocensoring.Rdata")
mdd.60<-min(hrH.plims[pAnyH.approx.60>=0.80])
mdd.80<-min(hrH.plims[pAnyH.approx.80>=0.80])
mdd.100<-min(hrH.plims[pAnyH.approx.100>=0.80])
# e for type-1 error
Edd.60<-max(hrH.plims[pAnyH.approx.60<=0.10])
Edd.80<-max(hrH.plims[pAnyH.approx.80<=0.10])
Edd.100<-max(hrH.plims[pAnyH.approx.100<=0.10])
# Look at thetaH around 0.80 (considering beneficial)
loc.beni <- min(which(hrH.plims>=0.75))
beni_hr <- round(hrH.plims[loc.beni],2)
