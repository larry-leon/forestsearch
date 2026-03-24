# Code for Matching

level.cov<-function(x){
l.x<-c(unique(x))
if(length(l.x)>=5){
cat("Covariate unique values=",c(length(l.x)),"\n")
 stop("Attemp to Match on covariate with more than 5 levels: Are you sure?")
}
return(l.x)
}

matching.fit<-function(Vps.names,fit.ps,data,ExactMatch.names=NULL,details=FALSE){
require(Matching)
Vps<-as.matrix(data[,c(Vps.names)])
if(is.null(fit.ps)) stop("PS fit is null")
X<-fit.ps$fitted

Y<-data$Y
Delta<-data$Delta

if(!is.null(ExactMatch.names)){
ExactMatch.covs<-as.matrix(data[,c(ExactMatch.names)])
pE.covs<-dim(ExactMatch.covs)[2]
# matching.type=0 --> get all three; otherwise only 1 of the matching approaches
# Match on up to 5 covariates
if(pE.covs>6) stop("Matching up to 6 covariates currently allowed [See matching.fit function]")
if(pE.covs==1){
W<-c(ExactMatch.covs)
}
if(pE.covs==2){
# First factor to match on
V.EM.1<-ExactMatch.covs[,1]
V.EM.2<-ExactMatch.covs[,2]
l.1<-level.cov(V.EM.1)
l.2<-level.cov(V.EM.2)
Ws<-expand.grid(l.1,l.2)
W<-V.EM.1 # This just initiates where to store
countw<-0
for(www in 1:nrow(Ws)){
thisw<-Ws[www,]
countw<-countw+1
W[which(V.EM.1==Ws[www,1] & V.EM.2==Ws[www,2])]<-countw
}
}
if(pE.covs==3){
# First factor to match on
V.EM.1<-ExactMatch.covs[,1]
V.EM.2<-ExactMatch.covs[,2]
V.EM.3<-ExactMatch.covs[,3]
l.1<-level.cov(V.EM.1)
l.2<-level.cov(V.EM.2)
l.3<-level.cov(V.EM.3)
Ws<-expand.grid(l.1,l.2,l.3)
W<-V.EM.1 # This just initiates where to store
countw<-0
for(www in 1:nrow(Ws)){
thisw<-Ws[www,]
countw<-countw+1
W[which(V.EM.1==Ws[www,1] & V.EM.2==Ws[www,2] & V.EM.3==Ws[www,3])]<-countw
}
}
if(pE.covs==4){
# First factor to match on
V.EM.1<-ExactMatch.covs[,1]
V.EM.2<-ExactMatch.covs[,2]
V.EM.3<-ExactMatch.covs[,3]
V.EM.4<-ExactMatch.covs[,4]
l.1<-level.cov(V.EM.1)
l.2<-level.cov(V.EM.2)
l.3<-level.cov(V.EM.3)
l.4<-level.cov(V.EM.4)
Ws<-expand.grid(l.1,l.2,l.3,l.4)
W<-V.EM.1 # This just initiates where to store
countw<-0
for(www in 1:nrow(Ws)){
thisw<-Ws[www,]
countw<-countw+1
W[which(V.EM.1==Ws[www,1] & V.EM.2==Ws[www,2] & V.EM.3==Ws[www,3] & V.EM.4==Ws[www,4])]<-countw
}
}
if(pE.covs==5){
# First factor to match on
V.EM.1<-ExactMatch.covs[,1]
V.EM.2<-ExactMatch.covs[,2]
V.EM.3<-ExactMatch.covs[,3]
V.EM.4<-ExactMatch.covs[,4]
V.EM.5<-ExactMatch.covs[,5]
l.1<-level.cov(V.EM.1)
l.2<-level.cov(V.EM.2)
l.3<-level.cov(V.EM.3)
l.4<-level.cov(V.EM.4)
l.5<-level.cov(V.EM.5)
Ws<-expand.grid(l.1,l.2,l.3,l.4,l.5)
W<-V.EM.1 # This just initiates where to store
countw<-0
for(www in 1:nrow(Ws)){
thisw<-Ws[www,]
countw<-countw+1
W[which(V.EM.1==Ws[www,1] & V.EM.2==Ws[www,2] & V.EM.3==Ws[www,3] & V.EM.4==Ws[www,4] & V.EM.5==Ws[www,5])]<-countw
}
}
if(pE.covs==6){
# First factor to match on
V.EM.1<-ExactMatch.covs[,1]
V.EM.2<-ExactMatch.covs[,2]
V.EM.3<-ExactMatch.covs[,3]
V.EM.4<-ExactMatch.covs[,4]
V.EM.5<-ExactMatch.covs[,5]
V.EM.6<-ExactMatch.covs[,6]
l.1<-level.cov(V.EM.1)
l.2<-level.cov(V.EM.2)
l.3<-level.cov(V.EM.3)
l.4<-level.cov(V.EM.4)
l.5<-level.cov(V.EM.5)
l.6<-level.cov(V.EM.6)
Ws<-expand.grid(l.1,l.2,l.3,l.4,l.5,l.6)
W<-V.EM.1 # This just initiates where to store
countw<-0
for(www in 1:nrow(Ws)){
thisw<-Ws[www,]
countw<-countw+1
W[which(V.EM.1==Ws[www,1] & V.EM.2==Ws[www,2] & V.EM.3==Ws[www,3] & V.EM.4==Ws[www,4] & V.EM.5==Ws[www,5] & V.EM.6==Ws[www,6])]<-countw
}
}
rr<-Matchby(Y=Y,Tr=Delta,X=X,M=1,ties=TRUE,by=W,print.level=0)
}

if(is.null(ExactMatch.names)){
rr<-Match(Y=Y,Tr=Delta,X=X,M=1,ties=TRUE)
}
# Weighted difference estimator
#cat("How many are matched:",c(rr$nobs),"\n")
#cat("Estimate allowing for ties:",c(rr$est),"\n")
# Verify calculation
#aa<-(Y[rr$index.treated]-Y[rr$index.control])*rr$weights
#mean(aa)/mean(rr$weights)
m.weights<-rr$weights
Y1.matched<-Y[rr$index.treated]
Y0.matched<-Y[rr$index.control]
m1.M<-mean(Y1.matched*m.weights)/mean(m.weights)
m0.M<-mean(Y0.matched*m.weights)/mean(m.weights)
delta.M<-m1.M-m0.M
if(round(delta.M-rr$est,8)!=0) stop("Error in Matching Estimator (BR.TE.fit)")
data0.M1<-rbind(data[rr$index.treated,],data[rr$index.control,])
data0.M1.unique<-rbind(data[unique(rr$index.treated),],data[unique(rr$index.control),])

if(details){
cat("m1,m0,delta=",c(m1.M,m0.M,delta.M),"\n")
cat("##########################################################################################","\n")
cat("Propensity-Score matching","\n")
if(!is.null(ExactMatch.names)){ 
cat("Matching Exactly on Covariates=",c(ExactMatch.names),"\n")
cat("Summary of weights","\n")
print(table(m.weights))
}
data.temp<-as.data.frame(cbind(Y,Delta,Vps))
data.temp$treat<-factor(data.temp$Delta,labels=c("Historical-Control","Experimental"))
label(data.temp$treat)<-c("Study")
data.matched<-data.temp[c(rr$index.treated,rr$index.control),]
# Check matching on PS factors (include outcome Y)
to.check<-as.formula(paste("treat ~ ", paste(colnames(cbind(Vps)), collapse= "+")))
summary.match<-summary(to.check, method="reverse", overall=FALSE, test=TRUE, data=data.matched)
print(summary.match,long=TRUE,prtest=c('P'))

# Restrict to unique ids
cat("Restricting Summary Demographics to unique subjects","\n")
data.matched.unique<-data.temp[c(unique(rr$index.treated),unique(rr$index.control)),]
to.check<-as.formula(paste("treat ~ ", paste(colnames(cbind(Vps)), collapse= "+")))
summary.match.unique<-summary(to.check, method="reverse", overall=FALSE, test=TRUE, data=data.matched.unique)
print(summary.match.unique,long=TRUE,prtest=c('P'))
}



#out<-list(rr=rr,m1.M=m1.M,m0.M=m0.M,delta.M=delta.M,data.matched=data.matched,summary.match=summary.match)
out<-list(rr=rr,m1.M=m1.M,m0.M=m0.M,delta.M=delta.M,data.matched=data0.M1,data.matched.unique=data0.M1.unique)
return(out)
}




x.truncate<-function(x,truncate){
x.quant<-quantile(x,c(truncate,1-truncate))
x.trunc<-x
x.trunc[which(x<=x.quant[1])]<-x.quant[1]
x.trunc[which(x>=x.quant[2])]<-x.quant[2]
return(x.trunc)
}


#####################
# Weights by Factors
#####################

check.weight<-function(W,X,Delta){
par(mfrow=c(2,2))

hist(W[Delta==1 & X==1],xlab="inv-psE",ylab="Exp arm: X=1")
hist(W[Delta==0 & X==1],xlab="inv-psC",ylab="Control arm: X=1")
hist(W[Delta==1 & X==0],xlab="inv-psE",ylab="Exp arm: X=0")
hist(W[Delta==0 & X==0],xlab="inv-psC",ylab="Control arm: X=0")

cat("Sum of weights: Exp",c(sum(W[Delta==1])),"\n")
cat("Sum of weights: Exp, X=1,0",c(sum(W[Delta==1 & X==1]),sum(W[Delta==1 & X==0])),"\n")
cat("N: Exp, X=1,0",c(length(W[Delta==1 & X==1]),length(W[Delta==1 & X==0])),"\n")
cat("Sum of weights: Con",c(sum(W[Delta==0])),"\n")
cat("Sum of weights: Con, X=1,0",c(sum(W[Delta==0 & X==1]),sum(W[Delta==0 & X==0])),"\n")
cat("N: Con, X=1,0",c(length(W[Delta==0 & X==1]),length(W[Delta==0 & X==0])),"\n")
}

check.dhats.ByX<-function(D,X,Delta){
par(mfrow=c(2,2))
hist(D[Delta==1 & X==1],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Exp arm: X=1")
hist(D[Delta==0 & X==1],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Control arm: X=1")
hist(D[Delta==1 & X==0],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Exp arm: X=0")
hist(D[Delta==0 & X==0],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Control arm: X=1")
N<-length(D)
#cat("Sum of dhats: Exp",c(sum(D[Delta==1])),"\n")
#cat("Sum of dhats: Exp, X=1,0",c(sum(D[Delta==1 & X==1]),sum(D[Delta==1 & X==0])),"\n")
cat("N: Exp, X=1,0",c(length(D[Delta==1 & X==1]),length(D[Delta==1 & X==0])),"\n")
cat("Contribution to dhat: Exp, X=1,0",c(sum(D[Delta==1 & X==1])/N,sum(D[Delta==1 & X==0])/N),"\n")
cat("Contribution to dhat: Exp",c(sum(D[Delta==1])/N),"\n")
#cat("Sum of dhats: Con",c(sum(D[Delta==0])),"\n")
#cat("Sum of dhats: Con, X=1,0",c(sum(D[Delta==0 & X==1]),sum(D[Delta==0 & X==0])),"\n")
cat("N: Con, X=1,0",c(length(D[Delta==0 & X==1]),length(D[Delta==0 & X==0])),"\n")
cat("Contribution to dhat: Con, X=1,0",c(sum(D[Delta==0 & X==1])/N,sum(D[Delta==0 & X==0])/N),"\n")
cat("Contribution to dhat: Con",c(sum(D[Delta==0])/N),"\n")
cat("Overall dhat=",c(mean(D)),"\n")
}


check.dhats.ByTreat<-function(D,Delta){
par(mfrow=c(1,2))
hist(D[Delta==1],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Exp arm")
hist(D[Delta==0],xlab="dhat:P(Y=1;Treat=1)-P(Y=1; Treat=0)",ylab="Control arm")
N<-length(D)
cat("Contribution to dhat: Exp",c(sum(D[Delta==1])/N),"\n")
cat("Contribution to dhat: Con",c(sum(D[Delta==0])/N),"\n")
cat("Overall dhat=",c(mean(D)),"\n")
}




# Vps.exact --> covariates to match exactly
# note: predict.gam changed to predict
BR.TE.fit<-function(data,Vor.names,Vps.names,Vps.exact=NULL,treat.inter=FALSE,family.or=c("binomial"),details=FALSE,wgt.details=FALSE,use.stable=FALSE,
                    wt.truncate=0.0,get.main=FALSE,naive.only=FALSE,get.robust=FALSE){

require(Matching)

if(!any(names(data)=="Delta")) stop("Treatment group = Delta needs to be in dataset")

if(treat.inter==FALSE & any(Vor.names!="Delta")){
#cat("Linear treatment term in OR modeling (i.e., no interactions with treatment and covariates)")
Vor.names<-c("Delta",Vor.names) # Add treatment to OR covariates
}
# Otherwise, full model specification via Vor.names with interaction (by treatment)
OR.model<- as.formula(paste("Y ~ ", paste(Vor.names, collapse= "+")))
PS.model<- as.formula(paste("Delta ~ ", paste(Vps.names, collapse= "+")))

DR.names<-c(Vor.names,"inv.ps")
# g.DR is added to OR model in double-robust estimation
DR.model<- as.formula(paste("Y ~ ", paste(DR.names, collapse= "+")))
# DR version 2:  Estimates separate parameters for PS weights (g.DR1--> Exp; g.DR2--> Control)
DR2.names<-c(Vor.names,"inv.psE","inv.psC")
# g.DR1 and g.DR2 is added to OR model in double-robust estimation
DR2.model<- as.formula(paste("Y ~ ", paste(DR2.names, collapse= "+")))
# GAM model to allow for non-linear psE and psC
DR3.names<-c(Vor.names,"s(psE)","s(psC)")
# g.DR1 and g.DR2 is added to OR model in double-robust estimation
DR3.model<- as.formula(paste("Y ~ ", paste(DR3.names, collapse= "+")))
# Non-linearity in ps
DR4.names<-c(Vor.names,"s(ps)")
DR4.model<- as.formula(paste("Y ~ ", paste(DR4.names, collapse= "+")))

if(naive.only){
get.main<-FALSE # Turn off "main" estimators
get.all<-FALSE # Turn of "ALL" estimators
Y<-data$Y
Delta<-data$Delta
m.naive.1<-mean(Y[Delta==1])
m.naive.0<-mean(Y[Delta==0])
}

# get.main --> Restrict to recommended ("main") estimators
if(get.main){
get.all<-FALSE
# Naive and DR estimators
Y<-data$Y
Delta<-data$Delta
m.naive.1<-mean(Y[Delta==1])
m.naive.0<-mean(Y[Delta==0])

# DR
fit.ps<-glm(PS.model,family="binomial",data=data)
pihat.Vps<-fit.ps$fitted
wt.1<-1/pihat.Vps
wt.0<-1/(1-pihat.Vps)
g.DR<-ifelse(Delta==1,wt.1,wt.0)

data$inv.ps<-g.DR


if(!get.robust){
fit.dr<-glm(DR.model,family=family.or,data=data)
# Set Delta=1 and g.DR=wt.1
data.1<-data
data.1$Delta<-1
data.1$inv.ps<-wt.1
e.DVor1<-predict.glm(fit.dr,data.1,type="response")
# Set Delta=0 and g.DR=wt.0
data.0<-data
data.0$Delta<-0
data.0$inv.ps<-wt.0
e.DVor0<-predict.glm(fit.dr,data.0,type="response")
DR.1<-mean(e.DVor1)
DR.0<-mean(e.DVor0)
}
if(get.robust){
require(robustbase)
fit.dr<-glmrob(DR.model,family=family.or,data=data)
# Set Delta=1 and g.DR=wt.1
data.1<-data
data.1$Delta<-1
data.1$inv.ps<-wt.1
e.DVor1<-predict(fit.dr,data.1,type="response")
# Set Delta=0 and g.DR=wt.0
data.0<-data
data.0$Delta<-0
data.0$inv.ps<-wt.0
e.DVor0<-predict(fit.dr,data.0,type="response")
DR.1<-mean(e.DVor1)
DR.0<-mean(e.DVor0)
}


# Include OR estimator
fit.or<-glm(OR.model,family=family.or,data=data)
# Predictions set-to receive treatment for all
# Set Delta=1
data.1<-data
data.1$Delta<-1
pred.or1<-predict.glm(fit.or,data.1,type="response")
# Predictions set-to receive control for all
# Set Delta=0
data.0<-data
data.0$Delta<-0
pred.or0<-predict.glm(fit.or,data.0,type="response")
OR.1<-mean(pred.or1)
OR.0<-mean(pred.or0)

# Include 1:1 matching allowing for ties!

get.matching<-matching.fit(Vps.names=Vps.names,data=data,
ExactMatch.names=Vps.exact,fit.ps=fit.ps,details=details)

m.match.1<-get.matching$m1.M
m.match.0<-get.matching$m0.M
delta.match<-get.matching$delta.M

data.matched<-get.matching$data.matched

# Fit separate OR models for model checking

fit.or1<-glm(OR.model,family=family.or,data=data[which(data$Delta==1),]) # Fit for Delta=1
fit.or0<-glm(OR.model,family=family.or,data=data[which(data$Delta==0),]) # Fit for Delta=0


}

if(!get.main & !naive.only) get.all<-TRUE

# All estimators
if(get.all){
Y<-data$Y
Delta<-data$Delta
m.naive.1<-mean(Y[Delta==1])
m.naive.0<-mean(Y[Delta==0])

###################
#  OR model
###################

fit.or<-glm(OR.model,family=family.or,data=data)

#print(OR.model)
#print(fit.or)

# Predictions set-to receive treatment for all
# Set Delta=1
data.1<-data
data.1$Delta<-1
pred.or1<-predict.glm(fit.or,data.1,type="response")

# Predictions set-to receive control for all
# Set Delta=0

data.0<-data
data.0$Delta<-0
pred.or0<-predict.glm(fit.or,data.0,type="response")

OR.1<-mean(pred.or1)
OR.0<-mean(pred.or0)

# Fit PS model

fit.ps<-glm(PS.model,family="binomial",data=data)
pihat.Vps<-fit.ps$fitted

# HT estimator

HT.1<-mean(c(Delta*Y)/pihat.Vps)
HT.0<-mean(c((1-Delta)*Y)/(1-pihat.Vps))

wt.1<-1/pihat.Vps
wt.0<-1/(1-pihat.Vps)

# Raw PS
data$ps<-ifelse(Delta==1,pihat.Vps,1-pihat.Vps)

if(use.stable){
pihat.null<-glm(Delta~1,family="binomial")$fitted
# Stabilized
wt.1<-pihat.null/pihat.Vps
wt.0<-(1-pihat.null)/(1-pihat.Vps)
}

g.DR<-ifelse(Delta==1,wt.1,wt.0)

if(wt.truncate>0){
g.DR<-x.truncate(g.DR,wt.truncate)
}

if(wgt.details){
cat("Summary of weights=","\n")
print(summary(g.DR))
cat("Summary of wt.1 (Experimental arm)=","\n")
print(summary(g.DR[which(Delta==1)]))
cat("Summary of wt.0 (Control arm)=","\n")
print(summary(g.DR[which(Delta==0)]))
}

# Append to data

data$inv.ps<-g.DR

if(!get.robust){
fit.dr<-glm(DR.model,family=family.or,data=data)
# OR Delta=1 fits
# Set Delta=1 and g.DR=wt.1
data.1<-data
data.1$Delta<-1
data.1$inv.ps<-wt.1

e.DVor1<-predict.glm(fit.dr,data.1,type="response")

# Set Delta=0 and g.DR=wt.0
data.0<-data
data.0$Delta<-0
data.0$inv.ps<-wt.0

e.DVor0<-predict.glm(fit.dr,data.0,type="response")

DR.1<-mean(e.DVor1)
DR.0<-mean(e.DVor0)

# Output these for diagnostics
e.DVor1.v1<-e.DVor1
e.DVor0.v1<-e.DVor0
}

if(get.robust){
require(robustbase)
fit.dr<-glmrob(DR.model,family=family.or,data=data)
# OR Delta=1 fits
# Set Delta=1 and g.DR=wt.1
data.1<-data
data.1$Delta<-1
data.1$inv.ps<-wt.1

e.DVor1<-predict(fit.dr,data.1,type="response")

# Set Delta=0 and g.DR=wt.0
data.0<-data
data.0$Delta<-0
data.0$inv.ps<-wt.0

e.DVor0<-predict(fit.dr,data.0,type="response")

DR.1<-mean(e.DVor1)
DR.0<-mean(e.DVor0)

# Output these for diagnostics
e.DVor1.v1<-e.DVor1
e.DVor0.v1<-e.DVor0
}


# DR version 2
g.DR1<-Delta*wt.1
g.DR2<-(1-Delta)*wt.0

# For Gam modeling
data$psE<-Delta*pihat.Vps
data$psC<-(1-Delta)*(1-pihat.Vps)

data$inv.psE<-g.DR1
data$inv.psC<-g.DR2

fit.dr2<-glm(DR2.model,family=family.or,data=data)

# Set Delta=1 and g.DR1=wt.1 and g.DR2=0
data.1<-data
data.1$Delta<-1.0
data.1$inv.psE<-wt.1
data.1$inv.psC<-0.0

e.DVor1<-predict.glm(fit.dr2,data.1,type="response")

# Set Delta=0 and g.DR1=0 and g.DR2=wt.0
data.0<-data
data.0$Delta<-0.0
data.0$inv.psE<-0.0
data.0$inv.psC<-wt.0

e.DVor0<-predict.glm(fit.dr2,data.0,type="response")

e.DVor1.v2<-e.DVor1
e.DVor0.v2<-e.DVor0

DR.1.v2<-mean(e.DVor1)
DR.0.v2<-mean(e.DVor0)

# Third version (Funk et al., 2010)
# This is an augmented version (AIPW)


# Do not use stabilized weights here
ipw.1<-1/pihat.Vps
ipw.0<-1/(1-pihat.Vps)

term1.a<-c(Delta*Y*ipw.1)
term1.b<-c((Delta-pihat.Vps)*ipw.1)
term1.c<-c(pred.or1)
e.DVor1<-term1.a-term1.b*term1.c
DR.1.v3<-mean(e.DVor1)

term0.a<-c((1-Delta)*Y*ipw.0)
term0.b<-c((Delta-pihat.Vps)*ipw.0)
term0.c<-c(pred.or0)
e.DVor0<-term0.a+term0.b*term0.c
DR.0.v3<-mean(e.DVor0)

e.DVor1.v3<-e.DVor1
e.DVor0.v3<-e.DVor0


# Fourth version
# Estimate OR models separately within treatment groups

# Only for non-interaction OR models

DR.1.v4<-DR.0.v4<-NULL
e.DVor1.v4<-e.DVor0.v4<-NULL
fit.or1<-fit.or0<-NULL
OR.1.v2<-OR.0.v2<-NULL

if(!treat.inter){
Vor.names2<-Vor.names[-c(1)] # Remove Delta from model
OR.model<- as.formula(paste("Y ~ ", paste(Vor.names2, collapse= "+")))

fit.or1<-glm(OR.model,family=family.or,data=data[which(data$Delta==1),]) # Fit for Delta=1
fit.or0<-glm(OR.model,family=family.or,data=data[which(data$Delta==0),]) # Fit for Delta=0

pred.or1<-predict.glm(fit.or1,data,type="response")
pred.or0<-predict.glm(fit.or0,data,type="response")


# OR version 2
OR.1.v2<-mean(pred.or1)
OR.0.v2<-mean(pred.or0)

term1.a<-c(Delta*Y*ipw.1)
term1.b<-c((Delta-pihat.Vps)*ipw.1)
term1.c<-c(pred.or1)
e.DVor1<-term1.a-term1.b*term1.c
DR.1.v4<-mean(e.DVor1)

term0.a<-c((1-Delta)*Y*ipw.0)
term0.b<-c((Delta-pihat.Vps)*ipw.0)
term0.c<-c(pred.or0)
e.DVor0<-term0.a+term0.b*term0.c
DR.0.v4<-mean(e.DVor0)

e.DVor1.v4<-e.DVor1
e.DVor0.v4<-e.DVor0

}


# Fifth version
# Allow non-linearities in weights in g.DR

fit.dr3<-gam(DR3.model,family=family.or,data=data)

# Set Delta=1 and g.DR1=wt.1 and g.DR2=0
data.1<-data
data.1$Delta<-1.0
data.1$psE<-pihat.Vps
data.1$psC<-0.0

e.DVor1<-predict(fit.dr3,data.1,type="response")

# Set Delta=0 and g.DR1=0 and g.DR2=wt.0
data.0<-data
data.0$Delta<-0.0
data.0$psE<-0.0
data.0$psC<-(1-pihat.Vps)


e.DVor0<-predict(fit.dr3,data.0,type="response")

e.DVor1.v5<-e.DVor1
e.DVor0.v5<-e.DVor0

DR.1.v5<-mean(e.DVor1)
DR.0.v5<-mean(e.DVor0)


# Sixth version

fit.dr4<-gam(DR4.model,family=family.or,data=data)

# Set Delta=1 and g.DR1=wt.1 and g.DR2=0
data.1<-data
data.1$Delta<-1.0
data.1$ps<-pihat.Vps

e.DVor1<-predict(fit.dr4,data.1,type="response")

# Set Delta=0 and g.DR1=0 and g.DR2=wt.0
data.0<-data
data.0$Delta<-0.0
data.0$ps<-1-pihat.Vps

e.DVor0<-predict(fit.dr4,data.0,type="response")

e.DVor1.v6<-e.DVor1
e.DVor0.v6<-e.DVor0

DR.1.v6<-mean(e.DVor1)
DR.0.v6<-mean(e.DVor0)

get.matching<-matching.fit(Vps.names=Vps.names,data=data,
ExactMatch.names=Vps.exact,fit.ps=fit.ps,details=details)

m.match.1<-get.matching$m1.M
m.match.0<-get.matching$m0.M
delta.match<-get.matching$delta.M
data.matched<-get.matching$data.matched

return(list(m1.DR=DR.1,m0.DR=DR.0,
m1.match=m.match.1,m0.match=m.match.0,delta.match=delta.match,
data.matched=data.matched,
e.DVor1.v1=e.DVor1.v1,e.DVor0.v1=e.DVor0.v1,
e.DVor1.v2=e.DVor1.v2,e.DVor0.v2=e.DVor0.v2,
e.DVor1.v3=e.DVor1.v3,e.DVor0.v3=e.DVor0.v3,
e.DVor1.v4=e.DVor1.v4,e.DVor0.v4=e.DVor0.v4,
e.DVor1.v5=e.DVor1.v5,e.DVor0.v5=e.DVor0.v5,
e.DVor1.v6=e.DVor1.v6,e.DVor0.v6=e.DVor0.v6,
m1.OR=OR.1,m0.OR=OR.0,m1.HT=HT.1,m0.HT=HT.0,
m1.OR2=OR.1.v2,m0.OR2=OR.0.v2,
m1.naive=m.naive.1,m0.naive=m.naive.0,
m1.DR2=DR.1.v2,m0.DR2=DR.0.v2,
m1.DR3=DR.1.v3,m0.DR3=DR.0.v3,
m1.DR4=DR.1.v4,m0.DR4=DR.0.v4,
m1.DR5=DR.1.v5,m0.DR5=DR.0.v5,
m1.DR6=DR.1.v6,m0.DR6=DR.0.v6,
delta.naive=m.naive.1-m.naive.0,
delta.DR=DR.1-DR.0,delta.OR=OR.1-OR.0,delta.HT=HT.1-HT.0,
delta.OR2=OR.1.v2-OR.0.v2,
delta.DR2=DR.1.v2-DR.0.v2,
delta.DR3=DR.1.v3-DR.0.v3,
delta.DR4=DR.1.v4-DR.0.v4,
delta.DR5=DR.1.v5-DR.0.v5,
delta.DR6=DR.1.v6-DR.0.v6,
wt=g.DR,
fit.or=fit.or,fit.ps=fit.ps,fit.dr=fit.dr,fit.dr2=fit.dr2,fit.dr3=fit.dr3,fit.dr4=fit.dr4,fit.or1=fit.or1,fit.or0=fit.or0))
}


if(get.main){
return(list(m1.naive=m.naive.1,m0.naive=m.naive.0,delta.naive=m.naive.1-m.naive.0,
m1.DR=DR.1,m0.DR=DR.0,delta.DR=DR.1-DR.0,
m1.OR=OR.1,m0.OR=OR.0,delta.OR=OR.1-OR.0,
fit.or=fit.or,fit.ps=fit.ps,fit.dr=fit.dr,fit.or1=fit.or1,fit.or0=fit.or0,
m1.match=m.match.1,m0.match=m.match.0,delta.match=delta.match,data.matched=data.matched))
}

if(naive.only){
return(list(m1.naive=m.naive.1,m0.naive=m.naive.0,delta.naive=m.naive.1-m.naive.0))
}
}


BR.TE.CI<-function(Vor.names,Vps.names,Vps.exact=NULL,family.or=c("binomial"),details=FALSE,draws=0,alpha=0.10,treat.inter=FALSE,use.stable=FALSE,wt.truncate=0.0,data,est.delta.names,est.m1.names,est.m0.names,
get.main=FALSE,naive.only=FALSE,seed.start=8316951,get.robust=FALSE,wgt.details=FALSE){

# naive.only --> only return un-adjusted estimates
# get.main --> only return "main" estimators tbd pending simulation results
# turn off get.main if naive.only
if(naive.only) get.main<-FALSE

calpha<-qnorm(1-alpha)                                                                                    

BR.estimates<-BR.TE.fit(Vor.names=Vor.names,Vps.names=Vps.names,Vps.exact=Vps.exact,get.robust=get.robust,wgt.details=wgt.details,
family.or=family.or,details=details,treat.inter=treat.inter,wt.truncate=wt.truncate,use.stable=use.stable,data=data,get.main=get.main,naive.only=naive.only)

if(is.null(est.delta.names)) est.delta.names<-c("delta.naive","delta.HT","delta.OR","delta.OR2","delta.DR","delta.DR2","delta.DR3","delta.DR4","delta.DR5","delta.DR6","delta.match")
if(is.null(est.m1.names)) est.m1.names<-c("m1.naive","m1.HT","m1.OR","m1.OR2","m1.DR","m1.DR2","m1.DR3","m1.DR4","m1.DR5","m1.DR6","m1.match")
if(is.null(est.m0.names)) est.m0.names<-c("m0.naive","m0.HT","m0.OR","m0.OR2","m0.DR","m0.DR2","m0.DR3","m0.DR4","m0.DR5","m0.DR6","m0.match")


BR.delta.obs<-as.numeric(BR.estimates[c(est.delta.names)])
BR.m1.obs<-as.numeric(BR.estimates[c(est.m1.names)])
BR.m0.obs<-as.numeric(BR.estimates[c(est.m0.names)])

names(BR.delta.obs)<-est.delta.names
names(BR.m1.obs)<-est.m1.names
names(BR.m0.obs)<-est.m0.names

if(draws==0) return(list(estimates=BR.estimates,BR.delta.obs=BR.delta.obs,BR.m1.obs=BR.m1.obs,BR.m0.obs=BR.m0.obs))

if(draws>0){
t.start<-proc.time()[1]

# Bootstrap difference estimators
# 1 naive version
# 1 HT version
# 2 OR versions
# 6 DR versions

# Differences
BR.delta.boots<-matrix(NA,nrow=draws,ncol=length(est.delta.names))
# Treatment
BR.m1.boots<-matrix(NA,nrow=draws,ncol=length(est.m1.names))
# Control
BR.m0.boots<-matrix(NA,nrow=draws,ncol=length(est.m0.names))

n<-length(data$Y)
data.id<-c(1:n)
data.obs<-data
for(dd in 1:draws){
set.seed(seed.start+dd*1000)
id.boot<-sample(data.id,size=n,replace=T)
data.boot<-data[id.boot,]
BR.fit.boots<-try(BR.TE.fit(Vor.names=Vor.names,Vps.names=Vps.names,Vps.exact=Vps.exact,
family.or=family.or,details=FALSE,treat.inter=treat.inter,wt.truncate=wt.truncate,
use.stable=use.stable,data=data.boot,get.main=get.main,naive.only=naive.only),TRUE)

if(!inherits(BR.fit.boots, "try-error")){

d.boot<-as.numeric(BR.fit.boots[c(est.delta.names)])
m1.boot<-as.numeric(BR.fit.boots[c(est.m1.names)])
m0.boot<-as.numeric(BR.fit.boots[c(est.m0.names)])

BR.delta.boots[dd,]<-c(d.boot)
BR.m1.boots[dd,]<-c(m1.boot)
BR.m0.boots[dd,]<-c(m0.boot)
}

if(details){
if(dd/draws==0.1) cat("10% of Bootstraps complete","\n")
if(dd/draws==0.2) cat("20% complete","\n")
if(dd/draws==0.3) cat("30% complete","\n")
if(dd/draws==0.4) cat("40% complete","\n")
if(dd/draws==0.5) cat("50% complete","\n")
if(dd/draws==0.6) cat("60% complete","\n")
if(dd/draws==0.7) cat("70% complete","\n")
if(dd/draws==0.8) cat("80% complete","\n")
if(dd/draws==0.9) cat("90% complete","\n")
}

}


colnames(BR.delta.boots)<-est.delta.names
colnames(BR.m1.boots)<-est.m1.names
colnames(BR.m0.boots)<-est.m0.names

SE.boots<-sqrt(apply(BR.delta.boots,2,var,na.rm=TRUE))
SE.m1.boots<-sqrt(apply(BR.m0.boots,2,var,na.rm=TRUE))
SE.m0.boots<-sqrt(apply(BR.m1.boots,2,var,na.rm=TRUE))

LB.delta<-BR.delta.obs-calpha*SE.boots
UB.delta<-BR.delta.obs+calpha*SE.boots

LB.m1<-BR.m1.obs-calpha*SE.m1.boots
UB.m1<-BR.m1.obs+calpha*SE.m1.boots

LB.m0<-BR.m0.obs-calpha*SE.m0.boots
UB.m0<-BR.m0.obs+calpha*SE.m0.boots

names(LB.delta)<-names(UB.delta)<-est.delta.names
names(LB.m1)<-names(UB.m1)<-est.m1.names
names(LB.m0)<-names(UB.m0)<-est.m0.names

if(details){
cat("# of bootstrap resamples=",c(sum(!is.na(BR.delta.boots[,c("delta.DR")]))),"\n")
t.end<-proc.time()[1]
t.min<-(t.end-t.start)/60
cat("# Time in minutes",c(t.min),"\n")
}

out<-list(estimates=BR.estimates,
BR.delta.obs=BR.delta.obs,
BR.m1.obs=BR.m1.obs,
BR.m0.obs=BR.m0.obs,
BR.delta.boots=BR.delta.boots,
SE.delta=SE.boots,
BR.m1.boots=BR.m1.boots,
BR.m0.boots=BR.m0.boots,
LB.delta=LB.delta,UB.delta=UB.delta,
LB.m1=LB.m1,UB.m1=UB.m1,
LB.m0=LB.m0,UB.m0=UB.m0,
data.matched=BR.estimates$data.matched,
fit.or=BR.estimates$fit.or,fit.ps=BR.estimates$fit.ps,
fit.dr=BR.estimates$fit.dr,fit.dr2=BR.estimates$fit.dr2,
fit.dr3=BR.estimates$fit.dr3)
return(out)
}
}
