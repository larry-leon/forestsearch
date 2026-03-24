
########################################################################################################################
# BR.boots and alpha.BR --> For boostrap estimation of DR alpha-level (alpha=alpha.BR) CI's;  default=turn off
# print.level=0 --> No printing; 2= most printing
# n0.max set to # of subject in historical control dataset ... assuming HC is fixed!
# If HC is evolving, then set n0.max to planned total # of HC subjects
# get.matching --> False then NO matching estimators are implemented
# M = number of matches; Default M=1 for 1:1 matching
# Matching.type=1 --> Only PS matching from "Match" package is implemented
# Matching.type=2 ---> PS matching ("Match" package) where exact matching is implemented for factor=Vps[,exact.index]
# For example, if index=c(1), then the 1st factor in Vps is matched exactly on;  Followed by PS matching
# Matching.type=3 ---> PS matching based on nearest neighborhood as implemented in "MatchIt" package
# Matching.type=0 ---> Each of 1-3 above is implemented (i.e., get all 3!)
##########################################################################################################################

ratio.zero<-function(x,y) { 
res <- x / y 
res[ is.na(res) ] <- 0 
return(res) } 


# Note:  This is exploratory --- In progress!!!!
check.post<-function(n,Y,m.boots,xlabel="'Adjusted P'",check.weighted=FALSE,parG=c(1,1)){
# "Estimator" bootstrap density
m.boots<-na.exclude(m.boots)
est.boots<-density(m.boots)
ee<-est.boots$x
a1<-parG[1]+Y
b1<-parG[2]+n-Y
phat<-Y/n
# Parametric "quasi" bootstrap
#m.quasi.boots<-rbinom(n=length(m.boots),size=n,prob=phat)/n
Y.boots<-rbinom(n=20000,size=n,prob=phat)
mY.boots<-sort(Y.boots/n)
# m.quasi.boots --> Sample Space for P's!
# Density of Y's corresponding to P's (Y.boots)
BB.Y<-sort(Y.boots)  # --> "Betas"=BB.Y/n
# Density of "Y/n" evaluated at observed (n,phat)
g.bhat.beta<-n*dbinom(x=BB.Y,size=n,prob=phat) # Bootstrap density!!!
# Represents the density from which the bootstrap samples are drawn
fit.bootY<-hist(mY.boots,probability=TRUE,xlab=xlabel,main="'P-Densities'")
# Posterior 
#density.post<-dbeta(x=xx,shape1=a1,shape2=b1)
density.post<-dbeta(x=mY.boots,shape1=a1,shape2=b1)
lines(mY.boots,density.post,type="l",lty=1,col="grey",lwd=4)
lines(ee,est.boots$y,col = "blue",lwd=0.5)  # Kernel estimate of "DR estimator" bootstrap density
lines(BB.Y/n,g.bhat.beta,col="black",lty=2,lwd=3)
rug(phat,col="black",lwd=4)
#rug(mean(m.boots),col="blue",lwd=4)

ymax<-max(density.post)

# Next, need the "likelihood": The density of bhat
# --> Take this from non-parametric "DR estimator" bootstrap
# We want to take ratio of g.beta.bhat and g.bhat.beta
# --> Need to evaluate at same points!
if(check.weighted){
g.beta.bhat<-est.boots$y # Evaluated at "ee" = "DR estimator" bootstrap points
# Okay, g.beta.bhat is a nonparametric density estimate of DR distribution
# Doesn't seem to work in this context
# So, let's try using the density from asymptotic normal based on SD estimated 
# from bootstraps
mu.like.asy<-mean(m.boots)
sd.like.asy<-sqrt(var(m.boots))
g.beta.bhat<-dnorm(x=ee,mean=phat,sd=sd.like.asy)
# Need to evaluate g.bhat.beta @ ee
g.bhat.beta.ee<-n*dbinom(x=round(n*ee),size=n,prob=phat)

plot(ee,g.beta.bhat,xlab=xlabel,main="'P-Densities'",lty=1,col="black",lwd=2)
lines(ee,g.bhat.beta.ee,lty=2,col="red",lwd=3)

R.beta<-ratio.zero(g.beta.bhat,g.bhat.beta.ee)
W.beta<-R.beta*dunif(ee)
cat("g.beta.bhat","\n")
print(summary(g.beta.bhat))
cat("g.bhat.beta","\n")
print(summary(g.bhat.beta))
cat("R.beta","\n")
print(summary(R.beta))
cat("W.beta","\n")
print(summary(W.beta))

lines(ee,W.beta,col="red",lwd=4,lty=5)
#cat("Length g.beta.bhat and g.bhat.beta.ee=",c(length(g.beta.bhat),length(g.bhat.beta.ee)),"\n")
cat("Sum pi*R=",sum(W.beta,na.rm=TRUE),"\n")
#norm.const<-sum(W.beta*diff(ee)) 
dx<-diff(ee)
intx<-W.beta
norm.const<-sum(dx*intx[-1],na.rm=TRUE)
cat("Sum pi(x)*R(x)dx=",c(norm.const),"\n")
#print(summary(g.beta.bhat))
#print(summary(g.bhat.beta.ee))
#print(summary(W.beta))
#print(summary(R.beta))
#W.beta<-W.beta/sum(W.beta,na.rm=TRUE)
W.beta<-W.beta/norm.const
lines(ee,W.beta,col="grey",lwd=4,lty=5)
#like.quasi<-n*dbinom(x=round(xx.quasi*n), size=n, prob=phat)
# Normal approximation for likelihood
#like.quasi<-dnorm(x=xx.quasi,mean=phat,sd=sqrt(phat*(1-phat)/n))
#boot.quasi<-fit.boot.quasi$y
#R.quasi<-ratio.zero(like.quasi,boot.quasi)
#W.density<-R.quasi*dunif(xx.quasi)
#lines(xx.quasi,W.density,col="red",lwd=3)
#W.density<-W.density/sum(W.density)
#lines(xx.quasi,W.density,col="grey",lwd=3)
}
return(list(ymax=ymax))
}


Est.stats<-function(est.fits,m1.name,m0.name,delta.name,n1,n0,n1.max,n0.max,parE,parS,print.level=1,delta.margin,theta.T,BR.alpha,BR.boots,check.sampling=FALSE){
m1.obs<-as.numeric(est.fits$BR.m1.obs[c(m1.name)])
m0.obs<-as.numeric(est.fits$BR.m0.obs[c(m0.name)])
delta.obs<-as.numeric(est.fits$BR.delta.obs[c(delta.name)])

if(BR.boots>0){
LB.delta<-as.numeric(est.fits$LB.delta[c(delta.name)])
UB.delta<-as.numeric(est.fits$UB.delta[c(delta.name)])
}
Y1.est<-round(n1*m1.obs)
Y0.est<-round(n0*m0.obs)
if(check.sampling){
if(BR.boots==0) stop("Need BR.boots>0")
#win.graph()
m1.boots<-as.numeric(est.fits$BR.m1.boots[,c(m1.name)])
m0.boots<-as.numeric(est.fits$BR.m0.boots[,c(m0.name)])
temp<-check.post(n=n1,Y=Y1.est,m.boots=m1.boots,parG=parE,xlabel=m1.name)
legend(min(m1.boots),0.75*temp$ymax,c("Parametric Boot","NonPar Boot","Posterior"),col=c("black","blue","grey"),lty=c(2,1,1),lwd=c(3,0.5,4),bty="n")
rm("temp")
temp<-check.post(n=n0,Y=Y0.est,m.boots=m0.boots,parG=parS,xlabel=m0.name)
legend(min(m0.boots),0.5*temp$ymax,c("Parametric Boot","NonPar Boot","Posterior"),col=c("black","blue","grey"),lty=c(2,1,1),lwd=c(3,0.5,4),bty="n")
rm("temp")
}
if(print.level>=3) cat("Experimental and Historical-Control sample sizes (n1,n0)=",c(n1,n0),"\n")

PredProb<-as.numeric(predprobDist(x=Y1.est, n=n1, Nmax=n1.max, xS=Y0.est, nS=n0, NmaxControl=n0.max,
delta=delta.margin, thetaT=theta.T, parE=parE,parS=parS))

PostProb<-as.numeric(postprobDist(x=Y1.est, n=n1, xS=Y0.est, nS=n0, delta=delta.margin, parE=parE,parS=parS))

if(print.level>=3){
cat("##########################################################################################","\n")
cat("Estimator=",c(delta.name),"\n")
cat("m1,m0,diff=",c(m1.obs,m0.obs,delta.obs),"\n")
if(BR.boots>0){
cat("Bootstrap resamples=",c(BR.boots),"\n")
cat("2-sided CI level %=",c(100*(1-2*BR.alpha)),"\n")

cat("Lower,Upper=",c(LB.delta,UB.delta),"\n")
cat("Note: Lower bound can be compared to 0 for 1-sided superiority comparison at alpha=",c(BR.alpha),"\n")
}
cat("Predicted Probability (PP)=",c(PredProb),"\n")
cat("PP evaluation of Delta>delta.margin comparison, delta.margin=",c(delta.margin),"\n")

cat("Posterior Probability=",c(PostProb),"\n")
cat("Posterior Probability of Delta>delta.margin comparison, delta.margin=",c(delta.margin),"\n")


}
return(list(PP=PredProb,Post=PostProb))
}



check.est.post<-function(est.fits,m1.name,m0.name,n1,n0,n1.max,n0.max,parE,parS){
m1.obs<-as.numeric(est.fits$BR.m1.obs[c(m1.name)])
m0.obs<-as.numeric(est.fits$BR.m0.obs[c(m0.name)])
Y1.est<-round(n1*m1.obs)
Y0.est<-round(n0*m0.obs)
m1.boots<-as.numeric(est.fits$BR.m1.boots[,c(m1.name)])
m0.boots<-as.numeric(est.fits$BR.m0.boots[,c(m0.name)])
temp<-check.post(n=n1,Y=Y1.est,m.boots=m1.boots,parG=parE,xlabel=m1.name)
legend(min(m1.boots),0.75*temp$ymax,c("Parametric Boot","NonPar Boot","Posterior"),col=c("black","blue","grey"),lty=c(2,1,1),lwd=c(3,0.5,4),bty="n")
rm("temp")
temp<-check.post(n=n0,Y=Y0.est,m.boots=m0.boots,parG=parS,xlabel=m0.name)
legend(min(m0.boots),0.75*temp$ymax,c("Parametric Boot","NonPar Boot","Posterior"),col=c("black","blue","grey"),lty=c(2,1,1),lwd=c(3,0.5,4),bty="n")
rm("temp")
}





CTC.adjusted<-function(Vor.names,Vps.names,Vps.exact=NULL,data,BR.boots=0,BR.alpha=0.2,wt.truncate=0.0,wgt.details=FALSE,
print.level=0,n1.max,n0.max=sum(data$Delta==0),delta.margin,theta.T=0.60,
get.main=FALSE,
treat.inter=FALSE,use.stable=FALSE,get.robust=FALSE,
est.delta.names=NULL,est.m1.names=NULL,est.m0.names=NULL,check.PPsampling=FALSE,parE=c(1,1),parS=c(1,1)){

if(n1.max<sum(data$Delta==1))  stop("n1.max cannot be smaller than the current sample size: Can be equal to current sample size, or the anticipated 'final' sample size")
  
t.start<-proc.time()[1]

# est.delta.names specifies names of estimators to calculate
# If null, then get all
if(!get.main){
if(is.null(est.delta.names)){
est.delta.names<-c("delta.naive","delta.HT","delta.OR","delta.OR2","delta.DR","delta.DR2","delta.DR3","delta.DR4","delta.DR5","delta.DR6","delta.match")
est.m1.names<-c("m1.naive","m1.HT","m1.OR","m1.OR2","m1.DR","m1.DR2","m1.DR3","m1.DR4","m1.DR5","m1.DR6","m1.match")
est.m0.names<-c("m0.naive","m0.HT","m0.OR","m0.OR2","m0.DR","m0.DR2","m0.DR3","m0.DR4","m0.DR5","m0.DR6","m0.match")
}
}

if(get.main){
  est.delta.names<-c("delta.naive","delta.OR","delta.DR","delta.match")
  est.m1.names<-c("m1.naive","m1.OR","m1.DR","m1.match")
  est.m0.names<-c("m0.naive","m0.OR","m0.DR","m0.match")
}  
  

# Baseline demographics!
if(print.level>1){
data.temp<-as.data.frame(data)
data.temp$treat<-factor(data.temp$Delta,labels=c("Historical-Control","Experimental"))
label(data.temp$treat)<-c("Study")
to.check<-as.formula(paste("treat ~ ", paste(c(Vps.names,"Y"), collapse= "+")))
summary.raw<-summary(to.check, method="reverse", overall=FALSE, test=TRUE, data=data.temp)
print(summary.raw,long=TRUE,prtest=c('P'))
}


PP.BR<-rep(NA,length(est.delta.names))
Post.BR<-rep(NA,length(est.delta.names))


br.fits<-BR.TE.CI(Vor.names=Vor.names,Vps.names=Vps.names,Vps.exact=Vps.exact,get.robust=get.robust,wgt.details=wgt.details,
                  get.main=get.main,
data=data,family.or=c("binomial"),details=(print.level>1),draws=BR.boots,alpha=BR.alpha,wt.truncate=wt.truncate,use.stable=use.stable,
treat.inter=treat.inter,est.delta.names=est.delta.names,est.m1.names=est.m1.names,est.m0.names=est.m0.names,naive.only=FALSE)

n1<-sum(data$Delta)
n0<-sum(1-data$Delta)

for(nm in 1:length(est.delta.names)){
get.est<-Est.stats(est.fits=br.fits,
m1.name=est.m1.names[nm],m0.name=est.m0.names[nm],delta.name=est.delta.names[nm],
n1=n1,n0=n0,n1.max=n1.max,n0.max=n0.max,delta.margin=delta.margin,theta.T=theta.T,
BR.boots=BR.boots,BR.alpha=BR.alpha,parE=parE,parS=parS,print.level=print.level,check.sampling=check.PPsampling)
PP.BR[nm]<-get.est$PP
Post.BR[nm]<-get.est$Post
}

names(PP.BR)<-est.delta.names

if(BR.boots>0){
Delta.BR<-c(as.numeric(br.fits$BR.delta.obs[c(est.delta.names)]))
names(Delta.BR)<-c(est.delta.names)
result<-list(PP.BR=PP.BR,Delta.BR=Delta.BR,br.fits=br.fits)
dhat.all<-c(br.fits$BR.delta.obs)
m1.all<-c(br.fits$BR.m1.obs)
m0.all<-c(br.fits$BR.m0.obs)
dhat.lower<-c(br.fits$LB.delta)
dhat.upper<-c(br.fits$UB.delta)
SE.dhat<-c(br.fits$SE.delta)
PP.all<-c(PP.BR)
Post.all<-c(Post.BR)

if(!get.main){
out.all<-as.data.frame(rbind(dhat.all,m1.all,m0.all,dhat.lower,dhat.upper,SE.dhat,PP.all,Post.all))
rownames(out.all)<-c("Delta","m1","m0","LB(Delta)","UB(Delta)","SE(Delta)","PredProb","PostProb")
colnames(out.all)<-c("Naive","HT","OR","OR2","DR","DR2","DR3","DR4","DR5","DR6","1:all")
}

if(get.main){
  out.all<-as.data.frame(rbind(dhat.all,m1.all,m0.all,dhat.lower,dhat.upper,SE.dhat,PP.all,Post.all))
  rownames(out.all)<-c("Delta","m1","m0","LB(Delta)","UB(Delta)","SE(Delta)","PredProb","PostProb")
  colnames(out.all)<-c("Naive","OR","DR","1:all")
}
 }


if(BR.boots==0){
Delta.BR<-c(as.numeric(br.fits$BR.delta.obs[c(est.delta.names)]))
names(Delta.BR)<-c(est.delta.names)
result<-list(PP.BR=PP.BR,Delta.BR=Delta.BR,br.fits=br.fits)
dhat.all<-c(br.fits$BR.delta.obs)
m1.all<-c(br.fits$BR.m1.obs)
m0.all<-c(br.fits$BR.m0.obs)
dhat.lower<-dhat.upper<-rep(NA,length(dhat.all))
SE.dhat<-c(br.fits$SE.delta)
PP.all<-c(PP.BR)
Post.all<-c(Post.BR)
out.all<-as.data.frame(rbind(dhat.all,m1.all,m0.all,dhat.lower,dhat.upper,SE.dhat,PP.all,Post.all))


if(!get.main){
  out.all<-as.data.frame(rbind(dhat.all,m1.all,m0.all,dhat.lower,dhat.upper,SE.dhat,PP.all,Post.all))
  rownames(out.all)<-c("Delta","m1","m0","LB(Delta)","UB(Delta)","SE(Delta)","PredProb","PostProb")
  colnames(out.all)<-c("Naive","HT","OR","OR2","DR","DR2","DR3","DR4","DR5","DR6","1:all")
}

if(get.main){
  out.all<-as.data.frame(rbind(dhat.all,m1.all,m0.all,dhat.lower,dhat.upper,SE.dhat,PP.all,Post.all))
  rownames(out.all)<-c("Delta","m1","m0","LB(Delta)","UB(Delta)","SE(Delta)","PredProb","PostProb")
  colnames(out.all)<-c("Naive","OR","DR","1:all")
}


}

t.end<-proc.time()[1]
t.min<-(t.end-t.start)/60
if(print.level>=1) cat("# Time in minutes",c(t.min),"\n")

if(print.level>=3){
print(out.all)
}

if(check.PPsampling & BR.boots>0){
# check sampling for naive & DR
par(mfrow=c(2,2))
# naive estimator
m1.name<-"m1.naive"
m0.name<-"m0.naive"
check.est.post(est.fits=br.fits,m1.name=m1.name,m0.name=m0.name,n1=n1,n0=n0,n1.max=n1.max,n0.max=n0.max,parE=parE,parS=parS)
m1.name<-"m1.DR"
m0.name<-"m0.DR"
check.est.post(est.fits=br.fits,m1.name=m1.name,m0.name=m0.name,n1=n1,n0=n0,n1.max=n1.max,n0.max=n0.max,parE=parE,parS=parS)
}

return(list(result=result,est.summary=out.all,br.fits=br.fits))
}
