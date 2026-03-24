############################################################################
# This function fits the standard Bang & Robins Double-Robust estimator
############################################################################
mu.DoubleRobust<-function(Y2,Z,c,R,V.q,V.h,V.pi,checking,cov.OR.toplot=1,cov.PS.toplot=1,
inter.ORPS=TRUE,family,draws,se.method="Bootstrap",abstol=0.005,
trim.HT=0.05,trim.DR=0.05,trim.OR=0.05){

# Function only valid for gaussian and logistic families

if(!inter.ORPS & draws>0) stop("Not correct for no intercept OR-PS Model")

family.valid<-is.element(family,c("gaussian","binomial"))
if(!family.valid) stop("Only gaussian or binomial families are allowed")

n<-length(Y2)
id<-(Z==c)
nc<-sum(id); Rc<-R[id]; V.qc<-as.matrix(V.q[id,]); V.hc<-as.matrix(V.h[id,]); V.pic<-as.matrix(V.pi[id,])
Y2c<-Y2[id]; Zc=Z[id]


if(checking) cat("Number of subjects in group",nc,"\n")

case<-(Rc==1)  # Completers

################### Outcome-Response (OR) Model -- Bang-Robins  ################
# e_q Model of Davidian, Tsiatis, Leon (DTL)
# which is the e_OR model of Bang & Robins (BR)
################################################################################

b.OR<-as.vector(glm(Y2c[case]~as.matrix(V.qc[case,]), family=family)$coefficients)

z<-cbind(rep(1,nc),V.qc) 
xbeta<-z%*%b.OR

if(family=="binomial") e.OR<-pi.inf(exp(xbeta)/(1+exp(xbeta))) 
if(family=="gaussian") e.OR<-xbeta

ehat.q<-e.OR  # DTL e_q - model

zall<-cbind(rep(1,n),V.q)  # Obtain fitted values for ALL subjects (i.e., in BOTH groups)
xbeta<-zall%*%b.OR
if(family=="binomial") ehat.q.all<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") ehat.q.all<-xbeta

# B-R OR model estimator

mu.BR.OR<-mean(e.OR)

#########################################################
#          Fit e_h -- only required for DTL             #
#########################################################

b.eh<-as.vector(glm(Y2c[case]~as.matrix(V.hc[case,]), family=family)$coefficients)
z<-cbind(rep(1,nc),V.hc) 
xbeta<-z%*%b.eh

if(family=="binomial") ehat.h<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") ehat.h<-xbeta

zall<-cbind(rep(1,n),V.h)  # Obtain fitted values for all subjects
xbeta<-zall%*%b.eh
if(family=="binomial") ehat.h.all<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") ehat.h.all<-xbeta

#################################################################
#                Fit pi -- missing data model                   #
#################################################################

fit.pi<-glm(Rc~as.matrix(V.pic), family=binomial)

if(checking) print(summary(fit.pi))

pihat<-fitted.values(fit.pi)

if(checking){
print(summary(pihat))
print(summary(1/pihat))
}

b.pi<-as.vector(fit.pi$coefficients)

#################################################################################
############       Plot Drop-Out and Outcome Response Model Fits ################
if(checking){
plot(V.pic[,cov.PS.toplot],Rc,xlab="Covariate",ylab="Drop-Out Model")
xx<-V.pic[,cov.PS.toplot]; id.sort<-order(xx); xx<-xx[id.sort]
lines(xx,pihat[id.sort],type="l",lty=1,col="blue")
lines(lowess(xx,Rc[id.sort]),type="l",lty=2,col="red")

plot(V.qc[case,cov.OR.toplot],Y2c[case],xlab="Cov1",ylab="Outcome Response Model")
xx<-V.qc[case,cov.OR.toplot]; id.sort<-order(xx); xx<-xx[id.sort]
yy<-Y2c[case]
lines(xx,ehat.q[id.sort],type="l",lty=1,col="blue")
lines(lowess(xx,yy[id.sort]),type="l",lty=2,col="red")
}
###################################################################################
################ Next, implement the Bang-Robins DR Estimator #####################
#          Propensity-Score Model is already fitted -- pihat0
###################################################################################

Y.HT<-Rc*Y2c/pihat

if(checking){
print(summary(Rc))
print(summary(Y2c))
}

mu.BR.HT<-mean(Y.HT)  # Horvitz-Thompson Estimator

# Outcome-Response / PS Model (ORPS) 
V.ORPS<-cbind(V.qc,1/pihat) # Adding propensity-score to OR model

if(!inter.ORPS) b.ORPS<-as.vector(glm(Y2c[case]~as.matrix(V.ORPS[case,])-1, family=family)$coefficients)
if(inter.ORPS) b.ORPS<-as.vector(glm(Y2c[case]~as.matrix(V.ORPS[case,]), family=family)$coefficients)

if(inter.ORPS) z<-as.matrix(cbind(rep(1,nc),V.ORPS))
if(!inter.ORPS) z<-V.ORPS

xbeta<-z%*%b.ORPS

if(family=="binomial") e.ORPS<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") e.ORPS<-xbeta

mu.BR.DR<-mean(e.ORPS)              # Bang-Robins Double-Robust Estimator

########################################
# Contribution from non-cases
########################################
mu.BR.DR.noncase<-sum(e.ORPS[!case])/nc
cat("Contribution from non-cases: n,n_nc,mean",c(nc,sum(!case),mu.BR.DR.noncase),"\n")

cat("Estimators (DR,OR,HT)",c(mu.BR.DR,mu.BR.OR,mu.BR.HT),"\n")

if(draws>0) Gmat<-matrix(rnorm(nc*draws,0,1),nrow=nc,ncol=draws)

mu.BR.HT.star<-rep(NA,draws)
mu.BR.DR.star<-rep(NA,draws)
mu.BR.OR.star<-rep(NA,draws)

mu.BR.DR.noncase.star<-rep(NA,draws)

################################################################
# Observed covariates -- Fix the observed covariates for 
# re-sampling                                           
################################################################

xx.PS.obs<-as.matrix(cbind(rep(1,nc),V.pic))
xx.ORPS.obs<-as.matrix(cbind(rep(1,nc),V.ORPS))  
e.ORPS.obs<-e.ORPS

xx.OR.obs<-as.matrix(cbind(rep(1,nc),V.qc))

var.BR.DR<-1
var.BR.HT<-1
var.BR.OR<-1

if(draws>0){

if(se.method=="Bootstrap"){ 
data.id<-c(1:nc)
data.obs<-list(Rc=Rc,Y2c=Y2c,V.pic=V.pic,V.qc=V.qc)
}

for(dd in 1:draws){

if(se.method=="Parzen"){

# First, resample PS estimator

gg<-Gmat[,dd]
xbeta<-xx.PS.obs%*%b.pi
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector(Rc-m.b)
temp<-xx.PS.obs*resid*gg                       # (n x p) matrix
q.star<-apply(temp,2,mean)

if(checking) cat("PS Model Parzen Resampling Draw=",c(dd), "\n")

b.start<-b.pi
fit.NR<-QLE.NR(b.start=b.pi,Y=Rc,X=as.matrix(xx.PS.obs),q.star=q.star,abstol=abstol,checking=checking,family=family)

if(!is.null(fit.NR)){
bstar.pi<-c(fit.NR$beta)
Qsoln<-fit.NR$norm
NR.iter<-fit.NR$iter
converged.PS<-TRUE
if(checking) cat("Good: PS Converged via N-R (Qsoln,iter)",c(Qsoln,NR.iter),"\n")
}

# If Newton-Raphson Does not converge use optimization routines via "rescue" 

if(is.null(fit.NR)){ 
b.start<-b.pi
temp<-rescue.resample.QLE(b.start=b.start,Y=Rc,X=xx.PS.obs,q.star=q.star,
family=family,abstol=abstol,checking=checking,draw=dd)
bstar.pi<-c(temp$bstar)
converged.PS<-temp$converged
}

# If still does not converge -- try again with b=0 as initial value
# and less stringent convergence criteria

if(!converged.PS){
#b.zero<-rep(0,length(b.pi))
b.start<-temp$b.combo
abstol2<-2*abstol # Decrease the convergence hurdle
if(checking) cat("Trying second iteration of rescue for PS Model (abstol=)",c(abstol2),"\n")
temp<-rescue.resample.QLE(b.start=b.start,Y=Rc,X=xx.PS.obs,q.star=q.star,
family=family, abstol=abstol2,checking=checking,draw=dd)
bstar.pi<-c(temp$bstar)
converged.PS<-temp$converged
}


# If PS model converged
if(converged.PS){
xbeta.star<-xx.PS.obs%*%bstar.pi
if(family=="binomial") pistar<-pi.inf(exp(xbeta.star)/(1+exp(xbeta.star)))
if(family=="gaussian") pistar<-xbeta.star
Y.HT<-Rc*Y2c/pistar
mu.BR.HT.star[dd]<-mean(Y.HT)-mu.BR.HT
}

# End PS Model resampling

# Next, the outcome response (OR) model
xbeta<-xx.OR.obs%*%b.OR
if(family=="binomial") m.b<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") m.b<-xbeta
resid<-as.vector(Y2c-m.b)
temp<-xx.OR.obs*Rc*resid*gg                       # (n x p) matrix
q.star.OR<-apply(temp,2,mean)

if(checking) cat("OR Model Parzen Resampling Draw=",c(dd),"\n")

b.start<-b.OR
fit.NR<-QLE.NR(b.start=b.OR,Y=Y2c,X=as.matrix(xx.OR.obs),q.star=q.star.OR,wt=Rc,
family=family, abstol=abstol,checking=checking)

if(!is.null(fit.NR)){
bstar.OR<-fit.NR$beta
Qsoln<-fit.NR$norm
NR.iter<-fit.NR$iter
converged.OR<-TRUE
if(checking) cat("Good: OR Converged via N-R (Qsoln,iter)",c(Qsoln,NR.iter),"\n")
}

if(is.null(fit.NR)){
b.start<-b.OR
temp<-rescue.resample.QLE(b.start=b.start,Y=Y2c,X=as.matrix(xx.OR.obs),
family=family,q.star=q.star.OR,wt=Rc,abstol=abstol,checking=checking,draw=dd)
bstar.OR<-temp$bstar
converged.OR<-temp$converged
}

if(converged.OR){
xbeta.star<-xx.OR.obs%*%bstar.OR
if(family=="binomial") ehat.OR.star<-pi.inf(exp(xbeta.star)/(1+exp(xbeta.star)))
if(family=="gaussian") ehat.OR.star<-xbeta.star
mu.BR.OR.star[dd]<-mean(ehat.OR.star)-mu.BR.OR
}

# End OR Model resampling

if(converged.PS){

if(checking) cat("PS Model Converged -- Resampling OR-PS Model Draw=",c(dd),"\n")
########################################################
#                 BR-DR estimator
########################################################
temp<-xx.ORPS.obs*as.vector(Rc*(Y2c-e.ORPS.obs)*gg)
term1<-apply(temp,2,mean)

xx.ORPS.star<-as.matrix(cbind(rep(1,nc),V.qc,1/pistar))
xx.diff<-xx.ORPS.star-xx.ORPS.obs
temp<-xx.diff*as.vector(Rc*(Y2c-e.ORPS.obs))
term2<-apply(temp,2,mean)

temp<-xx.ORPS.obs*as.vector(Rc*(pihat-pistar))
term3<-apply(temp,2,mean)
q.ORPS.star<-term1+term2+term3

b.start<-b.ORPS
fit.NR<-QLE.NR(b.start=b.start,Y=Y2c,X=as.matrix(xx.ORPS.obs),q.star=q.ORPS.star,wt=Rc,
family=family,abstol=abstol,checking=checking)

if(!is.null(fit.NR)){
b.ORPS.star<-fit.NR$beta
Qsoln<-fit.NR$norm
NR.iter<-fit.NR$iter
converged.ORPS<-TRUE
if(checking) cat("Great: OR-PS Converged via N-R (Qsoln,iter)",c(Qsoln,NR.iter),"\n")
}

if(is.null(fit.NR)){
b.start<-b.ORPS
temp<-rescue.resample.QLE(b.start=b.start,Y=Y2c,X=xx.ORPS.obs,q.star=q.ORPS.star,wt=Rc,
family=family, abstol=abstol,checking=checking,draw=dd)
b.ORPS.star<-temp$bstar
converged.ORPS<-temp$converged
}

# If this last guy still does not converge try again 
if(!converged.ORPS){
if(checking) cat("Trying second iteration of rescue for OR-PS Estimator","\n")
b.start<-temp$b.combo
abstol2<-2*abstol # Decrease the convergence hurdle

if(checking) cat("Initial value=",c(b.start),"\n")

temp<-rescue.resample.QLE(b.start=b.start,Y=Y2c,X=xx.ORPS.obs,q.star=q.ORPS.star,wt=Rc,
family=family, abstol=abstol2,checking=checking,draw=dd)
b.ORPS.star<-temp$bstar
converged.ORPS<-temp$converged

#logit.star<-optim(par=c(b.start),fn=QLE.star,X=xx.ORPS.obs,Y=Y2c,wt=Rc,q.star=q.ORPS.star,family=family,
#control=list(abstol=abstol), method="SANN")
#b.ORPS.star<-logit.star$par
#Qsoln<-logit.star$value
#converged.ORPS<-(logit.star$convergence==0 & Qsoln<=abstol)
#if(checking & converged.ORPS) cat("Great: OR-PS Converged via N-R (Qsoln)",c(Qsoln),"\n")

}

if(converged.ORPS){
xbeta.star<-xx.ORPS.obs%*%b.ORPS.star
if(family=="binomial") eORPS.star<-pi.inf(exp(xbeta.star)/(1+exp(xbeta.star)))
if(family=="gaussian") eORPS.star<-xbeta.star
mu.BR.DR.star[dd]<-mean(eORPS.star)-mu.BR.DR     
}

if(!checking) cat("Parzen sample=",c(dd),"\n")

} # End OR-PS Model resampling

} # End Parzen Resampling

if(se.method=="Bootstrap"){ # Start Bootstrap
id.boot<-sample(data.id,size=nc,replace=T)

Rc.boot<-data.obs$Rc[id.boot]
Y2c.boot<-data.obs$Y2c[id.boot]
Vpi.boot<-data.obs$V.pic[id.boot,]
Vq.boot<-data.obs$V.qc[id.boot,]

fit.pi.boot<-glm(Rc.boot~as.matrix(Vpi.boot), family=family)
pihat.boot<-fitted.values(fit.pi.boot)

################ Next, implement the Bang-Robins DR Estimator #####################
#          Propensity-Score Model is already fitted -- pihat0

Y.HT.boot<-Rc.boot*Y2c.boot/pihat.boot

mu.BR.HT.boot<-mean(Y.HT.boot)  # Horvitz-Thompson Estimator

# OR-PS model
VORPS.boot<-cbind(Vq.boot,1/pihat.boot) # Adding propensity-score to OR model
case.boot<-(Rc.boot==1)
b.ORPS.boot<-as.vector(glm(Y2c.boot[case.boot]~as.matrix(VORPS.boot[case.boot,]), family=family)$coefficients)

z<-as.matrix(cbind(rep(1,nc),VORPS.boot))
xbeta<-z%*%b.ORPS.boot
if(family=="binomial") eORPS.boot<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") eORPS.boot<-xbeta
mu.BR.DR.boot<-mean(eORPS.boot)              # Bang-Robins Double-Robust Estimator

mu.BR.HT.star[dd]<-mu.BR.HT.boot-mu.BR.HT
mu.BR.DR.star[dd]<-mu.BR.DR.boot-mu.BR.DR

mu.BR.DR.noncase.star[dd]<-(sum(eORPS.boot[!case.boot])/nc)-mu.BR.DR.noncase

b.OR.boot<-as.vector(glm(Y2c.boot[case.boot]~as.matrix(Vq.boot[case.boot,]), family=family)$coefficients)

z<-cbind(rep(1,nc),Vq.boot) 
xbeta<-z%*%b.OR.boot
if(family=="binomial") ehat.OR.boot<-pi.inf(exp(xbeta)/(1+exp(xbeta)))
if(family=="gaussian") ehat.OR.boot<-xbeta

# B-R OR model estimator

mu.BR.OR.boot<-mean(ehat.OR.boot)

mu.BR.OR.star[dd]<-mu.BR.OR.boot-mu.BR.OR

if(checking) cat("Bootstrap sample=",c(dd),"\n")

}  # Finish Bootstrap

}

if(checking){ 
par(mfrow=c(2,2)) # Plot histogram using trimming
names(mu.BR.DR.star)<-c("DR Estimator")
names(mu.BR.HT.star)<-c("HT Estimator")
names(mu.BR.OR.star)<-c("OR Estimator")
hist.trim(mu.BR.DR.star,trim=trim.DR)
hist.trim(mu.BR.HT.star,trim=trim.HT)
hist.trim(mu.BR.OR.star,trim=trim.OR)
}

par(mfrow=c(1,1))

#var.BR.HT<-var.trim(mu.BR.HT.star,trim=trim.HT) 
#var.BR.DR<-var.trim(mu.BR.DR.star,trim=trim.DR)
#var.BR.OR<-var.trim(mu.BR.OR.star,trim=trim.OR)

# Robust variance estimators

var.BR.HT<-(mad(mu.BR.HT.star))^2
var.BR.DR<-(mad(mu.BR.DR.star))^2
var.BR.OR<-(mad(mu.BR.OR.star))^2

}

# Let's check computational details
if(checking){
# Sum of complete case
sum.CC<-sum(Y2c[case])
sum.OR<-sum(e.OR[case])
sum.ORPS<-sum(e.ORPS[case])
cat("Group=",c(unique(Zc)),"\n")
cat("Total contribution from cases (CC,OR,ORPS)",c(sum.CC,sum.OR,sum.ORPS),"\n")
cat("Total contribution from non-cases (CC,OR,ORPS)",c(0,sum(e.OR[!case]),sum(e.ORPS[!case])),"\n")
}

print(mean(mu.BR.DR.noncase.star))

return(list(mu.BR.DR=mu.BR.DR,var.BR.DR=var.BR.DR,
mu.BR.HT=mu.BR.HT,var.BR.HT=var.BR.HT,
mu.BR.OR=mu.BR.OR,var.BR.OR=var.BR.OR,
mu.BR.DR.star=mu.BR.DR.star,
mu.BR.OR.star=mu.BR.OR.star,
mu.BR.HT.star=mu.BR.HT.star,
mu.BR.DR.noncase=mu.BR.DR.noncase,
mu.BR.DR.noncase.star=mu.BR.DR.noncase.star,
pihat=pihat,e.OR=e.OR, e.ORPS=e.ORPS,
Z=Zc,Y2=Y2c,R=Rc,
ehat.h=ehat.h,ehat.h.all=ehat.h.all,
ehat.q=ehat.q,ehat.q.all=ehat.q.all))
}
