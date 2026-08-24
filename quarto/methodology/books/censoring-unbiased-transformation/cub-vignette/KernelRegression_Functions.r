###########################################
#           Kernels 
###########################################
w.Uniform<-function(x){
1/2*(x>=-1 & x<=1)
}
w.Ep<-function(x){
3/4*(1-x^2)*(x>=-1 & x<=1)
}
w.Biweight<-function(x){
15/16*(1-x^2)*(x>=-1 & x<=1)
}
w.Symmetric<-function(x){
35/32*((1-x^2)^3)*(x>=-1 & x<=1)
}
w.TriCube<-function(x){
((1-abs(x)^3)^3)*(x>=-1 & x<=1)
}
w.BiCube<-function(x){
((1-abs(x)^2)^2)*(x>=-1 & x<=1)
}
w.BiQuad<-function(x){
(15/16*((1-x^2)^2))*(x>=-1 & x<=1)
}
w.M4<-function(x){
(105/64)*(1-5*x^2+7*x^4-3*x^6)*(x>=-1 & x<=1)
}
w.M6<-function(x){
(315/2048)*(15-140*x^2+378*x^4-396*x^6+143*x^8)*(x>=-1 & x<=1)
}
w.Twice<-function(x){
(2*exp(-(x^2)/2)/sqrt(2*pi))-(exp(-(x^2)/4)/sqrt(4*pi))
}

w.Kernel<-function(x,K.type){
switch(K.type,
Normal=dnorm(x,sd=1),
Epan=w.Ep(x),
Biweight=w.Biweight(x),
Symmetric=w.Symmetric(x),
Tricube=w.TriCube(x),
Uniform=w.Uniform(x),
Bicube=w.BiCube(x),
Biquad=w.BiQuad(x),
M4=w.M4(x),
M6=w.M6(x),
Twice=w.Twice(x))
}

#########################################################
# Computes Kernel Regression where Y may be multivariate
# but Z is univariate
####################################################
# x is ordered points at which to compute m()
# This gets the estimate at points, z.
mz.Kernel<-function(z,Z,Y,Kernel,bandwidth){
h<-bandwidth
if(!is.matrix(Y)) stop("Y should be matrix")
if(!is.vector(Z)) stop("Z should be vector")
p<-ncol(Y)
NW.z<-as.vector(w.Kernel(-(Z-z)/h,K.type=Kernel)/h) # Nadaraya-weight at z
f.z<-mean(NW.z)
if(f.z==0){
#cat("Zero density estimate at (z,h)", c(z,h),"\n")
m.z<-rep(0,length=p)
return(c(z,f.z,m.z))
}
if(f.z>0){
y.wz<-Y*NW.z  # Weight each column of Y by NW.z 
if(p>1) r.z<-apply(y.wz,2,mean)
if(p==1) r.z<-mean(y.wz)
m.z<-r.z/f.z
# Can also use m.z<-apply(Y,2,weighted.mean,w=NW.z)
return(c(z,f.z,m.z))
}
}

nomis<-function(x){x[!is.na(x)]}

m.Kernel<-function(z,Z,Y,Kernel="Uniform",bandwidth,is.sorted=F){
# Biquad is the quartic Kernel
# Epan is the Epanechnikov Kernel
K.valid<-is.element(Kernel,c("Normal","Epan","Biweight","Symmetric","Tricube","Uniform",
"Bicube","Biquad","M4","M6","Twice"))
if(!K.valid) stop("Kernel is not valid")
# First flag any missings
if(!is.vector(Z)) stop("Z should be univariate vector")
#k.z<-ncol(Z) 
#if(k.z!=1) stop("Z is not univariate")
z.mis<-is.na(Z)
y.mis<-is.na(Y)
missingflag<-sum(z.mis)+sum(y.mis)
if(missingflag>0) warning("missing values in (Z,Y): m.Kernel call")
Z<-as.vector(Z[!z.mis])
#Y<-as.matrix(Y[!y.mis])
Y<-as.matrix(apply(Y,2,nomis))
n<-nrow(Y)
# Sort data
if(!is.sorted){
id.z<-order(Z); Z<-Z[id.z]; Y<-as.matrix(Y[id.z,])
}
Kernel.est<-sapply(z,mz.Kernel,Z=Z,Y=Y,Kernel=Kernel,bandwidth=bandwidth)
# Note the first component returned is z, second is f.z, third + is m.z
zz<-unlist(Kernel.est[1,])
f.z<-unlist(Kernel.est[2,])
p<-ncol(Y)
if(p>1) m.z<-unlist(Kernel.est[c(3:(p+2)),])
if(p==1) m.z<-unlist(Kernel.est[3,])
#if(p>1) m.z<-matrix(temp,ncol=p,nrow=n,byrow=T)
return(list(at.points=zz,f.Kernel=f.z,m.Kernel=m.z))
}



rz.Kernel.star<-function(z,Z,D,Kernel,bandwidth){
h<-bandwidth
if(!is.vector(D)) stop("D should be vector")
if(!is.vector(Z)) step("Z should be vector")
NW.z<-as.vector(w.Kernel(-(Z-z)/h,K.type=Kernel)/h) # Nadaraya-weight at z
d.wz<-D*NW.z  
r.z<-mean(d.wz)
return(c(r.z))
}


r.Kernel.star<-function(G,z,Z,D,Kernel="Uniform",bandwidth,is.sorted=F){
# Biquad is the quartic Kernel
# Epan is the Epanechnikov Kernel
K.valid<-is.element(Kernel,c("Normal","Epan","Biweight","Symmetric","Tricube","Uniform",
"Bicube","Biquad","M4","M6","Twice"))
Y<-D*(G-1)
if(!K.valid) stop("Kernel is not valid")
# First flag any missings
if(!is.vector(Z)) stop("Z should be univariate vector")
if(!is.vector(Y)) stop("Y should be univariate vector")
z.mis<-is.na(Z)
y.mis<-is.na(Y)
missingflag<-sum(z.mis)+sum(y.mis)
if(missingflag>0) warning("missing values in (Z,Y): r.Kernel call")
Z<-as.vector(Z[!z.mis])
Y<-as.vector(Y[!y.mis])
n<-length(Y)
# Sort data
if(!is.sorted){
id.z<-order(Z); Z<-Z[id.z]; Y<-as.vector(Y[id.z])
}
r.Kernel<-unlist(lapply(z,rz.Kernel.star,Z=Z,D=Y,Kernel=Kernel,bandwidth=bandwidth))
return(c(mean(r.Kernel)*sqrt(n)))
}


ise.star<-function(z,Z,D,Kernel,bandwidth,draws){
n<-length(D)
G.mat<-matrix(rexp(n*draws),nrow=n,ncol=draws)
ise<-apply(G.mat,2,r.Kernel.star,z=z,Z=Z,D=D,Kernel=Kernel,bandwidth=bandwidth)
return(c(ise))
}


Hardle.data<-function(n,seed){
.Random.seed<-seed
x<-sort(runif(n))
m<-(sin(2*pi*x^3))^3
e<-rnorm(n)*sqrt(0.1)
y.1<-m+e
#y.2<-m+e*x
return(list(x=x,y=y.1,m=m,e=e))
}

Hardle.sims<-function(n){
x<-sort(runif(n))
m<-(sin(2*pi*x^3))^3
e<-rnorm(n)*sqrt(0.1)
y.1<-m+e
#y.2<-m+e*x
return(list(x=x,y=y.1,m=m,e=e))
}

binom.sims<-function(n){
x<-sort(runif(n))
m<-(sin(2*pi*x^3))^3
p<-exp(m)/(1+exp(m))
y<-rbinom(n,size=1,prob=p)
return(list(x=x,y=y,p=p))
}


binomial.sim<-function(n,draws,sims,bandwidth,Kernel="Biquad"){
dat.test<-binom.sims(n)
x<-as.vector(dat.test$x)
y<-as.matrix(dat.test$y)
m<-dat.test$p
xx<-sort(x)
test.kernel<-m.Kernel(z=xx,Z=x,Y=y,Kernel=Kernel,bandwidth=bandwidth)
id.x<-order(x)
m.true<-m[id.x]
yy<-y[id.x,]
x.obs<-xx
y.obs<-yy
m.obs<-m.true
m.fit<-test.kernel$m.Kernel
f.fit<-test.kernel$f.Kernel
mf.obs<-m.fit
D.1<-(yy-m.fit)/f.fit
Ise.star<-ise.star(z=xx,Z=x,D=as.vector(D.1),Kernel=Kernel,bandwidth=bandwidth,draws=draws)
dens.ise<-density(Ise.star,kernel="gaussian")
set.seed(831695)
ise.sims<-rep(NA,sims)
for(ss in 1:sims){
sim.dat<-binom.sims(n)
x<-as.vector(sim.dat$x)
y<-as.matrix(sim.dat$y)
m.sim<-sim.dat$p
xxx<-sort(x)
test.kernel<-m.Kernel(z=xxx,Z=x,Y=y,Kernel=Kernel,bandwidth=bandwidth)
id.x<-order(x)
m.sim<-m.sim[id.x]
m.fit<-test.kernel$m.Kernel
ise.sims[ss]<-sqrt(n)*mean(m.fit-m.sim)
}
dens.ise.sims<-density(ise.sims,kernel="gaussian")
xx.1<-dens.ise$x; yy.1<-dens.ise$y
xx.2<-dens.ise.sims$x; yy.2<-dens.ise.sims$y
result<-list(resamp.x=xx.1,resamp.dens=yy.1,sim.x=xx.2,sim.dens=yy.2,
x.obs=x.obs,y.obs=y.obs,m.obs=m.obs,mfit.obs=mf.obs)
return(result)
}



Hardle.sim<-function(n,draws,sims,bandwidth,Kernel="Biquad"){
seed<-c(13,22,28,13,47,0,45,38,15,37,28,2)
dat.test<-Hardle.data(n,seed=seed)
x<-as.vector(dat.test$x)
y<-as.matrix(dat.test$y)
m<-dat.test$m
xx<-sort(x)
test.kernel<-m.Kernel(z=xx,Z=x,Y=y,Kernel=Kernel,bandwidth=bandwidth)
id.x<-order(x)
m.true<-m[id.x]
yy<-y[id.x,]
x.obs<-xx
y.obs<-yy
m.obs<-m.true
m.fit<-test.kernel$m.Kernel
f.fit<-test.kernel$f.Kernel
mf.obs<-m.fit
D.1<-(yy-m.fit)/f.fit
Ise.star<-ise.star(z=xx,Z=x,D=as.vector(D.1),Kernel=Kernel,bandwidth=bandwidth,draws=draws)
dens.ise<-density(Ise.star,kernel="gaussian")
set.seed(831695)
ise.sims<-rep(NA,sims)
for(ss in 1:sims){
sim.dat<-Hardle.sims(n)
x<-as.vector(sim.dat$x)
y<-as.matrix(sim.dat$y)
m.sim<-sim.dat$m
xxx<-sort(x)
test.kernel<-m.Kernel(z=xxx,Z=x,Y=y,Kernel=Kernel,bandwidth=bandwidth)
id.x<-order(x)
m.sim<-m.sim[id.x]
m.fit<-test.kernel$m.Kernel
ise.sims[ss]<-sqrt(n)*mean(m.fit-m.sim)
}
dens.ise.sims<-density(ise.sims,kernel="gaussian")
xx.1<-dens.ise$x; yy.1<-dens.ise$y
xx.2<-dens.ise.sims$x; yy.2<-dens.ise.sims$y
result<-list(resamp.x=xx.1,resamp.dens=yy.1,sim.x=xx.2,sim.dens=yy.2,
x.obs=x.obs,y.obs=y.obs,m.obs=m.obs,mfit.obs=mf.obs)
return(result)
}


# Cross-validation for estimating v(x)=E[V|x]
CV.h<-function(bandwidth,V,x,Kernel,Fx,GCV=F,trim=0.05,status=F){
# For each i compute hat(v)_[-i]
n<-length(x)
Vx.i<-matrix(NA,nrow=nrow(V),ncol=ncol(V))
for(ii in 1:length(x)){
Vii<-as.matrix(V[-c(ii),])
xii<-x[-c(ii)]
x.out<-x[c(ii)] # We are predicting this guy, using xii and Vii
v.kernel<-m.Kernel(z=x.out,Z=as.vector(xii),Y=as.matrix(Vii),
Kernel=Kernel,bandwidth=bandwidth,is.sorted=T)$m.Kernel
Vx.i[ii,]<-c(v.kernel)
}
PE.h<-as.matrix(V-Vx.i)
xx.low<-min(x[Fx>=trim])
xx.up<-max(x[Fx<=(1-trim)])
w.x<-as.numeric(x>=xx.low & x<=xx.up)
w.GCV<-1.0
if(GCV){
K.0<-w.Kernel(0.0,K.type=Kernel)
w.GCV<-(1-((n*bandwidth)^(-1))*K.0)^(-2)
}
cv.h<-mean(PE.h^2*w.x)*w.GCV
if(status) print(c(bandwidth,cv.h))
#cv.h<-sqrt(t(PE.h)%*%PE.h)
#cv.h<-sum(apply((PE.h^2)*w.x,2,mean))
return(cv.h)
}
