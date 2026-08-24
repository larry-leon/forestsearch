# This is Hardle's simulated example pg. 124

rm(list=ls())

source("../R/fan-gijbels_transform_functions.R")
source("../R/kernel-regression_functions.R")
source("../R/km_functions_twosample_weighting.R")

hardle.data<-function(n,seed,censor=FALSE,cmin=-3,cmax=2){
  set.seed(seed)
    x<-sort(runif(n))
  m<-(sin(2*pi*x^3))^3
  e<-rnorm(n)*sqrt(0.1)
  y.1<-m+e
  y.2<-m+e*x
  y.1.cens<-y.2.cens<-NULL
  event1<-event2<-NULL
  if(censor){
    # log-scale censoring
    c<-runif(n=n,min=cmin,max=cmax)
    y.1.cens<-pmin(y.1,c)
    y.2.cens<-pmin(y.2,c)
    event1<-ifelse(y.1<=c,1,0)
    event2<-ifelse(y.2<=c,1,0)
  }
  
  return(list(x=x,y=cbind(y.1,y.2),y.cens=cbind(y.1.cens,y.2.cens),m=m,e=e,event=cbind(event1,event2)))
}


hardle.sims<-function(n){
  x<-sort(runif(n))
  m<-(sin(2*pi*x^3))^3
  e<-rnorm(n)*sqrt(0.1)
  y.1<-m+e
  y.2<-m+e*x
  return(list(x=x,y=cbind(y.1,y.2),m=m,e=e))
}


binom.sims<-function(n){
  x<-sort(runif(n))
  m<-(sin(2*pi*x^3))^3
  p<-exp(m)/(1+exp(m))
  #e<-rnorm(n)*sqrt(0.1)
  #y.1<-m+e
  #y.2<-m+e*x
  y<-rbinom(n,size=1,prob=p)
  return(list(x=x,y=y,p=p))
}


seed<-c(13,22,28,13,47,0,45,38,15,37,28,2)


dat.test<-binom.sims(200)

plot(dat.test$x,dat.test$y)
lines(dat.test$x,dat.test$p)


dat.test<-hardle.data(256,seed)

x<-as.vector(dat.test$x)
y<-as.matrix(dat.test$y)
m<-dat.test$m

xx<-sort(x)

## Get a few bandwidth Estimates

hh<-0.10
# Note Biquad is Hardle's quartic Kernel
test.kernel<-m.Kernel(z=xx,Z=x,Y=y,Kernel="Biquad",bandwidth=hh)

m.true<-dat.test$m
id.x<-order(x)
m.true<-m.true[id.x]
yy<-y[id.x,]

m.fit<-test.kernel$m.Kernel
f.fit<-test.kernel$f.Kernel

# h=0.01, h=0.05, h=0.10
m.0<-m.true
y.obs<-yy

#m.fit1<-m.fit; f.fit1<-f.fit
#m.fit2<-m.fit; f.fit2<-f.fit
#m.fit3<-m.fit; f.fit3<-f.fit


ymin<-min(c(yy,m.true,m.fit))
ymax<-max(c(yy,m.true,m.fit))

par(mfrow=c(1,2))
plot(xx,yy[,1],type="p",xlab="x",ylab="Hardle simulated data")
lines(xx,m.true,lty=1,lwd=1.5,col="blue")
lines(xx,m.fit[1,],lty=2,col="red")

plot(xx,yy[,2],type="p",xlab="x",ylab="Hardle simulated data")
lines(xx,m.true,lty=1,lwd=1.5,col="blue")
lines(xx,m.fit[2,],lty=2,col="red")



seed<-8316954
dat.test<-hardle.data(1000,seed=seed,censor=TRUE,cmin=-5,cmax=5)

x<-as.vector(dat.test$x)
y<-as.matrix(dat.test$y.cens)
event<-as.matrix(dat.test$event)
xx<-sort(x)

test.kernel<-m.Kernel(z=xx,Z=x,Y=y,Kernel="Biquad",bandwidth=hh)

m.true<-dat.test$m
id.x<-order(x)
m.true<-m.true[id.x]
yy<-y[id.x,]

m.fit<-test.kernel$m.Kernel
f.fit<-test.kernel$f.Kernel

# h=0.01, h=0.05, h=0.10
m.0<-m.true
y.obs<-yy

ymin<-min(c(yy,m.true,m.fit))
ymax<-max(c(yy,m.true,m.fit))

par(mfrow=c(1,2))

plot(xx,yy[,1],type="p",xlab="x",ylab="Hardle simulated data")
lines(xx,m.true,lty=1,lwd=1.5,col="blue")
lines(xx,m.fit[1,],lty=2,col="grey",lwd=4)

plot(xx,yy[,2],type="p",xlab="x",ylab="Hardle simulated data")
lines(xx,m.true,lty=1,lwd=1.5,col="blue")
lines(xx,m.fit[2,],lty=2,col="grey",lwd=4)


# Apply to heteroskedastic outcomes 
par(mfrow=c(1,1))

df<-data.frame(x=dat.test$x,time=dat.test$y.cens[,2],
               event.time=dat.test$event[,2])

data.new<-get.FG(df)

x<-as.vector(data.new$x)
y<-as.matrix(data.new$Y.FG)
xx<-sort(x)

test.kernel<-m.Kernel(z=xx,Z=x,Y=y,Kernel="Biquad",bandwidth=hh)

m.true<-dat.test$m
id.x<-order(x)
m.true<-m.true[id.x]
yy<-y[id.x,]

m.fit<-test.kernel$m.Kernel
f.fit<-test.kernel$f.Kernel

# h=0.01, h=0.05, h=0.10
m.0<-m.true
y.obs<-yy

ymin<-min(c(yy,m.true,m.fit))
ymax<-max(c(yy,m.true,m.fit))

plot(xx,yy,type="p",lwd=0.25,xlab="x",ylab="Hardle simulated data")
lines(xx,m.true,lty=1,lwd=1.5,col="blue")
lines(xx,m.fit,lty=2,col="grey",lwd=4)

print(summary(m.fit-m.true))



