
expit<-function(x){
  exp(x)/(1+exp(x))
}


x.truncate<-function(x,truncate){
  x.quant<-quantile(x,c(truncate,1-truncate))
  x.trunc<-x
  x.trunc[which(x<=x.quant[1])]<-x.quant[1]
  x.trunc[which(x>=x.quant[2])]<-x.quant[2]
  return(x.trunc)
}


get.ps.weights<-function(data,PS.model.fit,truncate=0){
  fit.ps<-glm(PS.model.fit,family="binomial",data=data)
  pihat.Vps<-fit.ps$fitted
  wt.1<-1/pihat.Vps
  wt.0<-1/(1-pihat.Vps)
  ipw.weights<-ifelse(data$Exposure==1,wt.1,wt.0)
  
  if(truncate==0) ipw.truncate<-ipw.weights
  if(truncate>0) ipw.truncate<-x.truncate(ipw.weights,truncate)
  
  pihat.null<-glm(Exposure~1,family="binomial",data=data)$fitted
  # Stabilized
  wt.1<-pihat.null/pihat.Vps
  wt.0<-(1-pihat.null)/(1-pihat.Vps)
  sw.weights<-ifelse(data$Exposure==1,wt.1,wt.0)
  return(list(ipw.weights=ipw.weights,sw.weights=sw.weights,ipw.trunc.weights=ipw.truncate))
}



# Modified from Personalized package

check.overlap<-function(data,type="histogram",wt.name="ipw",treat.name="Exposure", bins = 50L, alpha=0.5){
  pi.wt<-as.vector(data[,c(wt.name)])  
  trt<-data[,c(treat.name)]
  prop.scores<-data.frame(Treatment = as.factor(trt),prop.score = pi.wt)
  Treatment <- prop.score <- NULL
  
  if (type == "density") {
    
    pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) + 
      geom_density(alpha = alpha, colour = "grey50") + 
      geom_rug(aes(colour = Treatment)) + theme(legend.position = "bottom") + 
      ggtitle("Densities of propensity scores by treatment group") + 
      xlab("Propensity Score")
    
  }
  else if (type == "histogram") {
    pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) + 
      geom_histogram(bins = bins, alpha = alpha, position = "identity") + 
      geom_rug(aes(colour = Treatment)) + theme(legend.position = "bottom") + 
      ggtitle("Histograms of propensity scores by treatment group") + 
      xlab("Propensity Score")
  }
  
  pl.obj    
}



# The full data ("final")
# Note: assuming in ymd format
interim.cut<-function(data,time.name,event.name,rand.name,tau.date,details=FALSE){
Y.final<-data[,c(time.name)]
Delta.final<-data[,c(event.name)]
if(details) cat("# of events for full follow-up",c(sum(Delta.final)),"\n")
rand.date<-ymd(data[,c(rand.name)])
start.date<-min(rand.date)
F.tau<-as.numeric(tau.date-rand.date+1)/30.4375 # This is the follow-up per patient
Y.tau<-pmin(F.tau,Y.final)
Delta.tau<-ifelse(((Y.final<=F.tau) & Delta.final==1),1,0)
if(details) cat("# of events for interim",c(sum(Delta.tau)),"\n")
data.tau<-data
data.tau[,c(time.name)]<-Y.tau
data.tau[,c(event.name)]<-Delta.tau
# The interim data would be those who were enrolled by the analysis time tau.date
data.tau<-subset(data.tau,F.tau>0)
return(data.tau)
}


plot.band<-function(x,mean.value,lower,upper,show.axes=F,band=TRUE,ltype="l",lty=1,xlabel=NULL,ylabel=NULL,color="grey",ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
  plot(x[order(x)],mean.value[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
       ylab=ylabel,ylim=ylim)
  if(band) polygon(c(x[order(x)],rev(x[order(x)])),c(lower[order(x)],rev(upper[order(x)])),col=color,border=FALSE)
  lines(x[order(x)],mean.value[order(x)],lty=lty,lwd=2.5,type=ltype)
}


plot.band.two<-function(x,curve1,curve2,lower,upper,show.axes=F,
                        ltype="l",lty=1,xlabel=NULL,ylabel=NULL,ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
  plot(x[order(x)],curve1[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
       ylab=ylabel,ylim=ylim)
  polygon(c(x[order(x)],rev(x[order(x)])),
          c(lower[order(x)],rev(upper[order(x)])),col="lightgrey",border=F)
  lines(x[order(x)],curve1[order(x)],lty=lty,lwd=2.5,type=ltype)
  lines(x[order(x)],curve2[order(x)],lty=lty,lwd=2.5,type=ltype)
}


N.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error<=x))
}

R.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error>=x))
}

NA.CHR.Weighted<-function(time,Delta,W.n=rep(1,length(time)),W.d=rep(1,length(time)),
                          at.points=sort(time),se.type="greenwood",get.Stute=FALSE,tpoints.add=NULL){
  
  if(!is.null(tpoints.add)) at.points<-sort(c(unique(c(at.points,tpoints.add))))
  
  if(se.type!="greenwood" & se.type!="tsiatis") stop("Invalid se type -- greenwood or tsiatis allowed")
  #is.sorted<-(all(time==sort(time)))
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id<-order(time); time<-time[id]; Delta<-Delta[id]; W.n<-W.n[id]; W.d<-W.d[id]
  }
  risk<-unlist(lapply(as.list(at.points),R.Weighted,error=time,W=W.d))
  ###########################################################################
  ### Adaptive H process correspoinding to N-A rep. via integral wrt M.G ####
  Hmart.chf<-ifelse(risk>0,1/risk,0)
  ############################################################################
  counting<-unlist(lapply(as.list(at.points),N.Weighted,error=time,W=W.n*ifelse(Delta==1,1,0)))
  counting <- c(0, counting)
  dN<-diff(counting)
  dN.risk<-ifelse(risk>0,dN/risk,0.0)
  chf <- cumsum(dN.risk)
  var.chf<-cumsum(ifelse(risk>0,dN/(risk^2),0.0))
  S.KM <- cumprod(1-dN.risk)
  
  S.KM[which(S.KM<0)]<-0.0
  
  S.NA <- exp(-chf)
  var.NA<-(S.NA^2)*var.chf
  # Greenwood variance estimate
  if(se.type=="greenwood"){
    aa<-dN
    bb<-risk*(risk-dN)
    var.KM<-(S.KM^2)*cumsum(ifelse(risk>0,aa/bb,0.0))
    se.KM<-sqrt(var.KM)
  }
  if(se.type=="tsiatis"){
    var.KM<-(S.KM^2)*var.chf
    se.KM<-sqrt(var.KM)
  }
  result<-list(time=time,at.points=at.points,S.NA=S.NA,S.KM=S.KM,chf=chf,se.chf=sqrt(var.chf),
               se.NA=sqrt(var.NA),dN.risk=dN.risk,
               n.risk=risk,dN=dN,
               Hmart.chf=Hmart.chf,se.KM=se.KM)
  return(result)
}


Logrank<-function(time,Delta,X){
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id<-order(time); time<-time[id]; Delta<-Delta[id]; X<-X[id]
  }
  at.points<-sort(unique(c(time[Delta==1])))
  
  U0<-time[which(X==0)]
  D0<-Delta[which(X==0)]
  
  # Control group
  # Risk and Counting processes
  risk.z0<-unlist(lapply(as.list(at.points),R.Weighted,error=U0))
  counting<-unlist(lapply(as.list(at.points),N.Weighted,error=U0[D0==1]))
  N.z0<-counting
  dN.z0<- diff(c(0, counting))
  
  U1<-time[which(X==1)]
  D1<-Delta[which(X==1)]
  
  # Control group
  # Risk and Counting processes
  risk.z1<-unlist(lapply(as.list(at.points),R.Weighted,error=U1))
  counting<-unlist(lapply(as.list(at.points),N.Weighted,error=U1[D1==1]))
  N.z1<-counting
  dN.z1<- diff(c(0, counting))
  
  dN.pooled<-dN.z0+dN.z1
  risk.pooled<-risk.z0+risk.z1
  
  K<-(risk.z0*risk.z1)/(risk.pooled)
  
  term0<-sum(ifelse(risk.z0>0,(K/risk.z0)*dN.z0,0.0))
  term1<-sum(ifelse(risk.z1>0,(K/risk.z1)*dN.z1,0.0))
  
  lr<-term0-term1
  
  # variance
  h0<-ifelse(risk.z0==0,0,(K^2/risk.z0))
  h1<-ifelse(risk.z1==0,0,(K^2/risk.z1))
  dJ<-ifelse(risk.pooled==1,0,(dN.pooled-1)/(risk.pooled-1))
  dL<-ifelse(risk.pooled==0,0,dN.pooled/risk.pooled)
  sig2s<-(h0+h1)*(1-dJ)*dL
  sig2<-sum(sig2s)
  
  Z.lr<-lr/sqrt(sig2)
  result<-list(Z.lr=Z.lr,lr=lr,pval=1-pchisq(Z.lr^2,1))
  return(result)
}

