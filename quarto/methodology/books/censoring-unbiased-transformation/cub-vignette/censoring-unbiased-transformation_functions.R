
################ Functions ####################
Less.equal<-function(X){
  length(X[X<=X[1]])-1
}

#Ties<-function(X){
#length(X[X==X[1]])-1
#}

Jump.points<-function(X){
  cbind(sort(unique(X)),order(unique(X)))
}
Greater.equal<-function(X){
  length(X[X>=X[1]])-1
}


E.CDF<-function(Y){
  jump.points<-Jump.points(Y)
  print(jump.points)
  Temp<-matrix(rep(Y,times=length(jump.points)),nrow=length(jump.points),byrow=T)
  Time.mat<-cbind(jump.points,Temp)
  CDF<-apply(Time.mat,1,Less.equal)/length(Y)
  return(cbind(jump.points,CDF))
}


NA.CHR<-function(X,D){
  #######################################
  #jump.points<-Jump.points(X[D==1])[,1]
  # Want to calculate at all points
  # Already sorted by time
  #jump.points<-sort(X)
  #ids<-order(X)
  #censor<-D[ids]
  #z<-Z[ids,]
  jump.points<-X
  censor<-D
  ########################################
  Time.mat<-matrix(rep(X,times=length(jump.points)),nrow=length(jump.points),byrow=T)
  Times.mat<-cbind(jump.points,Time.mat)
  risk<-apply(Times.mat,1,Greater.equal)
  events<-X[D==1]
  Time.mat<-matrix(rep(events,times=length(jump.points)),nrow=length(jump.points),byrow=T)
  Times.mat<-cbind(jump.points,Time.mat)
  counting<-apply(Times.mat,1,Less.equal)
  
  n.risk<-risk
  counting<-c(0,counting)
  delta.counting<-diff(counting)
  risk.points<-risk
  DN.Risk<-(delta.counting/risk.points)
  n.event<-delta.counting
  chr.NA<-cumsum(DN.Risk)
  S.KM<-cumprod(1-DN.Risk)
  S.NA<-exp(-chr.NA)
  result<-list(S.NA=S.NA,S.KM=S.KM)
  return(result)
}



LS.int<-function(fatx,x){
  temp<-fatx[1:(length(fatx)-1)]
  d.x<-diff(x)
  cum.f<-x[1]+cumsum(d.x*temp)
  return(c(x[1],cum.f))
}



FG.transform<-function(time,censor,S.cens,g.gamma=0.0,tail=0.0){
  n.censored<-sum(censor==1)
  if(n.censored<=3){
    g.time<-time
  }
  else{
    # Calculate Integral
    # First modify S.cens
    S.mod<-ifelse(S.cens>g.gamma,S.cens,g.gamma)
    integrand.1<-1/S.mod
    
    #plot(time,S.mod,type="s",lty=1)
    #abline(h=0.02)
    
    int.1<-LS.int(integrand.1,time)
    
    num.theta<-int.1-time
    den.theta<-(time/S.mod)-int.1
    
    temp<-num.theta[censor==1]/den.theta[censor==1]
    theta<-min(temp)
    
    temp2<-ifelse(censor==1,(1+theta)*int.1,(1+theta)*int.1-(theta*time/S.mod))
    g.time<-ifelse(S.cens>=tail,temp2,time)
  }
  return(g.time)
}

get.FG<-function(df){
  df$id<-c(1:dim(df)[1])
  df.new<-df
  df1<-df[order(df$time),]
  id1<-df1$id
  tt<-df1$time
  cc<-1-df1$event.time
  temp<-NA.CHR(X=tt,D=cc)
  #plot(tt,temp$S.KM,type="s")
  g.time<-FG.transform(time=tt,censor=cc,S.cens=temp$S.NA,g.gamma=0.1,tail=0.1)
  id.match<-match(df$id,id1)
  df.new$Y.FG<-g.time[id.match]
  df.new$event.YFG<-rep(1,length(g.time))
  return(df.new)
}
