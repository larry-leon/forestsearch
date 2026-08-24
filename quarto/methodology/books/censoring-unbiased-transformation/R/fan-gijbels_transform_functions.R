
# Lebesque-Stieljes integral
LS.int<-function(fatx,x){
  temp<-fatx[1:(length(fatx)-1)]
  d.x<-diff(x)
  cum.f<-x[1]+cumsum(d.x*temp)
  return(c(x[1],cum.f))
}



FG.transform<-function(time,censor,S.cens,g.gamma=0.0,tail=0.0){
  n.censored<-sum(censor==1)
  if(n.censored<=3){
    warning("Number of censored observations less than 3, no transformation applied (return original times)")
    g.time<-time
  }
  else{
    # Calculate Integral
    # First modify S.cens
    S.mod<-ifelse(S.cens>g.gamma,S.cens,g.gamma)
    integrand.1<-1/S.mod
    
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
# K-M of censoring distribution (Events are now censorings)
temp<-NA.CHR.Weighted(time=tt,Delta=cc)
g.time<-FG.transform(time=tt,censor=cc,S.cens=temp$S.NA,g.gamma=0.1,tail=0.1)
id.match<-match(df$id,id1)
df.new$Y.FG<-g.time[id.match]
df.new$event.YFG<-rep(1,length(g.time))
return(df.new)
}
  
  
  





