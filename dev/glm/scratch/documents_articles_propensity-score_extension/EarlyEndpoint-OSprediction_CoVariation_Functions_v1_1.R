
response.prediction.OS<-function(X.model,Y.name="HR OS",data.response,X.new,pi.level=0.70,details=TRUE,show.r2=TRUE,xlab="Delta",ylab="Y"){
  
  resp.model<-as.formula(paste("Y ~ ", paste(X.model, collapse= "+")))
  resp.fit<-lm(resp.model,data=data.response)
  
  Y.pred<-predict(resp.fit, X.new, se.fit = TRUE,interval="prediction",level=pi.level)
  
  if(details){
    cat("Summary of Model Fit","\n")
    print(summary(resp.fit))
  }
  
  data.fit<-data.frame(delta=X.new[,1],hr.delta=Y.pred$fit)
  names(data.fit)<-c("Delta",c(Y.name),"Lower PI","Upper PI")
  
  if(details){
    Y.hat<-data.fit[,2]
    Y.low<-data.fit[,3]
    Y.up<-data.fit[,4]
    
    ymax<-max(c(Y.hat,Y.up,data.response$Y))
    ymin<-min(c(Y.low,Y.hat,data.response$Y))
    
    plot(data.response$d1,data.response$Y,cex=0.25,ylim=c(ymin,ymax),xlab=xlab,ylab=ylab)     # Co-variation
    lines(lowess(data.response$d1,data.response$Y),lwd=3,col="orange")
    
    lines(d0,Y.hat,lwd=3,col="blue")
    lines(d0,Y.low,lwd=3,col="blue",lty=3)
    lines(d0,Y.up,lwd=3,col="blue",lty=3)
    
    if(show.r2){
      r.sq<-summary(resp.fit)$r.squared
      stat.toshow<-vector('expression',1)
      stat.toshow[1]=substitute(expression(italic(R^2)== MYVALUE), 
                                list(MYVALUE = format(r.sq,dig=3)))[2]
      legend("top", legend = stat.toshow, bty = 'n')
    }
  }
  return(list(resp.fit=resp.fit,data.fit=data.fit))
}



#################################################################
# Observed vs Predictions for "cross-validation" (eg, bootstrap)
#################################################################

cv.prediction.OS<-function(data.obs,data.predict.cr,data.predict.orr,
                           data.predict.pfs3,data.predict.pfs6,Xmodel=X3.model){
  
  data.treat<-subset(data.obs,Delta==1)    
  data.control<-subset(data.obs,Delta==0)    
  
  d1.cr<-mean(data.treat$cr)-mean(data.control$cr)
  d1.orr<-mean(data.treat$orr)-mean(data.control$orr)
  
  # PFS
  
  KM.S1<-survfit(Surv(tpfs,event_pfs)~1,data=data.treat)
  KM.S0<-survfit(Surv(tpfs,event_pfs)~1,data=data.control)
  
  km1.3<-summary(KM.S1,c(3))$surv
  km0.3<-summary(KM.S0,c(3))$surv
  d1.pfs3<-km1.3-km0.3
  
  km1.6<-summary(KM.S1,c(6))$surv
  km0.6<-summary(KM.S0,c(6))$surv
  d1.pfs6<-km1.6-km0.6
  
  
  # CR prediction
  
  X.new<-data.frame(d1=d1.cr,d2=d1.cr^2,d3=d1.cr^3,d4=d1.cr^4)
  data.cr1<-response.prediction.OS(X.model=Xmodel,data.response=data.predict.cr,
                                   X.new=X.new,details=FALSE)$data.fit
  
  # ORR prediction
  
  X.new<-data.frame(d1=d1.orr,d2=d1.orr^2,d3=d1.orr^3,d4=d1.orr^4)
  data.orr1<-response.prediction.OS(X.model=Xmodel,data.response=data.predict.orr,
                                    X.new=X.new,details=FALSE)$data.fit
  
  
  # PFS3 prediction
  
  X.new<-data.frame(d1=d1.pfs3,d2=d1.pfs3^2,d3=d1.pfs3^3,d4=d1.pfs3^4)
  
  data.pfs31<-response.prediction.OS(X.model=Xmodel,data.response=data.predict.pfs3,
                                     X.new=X.new,details=FALSE)$data.fit
  
  # PFS6 prediction
  
  X.new<-data.frame(d1=d1.pfs6,d2=d1.pfs6^2,d3=d1.pfs6^3,d4=d1.pfs6^4)
  
  data.pfs61<-response.prediction.OS(X.model=Xmodel,data.response=data.predict.pfs6,
                                     X.new=X.new,details=FALSE)$data.fit
  
  
  data.obs.predictions<-rbind(data.cr1,data.orr1,data.pfs31,data.pfs61)
  rownames(data.obs.predictions)<-c("CR","ORR","PFS@3","PFS@6")
  
  # Observed OS
  coxfit<-summary(coxph(Surv(tos,event_os)~Delta,data=data.obs))
  #coxfit.adj<-summary(coxph(Surv(tpfs,event_pfs)~Delta+age+male+
  #stage4+ecog0+bulky+ref.LP+LP2+LPge3+ipi35,data=data.obs))
  
  fit.hr<-coxfit$conf.int[c(1,3,4)]
  #fit.hr.adj<-coxfit.adj$conf.int[1,c(1,3,4)]
  
  data.hr<-c(fit.hr)
  names(data.hr)<-c("HR","95% Lower CI","95% Upper CI")
  
  #data.hr<-rbind(fit.hr,fit.hr.adj)
  #colnames(data.hr)<-c("HR","95% Lower CI","95% Upper CI")
  #rownames(data.hr)<-c("UnAdjusted","Adjusted")
  
  bias.unadj<-c(data.obs.predictions[,2])-c(fit.hr[1])
  #bias.adj<-c(data.obs.predictions[,2])-c(fit.hr.adj[1])
  
  cover.unadj<-ifelse(data.obs.predictions[,3]<=fit.hr[1] & data.obs.predictions[,4]>=fit.hr[1],1,0)
  #cover.adj<-ifelse(data.obs.predictions[,3]<=fit.hr.adj[1] & data.obs.predictions[,4]>=fit.hr.adj[1],
  #1,0)
  
  bias.adj<-NULL; cover.adj<-NULL
  
  return(list(data.obs.predictions=data.obs.predictions,hr.obs=data.hr,
              bias.unadj=bias.unadj,bias.adj=bias.adj,cover.unadj=cover.unadj,
              cover.adj=cover.adj))
}
