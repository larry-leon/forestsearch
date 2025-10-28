# AFT dgm based on gbsg (survival package)
# Harm subgroup as found via forest-search procedure
# Simulate data in same format as original data

# Add treatment randomization (RR) and rand stratification


# cut.factor and cut.num previously in sim_aft_gbsg function files
cut.factor<-function(x,cutpoint=median(x),direction="LTE",labels=c("1","2")){
  if(direction=="GT"){
    cut.x<-factor(ifelse(x>cutpoint,1,0),levels=0:1,labels=labels)
  }
  if(direction=="LTE"){
    cut.x<-factor(ifelse(x<=cutpoint,1,0),levels=0:1,labels=labels)
  }
  return(cut.x)
}

cut.num<-function(x,cutpoint=median(x),direction="LTE"){
  if(direction=="GT"){
    cut.x<-ifelse(x>cutpoint,1,0)
  }
  if(direction=="LTE"){
    cut.x<-ifelse(x<=cutpoint,1,0)
  }
  return(cut.x)
}




dgm_aftm4_gbsg<-function(model="alt",k.treat=1,k.inter=1,details=TRUE,cens.type="weibull",parms_torand=FALSE){

  if (!("survival" %in% utils::installed.packages())) stop("gbsg dataset requires the survival library.")
  if(!is.element(cens.type,c("weibull","uniform"))) stop("Censoring type must be weibull or uniform")

  if(!(model %in% c("alt","null"))) stop("model must be either alt or null","\n")

  dfa<-gbsg

  dfa<-within(dfa,{
    y<-rfstime/30.4375
    id<-c(1:nrow(dfa))
    event<-ifelse(status==1,1,0)
    treat<-hormon
    z1<-ifelse(er<=quantile(er,c(z1_frac)),1,0)
    z2<-cut.num(age)
    z4<-cut.num(pgr)
    #z3<-ifelse(pgr<=32,1,0)
    #z3<-ifelse(pgr<=quantile(pgr,c(z3_frac)),1,0)
    z3<-ifelse(meno==0,1,0)
    z5<-cut.num(nodes)
    #z5<-ifelse(grade==3,1,0)
    zh<-treat*z1*z3
    flag.harm<-ifelse(z1==1 & z3==1,1,0)
    # Observed factors
    # For now, assume observe exact cuts
    v1<-as.factor(z1)
    v2<-as.factor(z2)
    v3<-as.factor(z3)
    v4<-as.factor(z4)
    v5<-as.factor(z5)
    # v6 and v7 are "noise"
    #v6<-as.factor(grade)
    # random binary facto
    #v7<-rbinom(nrow(dfa),size=1,prob=0.5)
    v6<-cut.factor(size)
    grade3<-ifelse(grade==3,1,0)
    v7<-as.factor(grade3)
  })

  # output subgroup identities
  # For forestsearch
  if(model!="null"){
    fs.harm.true<-c("v1.1","v3.1")
    # For virtual twins binary factors (e.g., v1=1) are identified as v1>=0.5
    # and v1=0 "v1< 0.5"
    vt.harm.true<-c("v1>=0.5","v3>=0.5")
    # For grf
    grf.harm.true<-c("v1=1","v3=1")
  }

  if(model=="null"){
    fs.harm.true<-vt.harm.true<-grf.harm.true<-NULL
    # Reset dfa$flag.harm=0
    dfa$flag.harm<-0
  }

  # Underlying covariates
  covs.true<-c("treat","z1","z2","z3","z4","z5")
  # If H subgroup, the model covariates will
  # have additional interaction
  # This indicates location in z.true
  loc.inter<-length(covs.true)+1

  # These are the factors available to the analyst
  covs.analysis<-c("v1","v2","v3","v4","v5","v6","v7")

  # Super population:  In sim() module will allow for sampling from df.super
  # if want (1) Different sample size and/or (2) Random covariates
  df.super<-dfa[,c("y","id","event","treat","flag.harm",covs.analysis)]

  # Note: hr.ratio.true is the causal hazard ratio (appended to df.sim below);
  # This is on patient level so that theta[hat(H)] can be calculated for any hat(H); see below.

  if(model!="null"){
    z.true<-as.matrix(dfa[,c(covs.true,"zh")])
    # Potential outcome versions
    z.true.1<-z.true
    z.true.1[,1]<-1 # Set treatment to treated (counterfactual)
    z.true.1[,loc.inter]<-z.true[,"z1"]*z.true[,"z3"]

    z.true.0<-z.true
    z.true.0[,1]<-0 # Set treatment to control
    z.true.0[,loc.inter]<-0
  }

  if(model=="null"){
    z.true<-as.matrix(dfa[,c(covs.true)])
    # Potential outcome versions
    z.true.1<-z.true
    z.true.1[,1]<-1

    z.true.0<-z.true
    z.true.0[,1]<-0
  }
  # Weibull model
  y<-dfa$y
  event<-dfa$event
  fit0<-survreg(Surv(y,event) ~ z.true, dist='weibull')

  names(fit0$coefficients)<-c("(Intercept)",colnames(z.true))

  #################
  # weibull scaling
  #################

  sig<-c(fit0$scale)
  mu<-c(coef(fit0)[1])
  gamma<-c(coef(fit0)[c(-1)])

  # Modify to rand version?
  if(parms_torand){
    # Set coefficients for treatment, z1,z2,z4, and z5 to
    # randomized analysis with meno==0
    df.rand<-subset(dfa,meno==0)
    fit.rand<-survreg(Surv(y,event) ~ treat+z1+z2+z3+z4+z5, dist='weibull',data=df.rand)
    gamma.rand<-c(coef(fit.rand)[c(-1)])
    # Only modify confounders
    gamma[c("z1","z2","z4","z5")]<-gamma.rand[c("z1","z2","z4","z5")]
  }


  if(model=="null"){
    gamma[1]<-k.treat*gamma[1] # Increase treatment strength
    # Modify impact of z3
    gamma[4]<-k.z3*gamma[4]

  }
  if(model!="null"){
    # Modify impact of z3
    gamma[4]<-k.z3*gamma[4]
    # Increase treatment and interaction effects
    gamma[1]<-k.treat*gamma[1]
    gamma[loc.inter]<-k.inter*gamma[loc.inter]
  }

  b.true<-gamma # Weibull log-scale
  b0<-c(-gamma)/sig

  lin.conf<-z.true%*%b.true
  lin1.conf<-z.true.1%*%b.true
  lin0.conf<-z.true.0%*%b.true

  # Potential outcome hazard "link" hr(1)/hr(0) for each subject
  # baseline cancels
  dfa$hlin.conf.1<-exp(z.true.1%*%b0)
  dfa$hlin.conf.0<-exp(z.true.0%*%b0)

  # Generate and return randomized data set
  # for "super population"
  # Generate large popln
  dfsamp<-dfa
  # potential outcome link (Weibull)
  dfsamp$lin1.conf<-lin1.conf
  dfsamp$lin0.conf<-lin0.conf
  # Randomize
  set.seed(8316951)
  nbig<-5000
  na<-round(nbig/2,0)
  nb<-nbig-na
  id_mix<-sample(c(1:nrow(dfsamp)),size=nbig,replace=TRUE)
  # Experimental
  dfsamp<-dfsamp[id_mix,]
  dfas<-dfsamp[c(1:na),]
  dfas$treat<-1
  # Set link to potential outcome for treat
  dfas$lin.conf.true<-dfas$lin1.conf
  # Control
  dfbs<-dfsamp[c((na+1):nbig),]
  dfbs$treat<-0
  dfbs$lin.conf.true<-dfbs$lin0.conf
  dfbig<-as.data.frame(rbind(dfas,dfbs))
  if(model=="alt"){
    # "Empirical versions"
    # Generate Weibull outcomes
    lin.conf<-dfbig$lin.conf.true
    epsilon<-log(rexp(nbig)) # Extreme value distribution
    logTs<-mu+sig*epsilon+lin.conf
    Ts<-exp(logTs)
    dfbig$Ts<-Ts
    dfbig$es<-1
    hr.H.true<-exp(coxph(Surv(Ts,es)~treat,data=subset(dfbig,flag.harm==1))$coefficients)
    hr.Hc.true<-exp(coxph(Surv(Ts,es)~treat,data=subset(dfbig,flag.harm==0))$coefficients)
    hr.causal<-exp(coxph(Surv(Ts,es)~treat,data=dfbig)$coefficients)
    size.H<-sum(ifelse(z.true[,"z1"]==1 & z.true[,"z3"]==1,1,0))
    size.Hc<-sum(ifelse(z.true[,"z1"]==0 | z.true[,"z3"]==0,1,0))
  }
  if(model=="null"){
    lin.conf<-dfbig$lin.conf.true
    epsilon<-log(rexp(nbig)) # Extreme value distribution
    logTs<-mu+sig*epsilon+lin.conf
    Ts<-exp(logTs)
    dfbig$Ts<-Ts
    dfbig$es<-1
    hr.H.true<-NA
    hr.causal<-exp(coxph(Surv(Ts,es)~treat,data=dfbig)$coefficients)
    hr.Hc.true<-hr.causal
  }

  if(details){
    cat("Super-population empirical harm and non-harm hazard ratios=",c(hr.H.true,hr.Hc.true),"\n")
    cat("Causal HR (empirical ITT)=",c(hr.causal),"\n")
  }

  df.super_rand<-dfbig

  # Censoring
  # Fit log-linear censoring models
  # If censoring is uniform, then cens.parms is null
  linC.conf<-NULL
  sigC<-NULL
  muC<-NULL
  # Otherwise:
  if(cens.type!="uniform"){
    zC.true<-as.matrix(df.super_rand[,c(covs.true)])
    if(cens.type=="weibull"){
      y<-df.super_rand$y
      event<-df.super_rand$event
      fit0C<-survreg(Surv(y,1-event) ~ zC.true, dist='weibull')
      sigC<-c(fit0C$scale)
      gammaC<-c(coef(fit0C)[c(-1)])
      muC<-c(coef(fit0C)[1])
      b0C<-c(-gammaC)/sigC
      # Back to weibull scale
      bC.true<--b0C*sig
      linC.conf<-zC.true%*%bC.true

      # Potential outcome versions
      zC.true.1<-zC.true
      zC.true.1[,1]<-1 # Set treatment to treated (counterfactual)
      zC.true.0<-zC.true
      zC.true.0[,1]<-0 # Set treatment to control

      linC1.conf<-zC.true.1%*%bC.true
      linC0.conf<-zC.true.0%*%bC.true
    }
  }

  df.super_rand<-within(df.super_rand,{
    hlin.ratio<-hlin.conf.1/hlin.conf.0
    linC1.conf<-linC1.conf
    linC0.conf<-linC0.conf
    h1.potential<-hlin.conf.1
    h0.potential<-hlin.conf.0
  })


  return(list(df.super_rand=df.super_rand,sigC=sigC,hr.H.true=hr.H.true,hr.Hc.true=hr.Hc.true,cens.type=cens.type,mu=mu,hr.causal=hr.causal,b.true=b.true,b_hr.true=b0,
              sig=sig,muC=muC,covs.analysis=covs.analysis,fs.harm.true=fs.harm.true,vt.harm.true=vt.harm.true,grf.harm.true=grf.harm.true,model=model))
}


# If n is specified then will sample from df.super

sim_aftm4_gbsg<-function(dgm,n=NULL,rand.ratio=1,simid=1,minc=NULL,maxc=NULL,maxFollow=Inf,muC.adj=0.0,drawtreat=TRUE){
  set.seed(8316951+1000*simid)
  cens.type<-dgm$cens.type

  if(cens.type=="uniform" & is.null(minc) & is.null(maxc)) stop("For uniform censoring, minc and maxc need to be specified")

  # "Super randomized" version
  # already randomized
  df.super<-dgm$df.super_rand
  # Simulated data will be appended
  if(is.null(n)) df.sim<-df.super
  if(!is.null(n)){
    # Randomly sample n
    n0<-round(n*(1/(1+rand.ratio)))
    n1<-n-n0
    ###############################################
    # Original, sampling from gbsg treatment arms
    # "retain confounding per study"
    # Only used if sampling from origina gbsg dataset
    ###############################################
    # Sample n1 from experimental of df.super
    if(!drawtreat){
      df1.super<-subset(df.super,treat==1)
      df0.super<-subset(df.super,treat==0)
      id1.draw<-sample(1:nrow(df1.super),size=n1,replace=TRUE)
      df1.sim<-df1.super[id1.draw,]
      id0.draw<-sample(1:nrow(df0.super),size=n0,replace=TRUE)
      df0.sim<-df0.super[id0.draw,]
      df.sim<-rbind(df0.sim,df1.sim)
    }
    # Randomly draw
    if(drawtreat){
      id_sample<-sample(c(1:nrow(df.super)),size=n,replace=TRUE)
      df.supernew<-df.super[id_sample,]
      # First 1:n1 as treatment
      #df1.sim<-df.super[c(1:n1),c("flag.harm","v1","v2","v3","v4","v5","v6","v7","lin1.conf","linC1.conf","lin0.conf","linC0.conf","hr.potential")]
      df1.sim<-df.supernew[c(1:n1),]
      df1.sim$treat<-1
      df1.sim$lin.conf.true<-df1.sim$lin1.conf
      df1.sim$linC.conf<-df1.sim$linC1.conf

      #df0.sim<-df.super[c((n1+1):n),c("flag.harm","v1","v2","v3","v4","v5","v6","v7","lin1.conf","linC1.conf","lin0.conf","linC0.conf","hr.potential")]
      df0.sim<-df.supernew[c((n1+1):n),]
      df0.sim$treat<-0
      df0.sim$lin.conf.true<-df0.sim$lin0.conf
      df0.sim$linC.conf<-df0.sim$linC0.conf
      df.sim<-as.data.frame(rbind(df0.sim,df1.sim))
    }
  } # draw n
  n<-nrow(df.sim)
  mu<-dgm$mu
  sig<-dgm$sig
  lin.conf<-df.sim$lin.conf.true
  # Extreme value error term
  epsilon<-log(rexp(n)) # Extreme value distribution
  logT.sim<-mu+sig*epsilon+lin.conf
  T.sim<-exp(logT.sim)

  muC<-dgm$muC+muC.adj
  sigC<-dgm$sigC
  linC.conf<-df.sim$linC.conf

  # Weibull
  if(cens.type=="weibull"){
    epsilonC<-log(rexp(n)) # Extreme value distribution
    logC.sim<-muC+sigC*epsilonC+linC.conf
    C.sim<-c(exp(logC.sim))
  }
  if(cens.type=="uniform") C.sim<-runif(n,min=minc,max=maxc)
  # Administratic censoring
  C.sim<-pmin(C.sim,maxFollow)
  event.sim<-ifelse(T.sim<=C.sim,1,0)
  y.sim<-pmin(T.sim,C.sim)
  # ID's need to be re-set
  df.sim<-within(df.sim,{
    y.sim<-y.sim
    event.sim<-event.sim
    t.sim<-T.sim
    id<-c(1:nrow(df.sim))
  })
  return(as.data.frame(df.sim))
}


# Root for hrH.obs=0
hrH.obj4<-function(kval,hrH.target,model,k.treat,cens.type){
  val<-dgm_aftm4_gbsg(k.inter=kval,model=model,k.treat=k.treat,details=FALSE,cens.type=cens.type)$hr.H.true
  return(c(val-hrH.target))
}

propH.obj4<-function(zq,pH.target){
  #val<-with(gbsg,mean(pgr<=quantile(pgr,c(z3_frac),1,0) & er<=quantile(er,zq)))
  val<-with(gbsg,mean(meno==0 & er<=quantile(er,zq)))
  return(c(val-pH.target))
}


get.dgm4<-function(mod.harm,N,k.treat,hrH.target,cens.type,out.loc=NULL,file.index=NULL,model.index=NULL,details=FALSE,parms_torand=FALSE,sol_tol=10^-8){
  k.inter<-NULL
  if(mod.harm!="null"){
    k.inter<-uniroot(hrH.obj4,c(-100,100),tol=sol_tol,model=mod.harm,k.treat=k.treat,
                     cens.type=cens.type,hrH.target=hrH.target)$root
  }
  dgm<-dgm_aftm4_gbsg(model=mod.harm,k.treat=k.treat,k.inter=k.inter,details=details,cens.type=cens.type,parms_torand=parms_torand)

  if(!is.null(out.loc)){
    out.file<-paste0(out.loc,model.index,sep="_")
    nsize<-paste0("N=",N)
    out.file<-paste0(out.file,nsize,sep="_")
    out.file<-paste0(out.file,mod.harm,sep="_")
    treat.id<-paste0("ktreat=",k.treat)
    out.file<-paste0(out.file,treat.id,sep="_")
  }

  # For bootstrap this is duplicate
  #if(mod.harm=="alt"){
  #hrH<-round(dgm$hr.H.true,2)
  #inter.id<-paste0("hrH=",hrH)
  #if(!is.null(out.loc))  out.file<-paste0(out.file,inter.id,sep="_")
  #}


  if(!is.null(out.loc)){
    out.file<-paste0(out.file,file.index)
    out.file<-paste0(out.file,".Rdata")
  }

  if(is.null(out.loc)) out.file <- NULL

  return(list(dgm=dgm,out.file=out.file))
}


get.dgm4.OC<-function(mod.harm,N,k.treat,hrH.target,cens.type,out.loc=NULL,file.index=NULL,model.index=NULL,details=FALSE,parms_torand=FALSE,sol_tol=10^-8){
  k.inter<-NULL
  if(mod.harm!="null"){
    k.inter<-uniroot(hrH.obj4,c(-100,100),tol=sol_tol,model=mod.harm,k.treat=k.treat,
                     cens.type=cens.type,hrH.target=hrH.target)$root
  }
  dgm<-dgm_aftm4_gbsg(model=mod.harm,k.treat=k.treat,k.inter=k.inter,details=details,cens.type=cens.type,parms_torand=parms_torand)

  if(!is.null(out.loc)){
    out.file<-paste0(out.loc,model.index,sep="_")
    nsize<-paste0("N=",N)
    out.file<-paste0(out.file,nsize,sep="_")
    out.file<-paste0(out.file,mod.harm,sep="_")
    treat.id<-paste0("ktreat=",k.treat)
    out.file<-paste0(out.file,treat.id,sep="_")
  }

  # For bootstrap this is duplicate
  # NOT for OC
  if(mod.harm=="alt"){
    hrH<-round(dgm$hr.H.true,2)
    inter.id<-paste0("hrH=",hrH)
    if(!is.null(out.loc))  out.file<-paste0(out.file,inter.id,sep="_")
  }


  if(!is.null(out.loc)){
    out.file<-paste0(out.file,file.index)
    out.file<-paste0(out.file,".Rdata")
  }

  if(is.null(out.loc)) out.file <- NULL

  return(list(dgm=dgm,out.file=out.file))
}





library(survival)
summary(gbsg)
library(cubature)

N <- 700
maxFollow<-84
cens.type<-"weibull"

outcome.name<-c("y.sim")
event.name<-c("event.sim")
id.name<-c("id")
treat.name<-c("treat")

cox.formula.sim<-as.formula(paste("Surv(y.sim,event.sim)~treat"))
cox.formula.adj.sim<-as.formula(paste("Surv(y.sim,event.sim)~treat+v1+v2+v3+v4+v5"))

# m1 -censoring adjustment
muC.adj<-log(1.5)
z1_frac <- 0.25

  k.z3 <- 1.0
  k.treat <- 0.9
  pH_super <- 0.125 # non-NULL re-defines z1_frac


  if(is.null(pH_super)){
    #pH_check<-with(gbsg,mean(pgr<=quantile(pgr,c(z3_frac),1,0) & er<=quantile(er,z1_frac)))
    pH_check<-with(gbsg,mean(meno==0 & er<=quantile(er,z1_frac)))
    cat("Underlying pH_super",c(pH_check),"\n")
  }
  # pH_super specified
  # If pH_super then override  z1_frac and find z1_frac to yield pH_super

  if(!is.null(pH_super)){
    # Approximate Z1 quantile to yield pH proportion
    z1_q<-uniroot(propH.obj4,c(0,1),tol=0.0001,pH.target=pH_super)$root
    #pH_check<-with(gbsg,mean(pgr<=quantile(pgr,c(z3_frac),1,0) & er<=quantile(er,z1_q)))
    pH_check<-with(gbsg,mean(meno==0 & er<=quantile(er,z1_q)))
    cat("pH",c(pH_check),"\n")
    rel_error<-(pH_super-pH_check)/pH_super
    if(abs(rel_error)>=0.1) stop("pH_super approximation relative error exceeds 10%")
    z1_frac<-z1_q
    cat("Underlying pH_super",c(pH_check),"\n")
  }

  mod.harm <- "alt"
  hrH.target <- 2.0

  out.loc <- NULL

  this.dgm <- get.dgm4.OC(mod.harm=mod.harm,N=N,k.treat=k.treat,model.index=model.index,sol_tol=10^-8,
                        hrH.target=hrH.target,cens.type=cens.type,out.loc=out.loc,file.index=file.index,details=TRUE,parms_torand=FALSE)

  dgm <- this.dgm$dgm






