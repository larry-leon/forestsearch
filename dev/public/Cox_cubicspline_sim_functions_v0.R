
require(randomizr)


SG_forest <- function(df,arrow_text=c("favors Treatment","Control"),E.name=c("Treat"),C.name=c("Control"),outcome.name,event.name,treat.name,footnote_text=NULL){
  title_text <- NULL
  # ITT 
  res_itt <- SG_HRtable2(dfa=df,df.OOB=NULL,fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,
                         sg_name="ITT",E.name=E.name,C.name=C.name)
  # Add separator header 
  zz <- res_itt  
  zz$Subgroup <- "                 "
  zz[,c(E.name,C.name)] <- ""
  zz[,c("est","low","hi","se")] <- NA
  
  dfsg <- subset(df,cps_1==0)
  sg_name <- c("CPS <1")
  res_sg1Q <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_1==1)
  sg_name <- c("CPS >=1")
  res_sg1B <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_5==0)
  sg_name <- c("CPS <5")
  res_sg2Q <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_5==1)
  sg_name <- c("CPS >=5")
  res_sg2B <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_10==0)
  sg_name <- c("CPS <10")
  res_sg3Q <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_10==1)
  sg_name <- c("CPS >=10")
  res_sg3B <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                          fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_15==1)
  sg_name <- c("CPS 1-<5")
  res_sg4 <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                         fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_510==1)
  sg_name <- c("CPS 5-<10")
  res_sg5 <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                         fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dfsg <- subset(df,cps_19==1)
  sg_name <- c("CPS 1-<10")
  res_sg6 <- SG_HRtable2(dfa=dfsg,df.OOB=NULL,
                         fa_SG=NULL,ntreat_fa=NA,ncontrol_fa=NA,outcome.name=outcome.name,event.name=event.name,treat.name=treat.name,sg_name=sg_name,E.name=E.name,C.name=C.name)
  
  dt <- rbind(res_itt,zz,res_sg1Q,res_sg1B,zz,res_sg2Q,res_sg2B,zz,res_sg3Q, res_sg3B, zz, res_sg4,res_sg5,res_sg6)
  
  sg_colors <- c("white")
  
  tm <- forest_theme(core=list(fg_params=list(hjust = 1, x = 0.9),
                               bg_params=list(fill = sg_colors)),
                     colhead=list(fg_params=list(hjust=0.5, x=0.5)),
                     footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"))
  
  dt$` ` <- paste(rep(" ", 20), collapse = " ")
  
  # Create a confidence interval column to display
  dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
                             sprintf("%.2f (%.2f to %.2f)", dt$est, dt$low, dt$hi))
  
  p <- forest(dt[,c(1:3, 8:9)],
              title=title_text,
              est = dt$est,
              lower = dt$low, 
              upper = dt$hi,
              sizes = dt$se,
              ci_column = 4,
              ref_line = 1,
              arrow_lab = arrow_text,
              xlim = c(0.25, 2.0),
              ticks_at = c(0.65, 1.0, 1.3),
              footnote = footnote_text, theme=tm)
  
  return(list(dt=dt,p=p))
}


get_dfCPS<-function(df){
df_cps<-within(df,{
  cps_10 <- ifelse(cps10=="Positive",1,0)
  cps_5 <- ifelse(cps >=5 | cps10=="Positive",1,0)
  cps_5[is.na(cps)] <- NA
  cps_1 <- ifelse(cps1=="Positive",1,0)
  cps_1[c(cps1=="")] <- NA
  temp_flag <- c(cps10=="Negative" | (cps10=="" & cps<10))
  temp_flag[is.na(cps)] <- FALSE
  cps_19 <- ifelse(cps1=="Positive" & temp_flag,1,0)
  cps_15 <- ifelse(cps_1==1 & cps_5==0,1,0)
  cps_510 <- ifelse(cps_5==1 & cps_10==0,1,0)
})
return(df_cps)
}


get_dgm <- function(df,knot=5,kappa=10,hrs=c(0.75,0.75,0.75),details=TRUE,masked=TRUE){

if(masked){
    df <- subset(df,mflag==1)  
    df <- within(df,{
      event <- event-en
      tte <- tte-yn
    })  
  }  
  
  # Spline at 1 knot
  dfa2 <- within(df,{
    z.treat <- z*treat
    z.k <- (z-knot)*ifelse(z>knot,1,0)
    z.k.treat <- z.k*treat
  })

fit.weibk <- survreg(Surv(tte,event) ~ treat+z+z.treat+z.k+z.k.treat, 
                     dist='weibull', data=dfa2)
  
# Censoring model independent of covariates

fitC.weib <- survreg(Surv(tte,1-event) ~ 1, dist='weibull', data=dfa2)
sigC <- c(fitC.weib$scale)
muC <- c(coef(fitC.weib)[1])
  
  # Increase strength of treatment terms
  fit.sims <- fit.weibk
  
  modify_coefficients <- TRUE
  
  # knot is point of decrease
  # kappa is "strongest decrease"
  
  loghr.0 <- log(hrs[1])
  loghr.knot <- log(hrs[2])
  loghr.kappa <- log(hrs[3])
  
  # Solve for gamma to satisfy
  # b0[1]+10*b0[3]+5*b0[5]=loghr.10
  # b0[1]+5*b0[3]=loghr.5
  # b0[1] = loghr.0
  
  sig <- c(fit.sims$scale)
  mu <- c(coef(fit.sims)[1])
  gamma <- c(coef(fit.sims)[c(-1)])
  
  if(details) print(c(-gamma/sig))
  
  if(modify_coefficients){
    # At the change-point
    a0 <- c(-gamma)/sig
    a0[1] <- loghr.0
    a0[3] <- (loghr.knot-loghr.0)/knot
    a0[5] <- (loghr.kappa-a0[1]-kappa*a0[3])/knot
    # back-solve gamma
    gamma <- -a0*sig
  }
  
  # Weibull log-scale
  b.true <- gamma
  # hazard rate scale parameters
  b0 <- c(-gamma)/sig
  
  if(details) print(b0)
  
  # predictions log(hr)
  dfp <- dfa2[order(dfa2$z,decreasing=FALSE),]

return(list(df_super=dfp,b.true=b.true,b0=b0,mu=mu,sig=sig,muC=muC,sigC=sigC))
}



# Include single covariate X (Eg; "ecogbl", "prior_line12")

draw_sim <- function(dgm,ss=1,details=FALSE,Ndraw=nrow(dgm$df_super),xname=c("ecogbl"),bx=0){

df_super <- dgm$df_super
b.true <- dgm$b.true
mu <- dgm$mu
sig <- dgm$sig
muC <- dgm$muC 
sigC <- dgm$sigC

set.seed(8316951+ss*1000)

  if(Ndraw!=nrow(df_super)){
  dfNew <- df_super  
  id_sample <- sample(c(1:nrow(dfNew)),size=Ndraw,replace=TRUE)
  df_super <- dfNew[id_sample,]
  }

N <- nrow(df_super)

# "treat" is a placeholder  
zmat <- as.matrix(df_super[,c("treat","z","z.treat","z.k","z.k.treat")])
x <- df_super[,c(xname)]

# Initiate as population with covariates 
# Outcomes and treatment assignment appended below
dfsim <- df_super 

strata <- df_super$stratum

# set treatment to 1
zmat.1 <- zmat
zmat.1[,"treat"] <- 1.0
zmat.1[,"z.treat"] <- zmat.1[,"z"]
zmat.1[,"z.k"] <- zmat.1[,"z.k"]
zmat.1[,"z.k.treat"] <- zmat.1[,"z.k"]

# set treatment to 0
zmat.0 <- zmat
zmat.0[,"treat"] <- 0.0
zmat.0[,"z.treat"] <- 0.0
zmat.0[,"z.k"] <- zmat.0[,"z.k"]
zmat.0[,"z.k.treat"] <- 0.0

# error term (log-scale [linear predictor scale])
  epsilon <- log(rexp(N))
  # linear predictor (setting treat=1)
  z.1 <- zmat.1[,"z"]
  eta1 <- mu + c(zmat.1%*%b.true)+x*bx
  # Excluding baseline hazard term
  phi1 <- -c(zmat.1%*%b.true+x*bx)/sig
  
  log.Y1 <- eta1 + sig*epsilon
  #strata1 <- strata

  z.0 <- zmat.0[,"z"]
  eta0 <- mu + c(zmat.0%*%b.true)+x*bx
  log.Y0 <- eta0 + sig*epsilon
  #strata0 <- strata

  # Excluding baseline hazard term
  phi0 <- -c(zmat.0%*%b.true+x*bx)/sig
  
  # Potential outcome difference
  loghr.po <- phi1-phi0
  
  # Randomize per stratar

  blocks <- strata
  
  # Randomize treatment
  Zr <- block_ra(blocks=blocks)
 
  log.Yr <- Zr*log.Y1 + (1-Zr)*log.Y0
  
  Yr <- exp(log.Yr)
  
  # At this point censoring is completely independent
  epsilonC <- log(rexp(N))
  log.YC <- muC + sigC*epsilonC
  
  log.YrC <- Zr*log.YC+(1-Zr)*log.YC 

  YrC <- exp(log.YrC)

    # Observed censored outcomes

  dfsim$event.sim <- ifelse(Yr<=YrC,1,0)
  dfsim$y.sim <- pmin(Yr,YrC)
  dfsim$treat.sim <- Zr
  dfsim$strata.sim <- strata
  dfsim$x <- x
    
  dfsim$loghr.po <- loghr.po
  dfsim$log.Y1 <- log.Y1
  dfsim$log.Y0 <- log.Y0
  
  if(details & ss <= 10) cat("% censored (first 10)=",mean(1-dfsim$event.sim),"\n")
 
  return(dfsim)
}

