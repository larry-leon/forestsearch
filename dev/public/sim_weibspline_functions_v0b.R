# Based on sim_functions_v4 from "Study F/working"

get_dgm_stratified <- function(df,knot=5,kappa=10,log.hrs=log(c(0.75,0.75,0.75)),details=FALSE,strata_tte=NULL,masked=TRUE){

  if(masked){
    df <- subset(df,mflag==1)  
    df <- within(df,{
      event <- event-en
      tte <- tte-yn
    })  
  }  
  # Spline at 1 knot
  dfa2<-within(df,{
    z.treat <- z*treat
    z.k <- (z-knot)*ifelse(z>knot,1,0)
    z.k.treat <- z.k*treat
  })

  # strata 
  if(!is.null(strata_tte)){  
    aa <- paste("strata(",eval(strata_tte),")")
    bb <- c("Surv(tte,event) ~ treat + z + z.treat + z.k + z.k.treat +")
    weib.formula <- as.formula(paste(bb,aa))
  }
  
  if(is.null(strata_tte)){
    weib.formula <- as.formula("Surv(tte,event) ~ treat + z + z.treat + z.k + z.k.treat")
  }
  
  fit.weibk <- survreg(weib.formula, dist='weibull', data=dfa2)
  
  # Censoring model independent of covariates
  
  fitC.weib <- survreg(Surv(tte,1-event) ~ 1, dist='weibull', data=dfa2)
  tauC <- c(fitC.weib$scale)
  muC <- c(coef(fitC.weib)[1])
  
  # Relative log(hr) scale
  # With stratification the hr's will depend on the stratum (scale) 
  # To simplify we will take median ("tau_approx" below)
  
  mu <- c(coef(fit.weibk)[1])
  gamma <- c(coef(fit.weibk)[c(-1)])
  tau <- c(fit.weibk$scale)  
  
  # Return as strataO --> outcome stratification
  # As opposed to randomization stratification which 
  # will be included in sim() function below
  
  if(!is.null(strata_tte)){
    # Outcome stratificaton
    strataO <- dfa2[,c(strata_tte)]
    # Extract scales corresponding to strataO
    aa <- names(tau)
    tau_ids <- unlist(lapply(strataO,function(x){grep(x,aa)}))
    tau.strataO <- tau[tau_ids]
    if(length(tau.strataO)!=nrow(dfa2)) stop("strata_tte not uniquely identified via matching")
    tau.approx <- median(tau.strataO)
    dfa2$tau.strataO <- tau.strataO
  }
  
  if(is.null(strata_tte)){
    strataO <- "All"  
    tau.strataO <- tau 
    tau.approx <- tau
    dfa2$tau.strataO <- tau
    }
  
  # Re-define to satisfy log(hrs) pattern 
  
  loghr.0 <- log.hrs[1]
  loghr.knot <- log.hrs[2]
  loghr.kappa <- log.hrs[3]
  
  # At the change-point
  # Weibull hazard-ratio parameters
  b0 <- c(-gamma)/tau.approx
  # solve for b1,b3,b5
  b0[1] <- loghr.0
  b0[3] <- (loghr.knot-b0[1])/knot
  b0[5] <- (loghr.kappa-b0[1]-kappa*b0[3])/(kappa-knot)
  
  # Re-define gamma on AFT (log(T)) parameterization 

  gamma.true <- -b0*tau.approx 

  # predictions log(hr)
  
  #dfp <- dfa2[order(dfa2$z,decreasing=FALSE),]
  
  dfp <- data.table::setorder(dfa2,z)
  
  if(details){
  # Plot across range of z
  zz <- seq(0,max(dfp$z),by=1)
  # psi(z)
  loghr.zz <- b0[1]+b0[3]*zz+b0[5]*(zz-knot)*ifelse(zz>knot,1,0)
  plot(zz,loghr.zz,type="s",lty=1,xlab="z",ylab="psi(z)")
  abline(h=c(log.hrs))
  abline(v=c(0,knot,kappa),lwd=0.5,col="blue",lty=2)
  abline(h=0,lwd=0.25,col="red",lty=1)
      }
  
  return(list(df_super=dfp,gamma.true=gamma.true,mu=mu,tau=tau,muC=muC,tauC=tauC,strata_tte=strata_tte,tau.approx=tau.approx))
}

# Include single covariate X (Eg; "ecogbl", "prior_line12")
# strata_tte
# strata_rand

# Add hrz_crit (=log(1.2)) to set biomarker HR threshold that is "acceptable"
# For example biomarkers for which the true log(hrs) are < 1.2 allowing for at most 20% increase

draw_sim_stratified <- function(dgm,ss=1,details=FALSE,Ndraw=nrow(dgm$df_super),strata_rand=c("strataNew"),xname=c("ecogbl"),bx=0,checking=FALSE,
                                hrz_crit=log(1.2), return_df=TRUE){

df_super <- dgm$df_super

var_names <- c(strata_rand,xname)
if(all(var_names %in% names(df_super)) != TRUE) stop("strata_rand and xname variables not in dgm$df_super")   

  # Outcomes
  # Regression parameters on AFT scale (gamma)
  gamma.true <- dgm$gamma.true
  mu <- dgm$mu
  tau <- dgm$tau
  # Censoring
  muC <- dgm$muC 
  tauC <- dgm$tauC
  
  strata_tte <- dgm$strata_tte
  
  set.seed(8316951+ss*1000)
  
  # Sampling from df_super if Ndraw differs from size of df_super
  # Otherwise simulated dataset dfsim is initiated as df_super
  if(Ndraw!=nrow(df_super)){
    dfNew <- data.table::copy(df_super)  
    id_sample <- sample(c(1:nrow(dfNew)),size=Ndraw,replace=TRUE)
    df_super <- dfNew[id_sample,]
  }
  
  # These are outcome taus for each subject
  tau.strataO <- df_super$tau.strataO
  
  N <- nrow(df_super)
  
  # "treat" is a placeholder  
  zmat <- as.matrix(df_super[,c("treat","z","z.treat","z.k","z.k.treat")])
  x <- df_super[,c(xname)]
  
  # Initiate as population with covariates 
  # Outcomes and treatment assignment appended below
  dfsim <- df_super 
  
  # Randomization stratification
  strataR <- df_super[,c(strata_rand)]
  
  # Outcome stratification
  if(!is.null(strata_tte)) strataO <- df_super[,c(strata_tte)]
  
  if(is.null(strata_tte))  strataO <- "All"  

  # Set treatment to 1
  zmat.1 <- zmat
  zmat.1[,"treat"] <- 1.0
  zmat.1[,"z.treat"] <- zmat.1[,"z"]
  zmat.1[,"z.k"] <- zmat.1[,"z.k"]
  zmat.1[,"z.k.treat"] <- zmat.1[,"z.k"]
  
  # Set treatment to 0
  zmat.0 <- zmat
  zmat.0[,"treat"] <- 0.0
  zmat.0[,"z.treat"] <- 0.0
  zmat.0[,"z.k"] <- zmat.0[,"z.k"]
  zmat.0[,"z.k.treat"] <- 0.0
  
  # error term (log-scale [linear predictor scale])
  epsilon <- log(rexp(N))
  
  # linear predictor (setting treat=1)
  z.1 <- zmat.1[,"z"]
  # log(Y) AFT parameterization
  # psi.true are AFT parameters 
  eta1 <- mu + c(zmat.1%*%gamma.true)+x*bx
  # Weibull log(hazard ratio) setting treat=1
  # and excluding mu and x*bx (since taking difference below)
  phi1 <- -c(zmat.1%*%gamma.true)/tau.strataO
  log.Y1 <- eta1 + tau.strataO*epsilon
  
  # Setting treat=0
  z.0 <- zmat.0[,"z"]
  eta0 <- mu + c(zmat.0%*%gamma.true)+x*bx
  log.Y0 <- eta0 + tau.strataO*epsilon
  phi0 <- -c(zmat.0%*%gamma.true)/tau.strataO
  
  # Potential outcome log(hr) difference
  loghr.po <- phi1-phi0
  
  # Randomize per strataR
    blocks <- strataR
  
  # Randomize treatment
  Zr <- block_ra(blocks=blocks)
  
  log.Yr <- Zr*log.Y1 + (1-Zr)*log.Y0
  
  Yr <- exp(log.Yr)
  
  # At this point censoring is completely independent
  epsilonC <- log(rexp(N))
  log.YC <- muC + tauC*epsilonC
  
  log.YrC <- Zr*log.YC+(1-Zr)*log.YC 
  
  YrC <- exp(log.YrC)
  
  # Observed censored outcomes
  
  dfsim$event.sim <- ifelse(Yr<=YrC,1,0)
  dfsim$y.sim <- pmin(Yr,YrC)
  dfsim$treat.sim <- Zr
  dfsim$strata.simR <- strataR
  dfsim$strata.simO <- strataO
  
  dfsim$x <- x
  
  dfsim$loghr.po <- loghr.po
  dfsim$log.Y1 <- log.Y1
  dfsim$log.Y0 <- log.Y0
  
  if(details & ss <= 10) cat("% censored =",mean(1-dfsim$event.sim),"\n")
  
  # If checking compare simulated estimates with super-population
  
  if(checking){
    cat("Stratification parm (taus) df_super",c(tau),"\n")
    # strata 
    aa <- paste("strata(",eval("strata.simO"),")")
    bb <- c("Surv(y.sim,event.sim) ~ treat.sim + z + z.treat + z.k + z.k.treat +")
    weib.formula <- as.formula(paste(bb,aa))
    fitit <- survreg(weib.formula, dist='weibull', data=dfsim)
    fittau <- c(fitit$scale)
    cat("Stratification parm (taus) simulated=",c(fittau),"\n")
    # Check loghr.po = (log.Y1-log.Y0)/tau.strata0
    dcheck <- loghr.po - (log.Y0-log.Y1)/tau.strataO
    cat("Max |loghr.po - (log.Y0-log.Y1)/tau| = ",c(max(abs(round(dcheck,12)))),"\n")
    }
  
# Return sorted by biomarker z
dfs <- data.table::setorder(dfsim,z)
  
if(!return_df){
res <- list()
ahr_empirical <- with(dfs,exp(mean(loghr.po)))
res$AHR <- ahr_empirical
# AHR by X (X=0, and X=1)  
dfs_x1 <- subset(dfs,x==1)
res$AHR_X1 <- with(dfs_x1,exp(mean(loghr.po)))
dfs_x0 <- subset(dfs,x==0)
res$AHR_X0 <- with(dfs_x0,exp(mean(loghr.po)))

if(details) cat("Overall empirical AHR=",c(ahr_empirical),"\n")
# Three analyses: Standard ITT un-adjusted, Strata by R, Include xname
aa <- paste("strata(",eval("strata.simR"),")")
bb <- c("Surv(y.sim,event.sim) ~ treat.sim")
coxmod1 <- as.formula(bb)
bb <- c("Surv(y.sim,event.sim) ~ treat.sim +")
coxmod2 <- as.formula(paste(bb,aa))
bb <- c("Surv(y.sim,event.sim) ~ treat.sim + x")
coxmod3 <- as.formula(bb)
fit1 <- summary(coxph(coxmod1, data=dfsim))$conf.int
fit2 <- summary(coxph(coxmod2, data=dfsim))$conf.int
fit3 <- summary(coxph(coxmod3, data=dfsim))$conf.int
if(details) cat("ITT: Un-adjusted, sR, sR+x",c(fit1[1],fit2[1],fit3[1]),"\n")
res$ITT_unadj <- fit1[1]
res$ITT_sR <- fit2[1]
res$ITT_sRx <- fit3[1]
df_x1 <- subset(dfsim,x==1)
fit <- summary(coxph(as.formula(coxmod1), data=df_x1))$conf.int
if(details) cat("X=1 Sub-population",c(fit[1]),"\n")
res$X_1 <- fit[1]
df_x0 <- subset(dfsim,x==0)
fit <- summary(coxph(as.formula(coxmod1), data=df_x0))$conf.int
if(details) cat("X=0 Sub-population",c(fit[1]),"\n")
res$X_0 <- fit[1]

cut.zero <- with(dfs,min(z[which(loghr.po < hrz_crit)]))
# consider "optimal"
dfs_opt <- subset(dfs,z>=cut.zero)
res$AHR_opt <- with(dfs_opt,exp(mean(loghr.po)))

# For zpoints calculates AHR for z > zpoints
# Eg, for z=2 calculate AHR for population z>2
#zpoints <- seq(min(dfs$z),quantile(dfs$z,0.75),by=1)
zpoints <- seq(min(dfs$z),max(dfs$z),by=1)
HR_zpoints <- rep(NA,length(zpoints))
HRminus_zpoints <- rep(NA,length(zpoints))

for(zindex in 1:length(zpoints)){
  zz <- zpoints[zindex]
  dfz <- subset(dfs,z>=zz)
  HR_zpoints[zindex] <- with(dfz,exp(mean(loghr.po)))
  dfz_minus <- subset(dfs,z<=zz)
  HRminus_zpoints[zindex] <- with(dfz_minus,exp(mean(loghr.po)))
}
res$zpoints <- zpoints
res$HR.zpoints <- HR_zpoints
res$HRminus.zpoints <- HRminus_zpoints
if(details){
par(mfrow=c(1,2))
  with(dfs,plot(z,loghr.po,type="s",lty=1,xlab="z",ylab="psi0(z)"))
  with(dfs,rug(z))
  abline(h=c(log.hrs))
  abline(h=c(log(ahr_empirical)),lwd=2,col="orange",lty=1)
  abline(v=c(cut.zero),lwd=1,col="blue",lty=2)
  plot(zpoints,HR_zpoints,type="s",xlab="z",ylab="AHR(z+)",lty=1,col="grey",lwd=3)
}
}
if(!return_df) return(res)
if(return_df)  return(dfs)
}

require("splines")

cox_cs_fit2 <- function(df,tte.name=c("os_time"),event.name=c("os_event"),treat.name=c("treat"),strata.name=c("stratum"),z.name=c("bm"),alpha=0.20,details=FALSE,boots=0,
                        xlab=c("z"),maxz = Inf, show_plot=TRUE,truebeta.name=NULL,ypadzero=0.01,ydel=0.25,cex_legend=1.0,zwindow=0.0,byz=1,
                        ylimit=NULL,main_title=NULL,cex_count=0.7){
  
  if(boots >0) stop("Bootstrap needs revision (Set boots=0)")  
  
  c_alpha <- qnorm(1-alpha/2)
  
  z <- df[,c(z.name)]
  
  Y <- df[,c(tte.name)]
  E <- df[,c(event.name)]
  Treat <- df[,c(treat.name)]
  
  if(!is.null(truebeta.name)){
    beta_true <- df[,c(truebeta.name)]  
  }
  
  # Set z values for predicted fit
  # Set stratar values as well (just needed for predictions)
  # However when summarizing the log(hazard ratio) these 
  # quantities do NOT depend on the stratam (or the referencing)
  
  # Here we use reference="zero", may add stratum later
  
  # larger stratum (placeholder)
  # Find largest stratum
  if(!is.null(strata.name)){
    stratar <- df[,c(strata.name)]
    tabs <- as.data.frame(with(df,table(df[,c(strata.name)])))
    idm <- which(tabs[,2]==max(tabs[,2]))
    stratar_profile <- as.character(c(tabs[idm,1]))
  }
  
  if(is.null(strata.name)){
    stratar <- rep("1",length(Y))
    stratar_profile <- c("1")
  }
  
  z_q <- c(quantile(z,c(0.75,0.80,0.90),na.rm=TRUE))
  
  if(details) cat("Z quantiles at 75%, 80%, 90%:",c(z_q),"\n")
  
  # log(HR) estimates as function of Z estimated at:
  
  z_profile <- seq(min(z),z_q[3],by=byz)
  # cutoff at maxz
  z_profile <- unique(pmin(z_profile,rep(maxz,length(z_profile))))
  
  # Standard Cox
  fit0 <- coxph(Surv(Y,E) ~ Treat + strata(stratar))
  bhat0 <- coef(fit0)[1]
  # Will add bhat0 to plots
  
  #fit <- coxph(Surv(Y,E) ~ Treat*ns(z, df=3) + strata(stratar))
  # Treatment profiles based on treatment and z-nonlin fit
  # Potential-outcome predictions (setting treatment) at z values for specific strata
  # potential-outcomes prediction: 
  # assign combo=x (x=0,1) across observed cps
  # Use predict()
  # Setting treat=1 (combo=0)
  #df1 <- data.frame(z=z_profile,Treat=1,stratar=stratar_profile)
  #ypred <- predict(fit, type="lp", newdata=df1, se=TRUE, reference="strata")
  #yhat1 <- ypred$fit
  #se1 <- ypred$se
  #yhat1.lower <- yhat1-c_alpha*se1
  #yhat1.upper <- yhat1+c_alpha*se1
  # control 
  #df0 <- data.frame(z=z_profile,Treat=0,stratar=stratar_profile)
  #ypred <- predict(fit, type="lp", newdata=df0, se=TRUE, reference="strata")
  #yhat0 <- ypred$fit
  #se0 <- ypred$se
  #yhat0.lower <- yhat0-c_alpha*se0
  #yhat0.upper <- yhat0+c_alpha*se0
  # Look at 'beta':=log(ratio) ("potential outcome", po) of the above
  #beta_po <- yhat1-yhat0
  #sebeta_po <- sqrt(se1^2+se0^2)
  #beta_po_lower <- beta_po-c_alpha*sebeta_po
  #beta_po_upper <- beta_po+c_alpha*sebeta_po
  #ymax <- max(beta_po_upper)
  #ymin <- min(beta_po_lower)
  
  # Create ns() basis manually
  z_basis <- ns(z,df=3)
  z_mat <- as.matrix(z_basis)
  z_inter_mat <- Treat*z_mat
  cs_covs <- cbind(Treat,z_mat,z_inter_mat)
  fit <- coxph(Surv(Y,E) ~ cs_covs+strata(stratar))
  bhat <- fit$coefficients
  cov_bhat <- vcov(fit)
  
  # For predictions, need to construct ns() at desired cps values based on data cps basis
  # prediction cps
  z_ns <- predict(z_basis, z_profile)
  z_mat_profile <- as.matrix(z_ns) # Fixed (same for both arms)
  # Setting treatment to "treated"=1
  Trt1 <- rep(1,length(z_profile))
  z_inter_mat1 <- Trt1*z_mat_profile
  cs_covs_1 <- cbind(Trt1,z_mat_profile,z_inter_mat1)
  yhat1 <- c(cs_covs_1%*%bhat)
  
  # Setting treatment to control
  Trt0 <- rep(0,length(z_profile))
  z_inter_mat0 <- Trt0*z_mat_profile
  cs_covs_0 <- cbind(Trt0,z_mat_profile,z_inter_mat0)
  yhat0 <- c(cs_covs_0%*%bhat)
  
  
  # SEs for difference: 
  # calculate difference est. directly
  # Differences do not depend on "cps_mat"
  
  dhat <- bhat[1] + z_inter_mat1%*%c(bhat[-c(1,2,3,4)])
  
  # Calculate SEs
  # Extract cov(bhat) components involved in above difference
  # NOTE that SEs are lower for direct calculation of the difference
  # compared to above use of predict function since the parameters
  # for "z" are not needed in difference estimates (they cancel)
  # We will verify out asymptotic SEs with bootstrap below
  
  cov1_d <- cov_bhat[1,c(1,5,6,7)]
  cov5_d <- cov_bhat[5,c(1,5,6,7)]
  cov6_d <- cov_bhat[6,c(1,5,6,7)]
  cov7_d <- cov_bhat[7,c(1,5,6,7)]
  cov_bhat_diff <- cbind(cov1_d,cov5_d,cov6_d,cov7_d)
  zd <- cbind(rep(1,length(z_profile)),z_inter_mat1)
  # Asymptotic variance 
  var_dhat <- diag(zd%*%cov_bhat_diff%*%t(zd))
  se_dhat <- sqrt(var_dhat)
  
  dhat_lower <- dhat-c_alpha*se_dhat
  dhat_upper <- dhat+c_alpha*se_dhat
  
  ymax <- max(dhat_upper)
  ymin <- min(dhat_lower)
  
  
  # Check with bootstrap
  if(boots > 0){
    df1 <- data.frame(z=z_profile,Treat=1,stratar=stratar_profile)
    df0 <- data.frame(z=z_profile,Treat=0,stratar=stratar_profile)
    # Note: Point estimates are same but SEs a bit different
    # Let's try bootstrap to confirm
    set.seed(8316951)
    dhats_boot <- matrix(NA,nrow=length(z_profile),ncol=boots)
    N <- length(Y)
    for(bb in 1:boots){
      id_boot <- sample(c(1:N),size=N,replace=TRUE)
      dfB <- df[id_boot,]
      
      zb <- dfB[,c(z.name)]
      Yb <- dfB[,c(tte.name)]
      Eb <- dfB[,c(event.name)]
      Treatb <- dfB[,c(treat.name)]
      
      if(!is.null(strata.name)){
        stratarb <- dfB[,c(strata.name)]
      }
      if(is.null(strata.name)){
        stratarb <- rep("1",length(Yb))
      }
      fit_boot <- coxph(Surv(Yb,Eb) ~ Treatb*ns(zb, df=3) + strata(stratarb))
      ypreda <- predict(fit_boot, type="lp", newdata=df1, se=TRUE, reference="zero")
      ypredb <- predict(fit_boot, type="lp", newdata=df0, se=TRUE, reference="zero")
      dhat_boot <-ypreda$fit - ypredb$fit 
      dhats_boot[,bb] <- dhat_boot
    }
    var_dhats <- apply(dhats_boot,1,var)
    dhat_lower_boot <- dhat-c_alpha*sqrt(var_dhats)
    dhat_upper_boot <- dhat+c_alpha*sqrt(var_dhats)
  }
  
  if(show_plot){
    #plot(z_profile,dhat,ylim=c(ymin-0.1,ymax+0.1),type="l",lty=1,col="grey",lwd=2,xlab=xlab,ylab="log(hazard ratio) approximation")
    #abline(h=0.0,col="orange",lty=1,lwd=0.25)
    #abline(h=bhat0,col="blue",lty=1,lwd=0.5)
    #lines(z_profile,dhat_lower,lty=2,col="grey",lwd=1)
    #lines(z_profile,dhat_upper,lty=2,col="grey",lwd=1)
    #if(boots>0){
    #lines(z_profile,dhat_lower_boot,lty=2,col="black",lwd=1)
    #lines(z_profile,dhat_upper_boot,lty=2,col="black",lwd=1)
    #}
    #rug(z_profile)
    #legend("topright", c("log(HR) est.","95% lower/upper bound","Cox primary"),lty=c(1,2,1), lwd=c(2,1,0.5), col=c("grey","grey","blue"), bty='n')
    
    
    # Include counts of cps along x-axis
    # counts at z_profile
    # Around zwindow +/-
    counts_profile <- c(unlist(lapply(z_profile,function(x){sum(z>=(x-zwindow) & z<=(x+zwindow),na.rm=TRUE)})))
    # Add enough room at bottom for counts by ydel
    
    if(!is.null(truebeta.name)){
      ymin <- min(ymin,beta_true)
      ymax <- max(ymax,beta_true)
    }
    
    yymin <- ymin-ydel
    yymax <- ymax+ydel
    yystart <- round(log(0.5),1)
    yyend <- round(yymax,1)
    ypoints <- c(seq(yystart,yyend,length=10))
    yminzero <- round(ymin,1)-ypadzero
    ypoints <- sort(c(round(ypoints,1)))
    
    if(is.null(ylimit)) ylimit <- c(yymin,yymax)
    
    plot(z_profile,dhat,ylim=ylimit,type="l",lty=2,col="black",
         lwd=3,xlab=xlab,ylab="log(hazard ratio)",axes=FALSE)
    
    if(!is.null(truebeta.name)) lines(z,beta_true,type="l",lty=1,col="orange",lwd=2) 
    
    if(boots>0){
      lines(z_profile,dhat_lower_boot,lty=1,col="grey",lwd=2)
      lines(z_profile,dhat_upper_boot,lty=1,col="gry",lwd=2)
    }
    abline(h=0.0,col="red",lty=2,lwd=0.25)
    abline(h=bhat0,col="blue",lty=2,lwd=2)
    abline(h=log(0.80),col="red",lty=1,lwd=1)
    lines(z_profile,dhat_lower,lty=2,col="black",lwd=1)
    lines(z_profile,dhat_upper,lty=2,col="black",lwd=1)
    legend("top", c("log(HR)","80% CI"),lty=c(1,2), lwd=c(2,1), col=c("black","black"), bty='n',cex=cex_legend)
    legend("topright",c("log(0.80)","Cox primary"),lty=c(1,2),lwd=c(1,2),col=c("red","blue"),bty="n",cex=cex_legend)
    axis(2,at=c(ypoints),las=1)
    axis(1,at=c(z_profile),labels=c(z_profile))
    # horizontal line at ymin
    abline(h=yminzero,lty=1,col=1)
    box()
    # between yymin and yminzero
    yycounts <- round(yymin+(yminzero-yymin)/2,2)
    text(c(z_profile),c(yycounts),c(counts_profile),col="black",cex=cex_count)
  }
  if(!is.null(main_title)){
    title(main=c(main_title))
  } 
  
  
  return(list(z_profile=z_profile,dhat=dhat,dhat_lower=dhat_lower,dhat_upper=dhat_upper,counts_profile=counts_profile))
}