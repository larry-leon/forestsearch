

get_dgm_stratified <- function(df,knot=5,kappa=10,hrs=c(0.75,0.75,0.75),details=TRUE,strata_tte=NULL,confounders=NULL){
  
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
sigC <- c(fitC.weib$scale)
muC <- c(coef(fitC.weib)[1])
  
  # Increase strength of treatment terms
  fit.sims <- fit.weibk

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
  
# At the change-point
a0 <- c(-gamma)
a0[1] <- loghr.0
a0[3] <- (loghr.knot-loghr.0)/knot
a0[5] <- (loghr.kappa-a0[1]-kappa*a0[3])/knot
# back-solve gamma
gamma <- -a0

    
  # Weibull log-scale
  b.true <- gamma

  # predictions log(hr)
  dfp <- dfa2[order(dfa2$z,decreasing=FALSE),]

return(list(df_super=dfp,b.true=b.true,mu=mu,sig=sig,muC=muC,sigC=sigC,strata_tte=strata_tte))
}

# Include single covariate X (Eg; "ecogbl", "prior_line12")
# strata_tte
# strata_rand

draw_sim_stratified <- function(dgm,ss=1,details=FALSE,Ndraw=nrow(dgm$df_super),strata_rand=c("strataNew"),xname=c("ecogbl"),bx=0,checking=FALSE){

df_super <- dgm$df_super
# Outcomes
b.true <- dgm$b.true
mu <- dgm$mu
sig <- dgm$sig
# Censoring
muC <- dgm$muC 
sigC <- dgm$sigC

strata_tte <- dgm$strata_tte

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

# Randomization stratification
strataR <- df_super[,c(strata_rand)]
# Outcome stratificaton
strataO <- df_super[,c(strata_tte)]

# Extract scales corresponding to strataO

aa <- names(sig)
sig_ids <- unlist(lapply(strataO,function(x){grep(x,aa)}))
sig_strataO <- sig[sig_ids]

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
# excluding baseline factor 
phi1 <- -c(zmat.1%*%b.true+x*bx)/sig_strataO
log.Y1 <- eta1 + sig_strataO*epsilon
  

z.0 <- zmat.0[,"z"]
eta0 <- mu + c(zmat.0%*%b.true)+x*bx
log.Y0 <- eta0 + sig_strataO*epsilon
phi0 <- -c(zmat.0%*%b.true+x*bx)/sig_strataO
# Potential outcome difference
loghr.po <- phi1-phi0
  

  # Randomize per strataR

  blocks <- strataR
  
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
  dfsim$strata.simR <- strataR
  dfsim$strata.simO <- strataO
  
  dfsim$x <- x
    
  dfsim$loghr.po <- loghr.po
  dfsim$log.Y1 <- log.Y1
  dfsim$log.Y0 <- log.Y0
  
  if(details & ss <= 10) cat("% censored (first 10)=",mean(1-dfsim$event.sim),"\n")
 
  # If checking compare simulated estimates with super-population
  
  
  if(checking){
  cat("Stratification parm (sigs) df_super",c(sig),"\n")
    # strata 
      aa <- paste("strata(",eval("strata.simO"),")")
      bb <- c("Surv(y.sim,event.sim) ~ treat.sim + z + z.treat + z.k + z.k.treat +")
      weib.formula <- as.formula(paste(bb,aa))
     fitit <- survreg(weib.formula, dist='weibull', data=dfsim)
    fitsig <- c(fitit$scale)
    cat("Stratification parm (sigs) simulated=",c(fitsig),"\n")
  }
  
  
  return(dfsim)
}





require("splines")

cox_cs_fit2 <- function(df,tte.name=c("os_time"),event.name=c("os_event"),treat.name=c("treat"),strata.name=c("stratum"),z.name=c("bm"),alpha=0.20,details=FALSE,boots=0,
                       xlab=c("z"),maxz = Inf, show_plot=TRUE,truebeta.name=NULL,ypadzero=0.01,ydel=0.175,cex_legend=1.0,zwindow=0.0,byz=1,ylimit=NULL,main_title=NULL,cex_count=0.7){
  
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
