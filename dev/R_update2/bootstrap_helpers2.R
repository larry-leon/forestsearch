#' Count ID Occurrences in Bootstrap Sample
#'
#' Counts the number of times an ID appears in a bootstrap sample.
#'
#' @param x ID value.
#' @param dfb Data frame of bootstrap sample.
#' @return Integer count of occurrences.
#' @export

count.id <- function(x,dfb){
  sum(dfb$id==x)
  }


#' Confidence Interval for Estimate
#'
#' Calculates confidence interval for an estimate, optionally on log(HR) scale.
#'
#' @param x Numeric estimate.
#' @param sd Numeric standard deviation.
#' @param alpha Numeric significance level (default: 0.025).
#' @param scale Character. "hr" or "1/hr".
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return List with length, lower, upper, sd, and estimate.
#' @importFrom stats qnorm
#' @export

ci_est <- function(x, sd, alpha = 0.025, scale = "hr", est.loghr = TRUE) {
  # Input validation
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || is.infinite(x)) {
    stop("'x' must be a single, finite numeric value.")
  }
  if (!is.numeric(sd) || length(sd) != 1 || is.na(sd) || is.infinite(sd) || sd < 0) {
    stop("'sd' must be a single, non-negative, finite numeric value.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 0.5) {
    stop("'alpha' must be a single numeric value between 0 and 0.5.")
  }
  if (!scale %in% c("hr", "1/hr")) {
    stop("'scale' must be either 'hr' or '1/hr'.")
  }
  if (!is.logical(est.loghr) || length(est.loghr) != 1 || is.na(est.loghr)) {
    stop("'est.loghr' must be a single logical value.")
  }

  c_alpha <- qnorm(1 - alpha)
  c_low <- x - c_alpha * sd
  c_up <- x + c_alpha * sd
  est <- x
  new_low <- c_low
  new_up <- c_up
  out_sd <- sd

  if (scale == "hr" && est.loghr) {
    est <- exp(x)
    out_sd <- exp(x) * sd
    new_low <- exp(c_low)
    new_up <- exp(c_up)
  }
  if (scale == "1/hr" && est.loghr) {
    est <- exp(-x)
    out_sd <- exp(-x) * sd
    new_low <- exp(-c_up)
    new_up <- exp(-c_low)
  }
  length <- new_up - new_low

  # Return as named list
  return(list(
    length = length,
    lower = new_low,
    upper = new_up,
    sd = out_sd,
    est = est
  ))
}


#' Fit Cox Model for Subgroup
#'
#' Fits a Cox model for a subgroup and returns estimate and standard error.
#'
#' @param df_sg Data frame for subgroup.
#' @param cox.formula Cox model formula.
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return List with estimate and standard error.
#' @importFrom survival coxph
#' @export

get_Cox_sg <- function(df_sg, cox.formula, est.loghr = TRUE) {
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)
  if (sum(!is.na(check2)) != length(names_tocheck)) stop("df_sg dataset NOT contain cox.formula variables")
  # Fit Cox model with robust standard errors
  fit <- summary(coxph(cox.formula, data = df_sg, robust = TRUE))$coefficients
  # log(hr) parameters
  if (est.loghr) {
    bhat <- c(fit[, "coef"])
    est_obs <- bhat
    se_obs <- c(fit[, "robust se"])
  }
  # Otherwise, hr
  if (!est.loghr) {
    bhat <- c(fit[, "coef"])
    est_obs <- exp(bhat)
    sebhat <- c(fit[, "robust se"])
    se_obs <- est_obs * sebhat
  }
  return(list(est_obs = est_obs, se_obs = se_obs))
}


#' Coverage Indicator for Confidence Interval
#'
#' Checks if a target value is covered by a confidence interval.
#'
#' @param lower Numeric lower bound.
#' @param upper Numeric upper bound.
#' @param target Numeric target value (default: 0).
#' @return List with coverage indicator and interval length.
#' @export

ci_cover<-function(lower,upper,target=0){
  if(length(target)==1){
    cover<-ifelse(lower<=target & upper>=target,1,0)
    LC<-upper-lower
    return(list(cover=cover,LC=LC))
  }
if(length(target)>1){
LC<-upper-lower
covers<-NULL
for(tt in 1:length(target)){
cover<-ifelse(lower<=target[tt] & upper>=target[tt],1,0)
covers<-c(covers,cover)
}
return(list(cover=covers,LC=LC))
  }
  }

#' Bootstrap Confidence Intervals for Two Targets
#'
#' Calculates confidence intervals for two bootstrap targets.
#'
#' @param Q1 Numeric vector of first target estimates.
#' @param Q2 Numeric vector of second target estimates.
#' @param ystar Matrix of bootstrap samples.
#' @return List with estimates, standard errors, and confidence intervals.
#' @export

getCIs<-function(Q1,Q2,ystar){
# Target 1
est<-get_targetEst(x=Q1,ystar=ystar)
# If est.loghr=TRUE then on log(HR) scale
q1<-est$target_est
#se1<-est$sehat
se1_new<-est$sehat_new
se1<-est$sehat
rm("est")
# use SE new
cest<-ci_est(x=q1,sd=se1_new)
H1_lower<-cest$lower
H1_upper<-cest$upper
rm("cest")

# Target 2
est<-get_targetEst(x=Q2,ystar=ystar)
# Call this H.bc
q2<-est$target_est
se2<-est$sehat
se2_new<-est$sehat_new
rm("est")
cest<-ci_est(x=q2,sd=se2_new)
H2_lower<-cest$lower
H2_upper<-cest$upper
rm("cest")
return(list(q1=q1,se1=se1,se1_new=se1_new,H1_lower=H1_lower,H1_upper=H1_upper,
            q2=q2,se2=se2,se2_new=se2_new,H2_lower=H2_lower,H2_upper=H2_upper))
}

#' Prepare Data for Bias Plot
#'
#' Prepares a data frame for plotting bias in bootstrap estimates.
#'
#' @param res Data frame of results.
#' @param dgm Data-generating mechanism (truth) for simulation.
#' @return Data frame for bias plot.
#' @importFrom data.table data.table
#' @export

get_dfPlot<-function(res,dgm){
  hrH_true<-dgm$hr.H.true
  hrHc_true<-dgm$hr.Hc.true
  if(est.loghr & est.scale=="loghr"){
    hrH_true<-log(dgm$hr.H.true)
    hrHc_true<-log(dgm$hr.Hc.true)
  }
  res_new<-within(res,{
    b1H_1<-100*(H1.bc-H_true)/H_true
    b1H_2<-100*(H1.bc-hatH_causal)/hatH_causal
    b1H_3<-100*(H1.bc-hrH_true)/hrH_true

    b2H_1<-100*(H2.bc-H_true)/H_true
    b2H_2<-100*(H2.bc-hatH_causal)/hatH_causal
    b2H_3<-100*(H2.bc-hrH_true)/hrH_true

    b0H_1<-100*(H_obs-H_true)/H_true
    b0H_2<-100*(H_obs-hatH_causal)/hatH_causal
    b0H_3<-100*(H_obs-hrH_true)/hrH_true
  }
  )

  df_bc<-NULL
  # hr(H_true) target
  res_new<-data.table(res_new)
  df_res<-res_new[,c("b1H_1")]
  df_res$est<-"BC_1"
  df_res$target<-"H_true"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b2H_1")]
  df_res$est<-"BC_2"
  df_res$target<-"H_true"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b0H_1")]
  df_res$est<-"Obs"
  df_res$target<-"H_true"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  # hatH_causal target
  df_res<-res_new[,c("b1H_2")]
  df_res$est<-"BC_1"
  df_res$target<-"H_po"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b2H_2")]
  df_res$est<-"BC_2"
  df_res$target<-"H_po"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b0H_2")]
  df_res$est<-"Obs"
  df_res$target<-"H_po"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  # H_causal fixed
  df_res<-res_new[,c("b1H_3")]
  df_res$est<-"BC_1"
  df_res$target<-"H_fixed"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b2H_3")]
  df_res$est<-"BC_2"
  df_res$target<-"H_fixed"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)

  df_res<-res_new[,c("b0H_3")]
  df_res$est<-"Obs"
  df_res$target<-"H_fixed"
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  return(df_bc)
}

#' Variance Summary for Bootstrap Results
#'
#' Summarizes variance and coverage for bootstrap results.
#'
#' @param res Data frame of results.
#' @return Data frame with summary statistics.
#' @export

var_summary<-function(res){
  df_var<-NULL
  # BC1 SD and Avg(est(sd))
  aa<-with(res,sqrt(var(H1.bc,na.rm=TRUE)))
  bb<-with(res,mean(seH1.bc,na.rm=TRUE))
  cc<-with(res,mean(H1.bc,na.rm=TRUE))
  dd<-with(res,mean(H1_cover1,na.rm=TRUE))
  ee<-with(res,mean(H1_cover2,na.rm=TRUE))
  ff<-with(res,mean(H1_cover3,na.rm=TRUE))
  gg<-with(res,mean(H1_length,na.rm=TRUE))
  df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

  aa<-with(res,sqrt(var(H2.bc,na.rm=TRUE)))
  bb<-with(res,mean(seH2.bc,na.rm=TRUE))
  cc<-with(res,mean(H2.bc,na.rm=TRUE))
  dd<-with(res,mean(H2_cover1,na.rm=TRUE))
  ee<-with(res,mean(H2_cover2,na.rm=TRUE))
  ff<-with(res,mean(H2_cover3,na.rm=TRUE))
  gg<-with(res,mean(H2_length,na.rm=TRUE))
  df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

  aa<-with(res,sqrt(var(Hc1.bc,na.rm=TRUE)))
  bb<-with(res,mean(seHc1.bc,na.rm=TRUE))
  cc<-with(res,mean(Hc1.bc,na.rm=TRUE))
  dd<-with(res,mean(Hc1_cover1,na.rm=TRUE))
  ee<-with(res,mean(Hc1_cover2,na.rm=TRUE))
  ff<-with(res,mean(Hc1_cover3,na.rm=TRUE))
  gg<-with(res,mean(Hc1_length,na.rm=TRUE))
  df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

  aa<-with(res,sqrt(var(Hc2.bc,na.rm=TRUE)))
  bb<-with(res,mean(seHc2.bc,na.rm=TRUE))
  cc<-with(res,mean(Hc2.bc,na.rm=TRUE))
  dd<-with(res,mean(Hc2_cover1,na.rm=TRUE))
  ee<-with(res,mean(Hc2_cover2,na.rm=TRUE))
  ff<-with(res,mean(Hc2_cover3,na.rm=TRUE))
  gg<-with(res,mean(Hc2_length,na.rm=TRUE))
  df_var<-rbind(df_var,c(cc,aa,bb,dd,ee,ff,gg))

  colnames(df_var)<-c("Est","SD","Est(SD)","C1","C2","C3","L")
  rownames(df_var)<-c("H1_bc","H2_bc","Hc1_bc","Hc2_bc")
  return(round(df_var,digits=3))
}


#' Confidence Interval for Cox Model Estimate
#'
#' Calculates confidence interval and coverage for Cox model estimate.
#'
#' @param df Data frame.
#' @param est Character. Name of estimate column.
#' @param se Character. Name of standard error column.
#' @param target Numeric or character. Target value or column name.
#' @param alpha Numeric significance level (default: 0.025).
#' @param digits Integer. Number of digits for rounding (default: 3).
#' @return List with lower bound, upper bound, and coverage indicator.
#' @export

getci_Cox<-function(df,est,se,target,alpha=0.025,digits=3){
  a<-df[est]
  b<-df[se]
  if(!is.numeric(target)){
    c<-df[target]
  }
  if(is.numeric(target)){
    c<-target
  }
  log_lb<-log(a)-qnorm(0.975)*(b/a)
  log_ub<-log(a)+qnorm(0.975)*(b/a)
  lb<-exp(log_lb)
  ub<-exp(log_ub)
  cov<-ifelse(c>=lb & c<=ub,1,0)
  return(list(lb=round(lb,3),ub=round(ub,3),cover=cov))
}

#' Summary Statistic for Data Frame Column
#'
#' Returns mean and standard deviation for a column, formatted as a string.
#'
#' @param df Data frame.
#' @param name Character. Column name.
#' @param sigdig Integer. Number of significant digits (default: 2).
#' @param includeSD Logical. Include "SD=" in output (default: FALSE).
#' @param showSD Logical. Show SD in output (default: TRUE).
#' @return Character string with mean and SD.
#' @export

SummaryStat<-function(df,name,sigdig=2,includeSD=FALSE,showSD=TRUE){
  if(is.data.table(df)){
    df<-as.data.frame(df)
  }
  # include "SD=" in output
  if(!includeSD){
    temp<-na.omit(df[,c(name)])
    m<-round(mean(temp),sigdig)
    sig<-round(sqrt(var(temp)),sigdig)
    # Check if binary
    if(all(as.numeric(temp) %in% c(0,1))){
      sig<-round(sqrt(m*(1-m)/length(temp)),sigdig)
    }
    if(sigdig<=4 & m<=0.001) m<-"<0.001"
    if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
    out<-paste0(m," (")
    out<-paste0(out,sig)
    out<-paste0(out,")")
  }
  if(includeSD){
    temp<-na.omit(df[,c(name)])
    m<-round(mean(temp),sigdig)
    sig<-round(sqrt(var(temp)),sigdig)
    # Check if binary
    if(all(as.numeric(temp) %in% c(0,1))){
      sig<-round(sqrt(m*(1-m)/length(m)),sigdig)
    }
    if(sigdig<=4 & m<=0.001) m<-"<0.001"
    if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
    out<-paste0(m," (SD=")
    out<-paste0(out,sig)
    out<-paste0(out,")")
  }
  if(!showSD) out<-c(m)
  return(out)
}

#' Difference in Rate Between Two Data Frames
#'
#' Calculates the difference in mean rate for a column between two data frames.
#'
#' @param df1 First data frame.
#' @param df2 Second data frame.
#' @param name Character. Column name.
#' @param sigdig Integer. Number of significant digits (default: 2).
#' @return Numeric value of difference in rate (percent).
#' @export

DiffRate<-function(df1,df2,name,sigdig=2){
  if(is.data.table(df1)){
    df1<-as.data.frame(df1)
    df2<-as.data.frame(df2)
  }
  temp<-na.omit(df1[,c(name)])
  m1<-mean(temp)
  temp<-na.omit(df2[,c(name)])
  m2<-mean(temp)
  out<-round(100*(m2-m1),1)
  return(out)
}

#' Power Calculation for Subgroup Size
#'
#' Calculates power (rejection rate) for subgroups above a minimum size.
#'
#' @param df Data frame.
#' @param minsize Integer. Minimum subgroup size (default: 0).
#' @param sigdig Integer. Number of significant digits (default: 3).
#' @return Numeric value of power (rejection rate).
#' @export

pow_size<-function(df,minsize=0,sigdig=3){
  rej<-with(df,mean(rej12 & size.Hc>=minsize))
  return(round(rej,sigdig))
}

