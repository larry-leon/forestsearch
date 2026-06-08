# Standalone cross-check: base-R Breslow Cox only (no forestsearch needed).
# gbsg ships in the survival package.
data(cancer, package = "survival")   # loads gbsg, rotterdam, etc.
stopifnot(exists("gbsg"))
gbsg$Y <- gbsg$rfstime / 30.4375          # months
gbsg$Event <- gbsg$status
gbsg$Treat <- gbsg$hormon                 # hormonal therapy indicator

# ---- base-R Breslow Cox (single binary covariate) ----
cox_bin <- function(y, d, z) {
  et <- sort(unique(y[d == 1])); if (!length(et)) return(NULL)
  Y1 <- sapply(et, function(t) sum(z[y >= t] == 1)); Y0 <- sapply(et, function(t) sum(z[y >= t] == 0))
  s1 <- sapply(et, function(t) sum(z[y == t & d == 1])); dk <- sapply(et, function(t) sum(y == t & d == 1))
  if (sum(dk) < 2 || all(Y1 == 0) || all(Y0 == 0)) return(NULL)
  U <- function(b){eb<-exp(b); zb<-(Y1*eb)/(Y1*eb+Y0); sum(s1-dk*zb)}
  rt <- tryCatch(uniroot(U, c(-15,15), extendInt="yes", tol=1e-9), error=function(e) NULL); if (is.null(rt)) return(NULL)
  b <- rt$root; eb<-exp(b); zb<-(Y1*eb)/(Y1*eb+Y0); info<-sum(dk*zb*(1-zb))
  if (!is.finite(b)||info<=0) return(NULL)
  dLam0<-dk/(Y1*eb+Y0); zbar_at<-approx(et,zb,xout=y,rule=2)$y
  L <- d*(z-zbar_at) - sapply(seq_along(y), function(i){idx<-which(et<=y[i]); exp(b*z[i])*sum((z[i]-zb[idx])*dLam0[idx])})
  list(beta=b, info=info, dfbeta=L/info, n=length(y), d=sum(d))
}
true_rate <- function(y,d,z,hr.c=1.0,R=2000,seed=1){
  set.seed(seed); thr<-log(hr.c); ns<-0L; nv<-0L; n<-length(y)
  for(r in 1:R){A<-rbinom(n,1,.5)==1
    if(sum(A)<5||sum(!A)<5||sum(d[A])<2||sum(d[!A])<2) next
    fa<-cox_bin(y[A],d[A],z[A]); fb<-cox_bin(y[!A],d[!A],z[!A]); if(is.null(fa)||is.null(fb)) next
    nv<-nv+1L; if(fa$beta>thr && fb$beta>thr) ns<-ns+1L}
  if(nv>0) ns/nv else NA_real_
}
approx_rate <- function(fit,hr.c=1.0,draws=4000,seed=1){
  thr<-log(hr.c); dfb<-fit$dfbeta; sD<-sqrt(sum(dfb^2)); delta<-fit$beta-thr
  closed<-max(0,2*pnorm(delta/sD)-1)
  set.seed(seed); G<-matrix(2*rbinom(draws*fit$n,1,.5)-1,nrow=draws); D<-as.numeric(G%*%dfb)
  c(HR=exp(fit$beta), sigma_D=sD, closed=closed, mc=mean((fit$beta-abs(D))>thr))
}
# the "flag" pattern: one entry point, choose direct splitting or approximation
consistency_rate <- function(y,d,z,hr.c=1.0,method=c("split","resample"),R=2000,draws=4000,seed=1){
  method<-match.arg(method)
  if(method=="split") return(true_rate(y,d,z,hr.c,R,seed))
  fit<-cox_bin(y,d,z); if(is.null(fit)) return(NA_real_); unname(approx_rate(fit,hr.c,draws,seed)["closed"])
}

cat("=== GBSG: paper's H-hat subgroup  Estrogen (er) <= 0 ===\n")
sub <- subset(gbsg, er <= 0)
fit <- cox_bin(sub$Y, sub$Event, sub$Treat)
tr  <- true_rate(sub$Y, sub$Event, sub$Treat, hr.c=1.0, R=4000, seed=7)
ar  <- approx_rate(fit, hr.c=1.0, draws=8000, seed=7)
cat(sprintf("  n=%d  events=%d  Cox HR (Breslow)=%.3f  robust sigma_D=%.3f\n",
            fit$n, fit$d, exp(fit$beta), ar["sigma_D"]))
cat(sprintf("  consistency rate:  direct-split=%.3f   approx-closed=%.3f   approx-mc=%.3f\n",
            tr, ar["closed"], ar["mc"]))
cat(sprintf("  flag check: consistency_rate(method='split')=%.3f  method='resample'=%.3f\n",
            consistency_rate(sub$Y,sub$Event,sub$Treat,method="split",R=4000,seed=7),
            consistency_rate(sub$Y,sub$Event,sub$Treat,method="resample",seed=7)))

cat("\n=== Candidate sweep: single-factor GBSG subgroups (screening HR>=1.25 context) ===\n")
gbsg$age_hi  <- as.integer(gbsg$age  > median(gbsg$age))
gbsg$nodes_hi<- as.integer(gbsg$nodes> median(gbsg$nodes))
gbsg$pgr_lo  <- as.integer(gbsg$pgr <= median(gbsg$pgr))
gbsg$size_hi <- as.integer(gbsg$size > median(gbsg$size))
cands <- list("er<=0"=gbsg$er<=0, "grade==3"=gbsg$grade==3, "meno==1"=gbsg$meno==1,
              "nodes>med"=gbsg$nodes_hi==1, "pgr<=med"=gbsg$pgr_lo==1,
              "size>med"=gbsg$size_hi==1, "age>med"=gbsg$age_hi==1)
cat(sprintf("%-12s %5s %6s %7s | %8s %8s %8s\n","subgroup","n","HR","sigmaD","split","closed","mc"))
cat(strrep("-",62),"\n")
for(nm in names(cands)){
  s<-gbsg[cands[[nm]],]; f<-cox_bin(s$Y,s$Event,s$Treat); if(is.null(f)) next
  t<-true_rate(s$Y,s$Event,s$Treat,1.0,R=3000,seed=7); a<-approx_rate(f,1.0,6000,seed=7)
  cat(sprintf("%-12s %5d %6.2f %7.3f | %8.3f %8.3f %8.3f\n",
              nm, f$n, exp(f$beta), a["sigma_D"], t, a["closed"], a["mc"]))
}
