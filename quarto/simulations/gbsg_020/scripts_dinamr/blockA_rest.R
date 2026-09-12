suppressMessages(library(forestsearch)); options(width=230)
R <- "results/"
wil <- function(p,n,z=stats::qnorm(0.975)){c<-(p+z^2/(2*n))/(1+z^2/n);h<-z*sqrt(p*(1-p)/n+z^2/(4*n^2))/(1+z^2/n);c(c-h,c+h)}
cells <- list(c(1.50,500),c(1.50,1000),c(1.50,1500),c(1.75,500),c(1.75,1000),c(1.75,1500))
lab <- function(x) sprintf("HR %.2f n %-4d", x[1], as.integer(x[2]))
dfile <- function(x) sprintf("%sdina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_nb20_dinamr_combined_1_2000.rds", R, round(100*x[1]), as.integer(x[2]))
ffile <- function(x) sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*x[1]), as.integer(x[2]),
                             if (abs(x[1]-1.75)<1e-9) "tier2" else "p12ext")
ev <- function(r) r[r$detected %in% 1L & is.finite(r$betaHhat_H) & is.finite(r$betaHhat_Hc) &
                    is.finite(r$fld_H_lo1s) & is.finite(r$fld_Hc_up1s_s) & is.finite(r$mr_H_lo), , drop=FALSE]

cat("\n#### 2. DINA BESIDE THE FS COMPARATOR -- coverage, Wilson, IJ misses by side, detection\n")
cat("####    CONFOUND: different identifier AND different family construction; FS at 12.4% is\n")
cat("####    maxeffCons at eps 0.10 (p12ext / tier2), DINA is effMaxSG at eps 0.20.  Descriptive.\n\n")
out <- do.call(rbind, lapply(cells, function(x) do.call(rbind, lapply(c("DINA","FS"), function(who) {
  r <- readRDS(if (who=="DINA") dfile(x) else ffile(x))$results; d <- ev(r)
  cH <- mean(d$betaHhat_H >= d$fld_H_lo1s); wH <- wil(cH,nrow(d))
  cC <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s); wC <- wil(cC,nrow(d))
  ij <- mean(d$betaHhat_H>=d$mr_H_lo & d$betaHhat_H<=d$mr_H_hi); wI <- wil(ij,nrow(d))
  data.frame(cell=lab(x), engine=who, detection=mean(r$detected %in% 1L), n_eval=nrow(d),
             K_med=stats::median(r$n_family, na.rm=TRUE),
             fldH_cov=cH, fldH_lo=wH[1], fldH_hi=wH[2],
             fldHc_s_cov=cC, fldHc_s_lo=wC[1], fldHc_s_hi=wC[2],
             ij2_H=ij, ij2_lo=wI[1], ij2_hi=wI[2],
             miss_below=mean(d$betaHhat_H<d$mr_H_lo), miss_above=mean(d$betaHhat_H>d$mr_H_hi),
             ij2_Hc=mean(d$betaHhat_Hc>=d$mr_Hc_lo & d$betaHhat_Hc<=d$mr_Hc_hi)) }))))
print(out, row.names=FALSE, digits=4)

cat("\n\n#### 3. STRATIFIED: field LOWER-BOUND coverage and RETAINED BIAS on the harm block\n")
cat("####    by n_family stratum and by p-hat bin, with counts.\n")
cat("####    Bins = within-cell empirical TERTILES (the FS summaries' own tert(): quantile\n")
cat("####    c(0,1/3,2/3,1), include.lowest=TRUE).  K=1 and K<=5 OVERLAP the tertiles.\n")
tert_fs <- function(v){ if(!any(is.finite(v))) return(rep(NA_integer_,length(v)))
  br <- stats::quantile(v,c(0,1/3,2/3,1),na.rm=TRUE,names=FALSE)
  if(length(unique(br))<2L) return(rep(1L,length(v)))
  cut(v,breaks=unique(br),include.lowest=TRUE,labels=FALSE) }
row1 <- function(d,cell,st,ov){ if(!nrow(d)) return(NULL)
  cv <- mean(d$betaHhat_H>=d$fld_H_lo1s); w <- wil(cv,nrow(d))
  data.frame(cell=cell, stratum=st, overlapping=ov, count=nrow(d),
             K_med=stats::median(d$n_family), p_hat_med=stats::median(d$p_hat_H),
             fldH_cov=cv, wilson_lo=w[1], wilson_hi=w[2],
             retained_bias_log=mean(log(d$fld_H_est2)-log(d$betaHhat_H)),
             naive_bias_log=mean(log(d$nv_H_est)-log(d$betaHhat_H)),
             mean_fld_H_se=mean(d$fld_H_se)) }
for (by in c("K","p")) {
 cat(sprintf("\n--- by %s ---\n", if(by=="K") "n_family stratum" else "p-hat bin"))
 S <- do.call(rbind, lapply(cells, function(x){ d <- ev(readRDS(dfile(x))$results); cl <- lab(x)
   v <- if(by=="K") d$n_family else d$p_hat_H; g <- tert_fs(v)
   rr <- lapply(sort(unique(g[!is.na(g)])), function(t){ dk <- d[!is.na(g)&g==t,]
     row1(dk, cl, sprintf("%s T%d [%s, %s]", if(by=="K")"K" else "p-hat", t,
          format(min(v[!is.na(g)&g==t]),digits=4), format(max(v[!is.na(g)&g==t]),digits=4)), FALSE) })
   if (by=="K") rr <- c(rr, list(row1(d[d$n_family==1L,],cl,"K = 1  (overlapping)",TRUE),
                                 row1(d[d$n_family<=5L,],cl,"K <= 5 (overlapping)",TRUE)))
   do.call(rbind, c(rr, list(row1(d,cl,"all detected (overlapping)",TRUE)))) }))
 print(S, row.names=FALSE, digits=4) }

cat("\n\n#### 4. K-STRATUM x p-HAT-BIN COUNTS, per cell\n")
for (x in cells) { d <- ev(readRDS(dfile(x))$results)
  cat(sprintf("\n%s  (n_eval %d)\n", lab(x), nrow(d)))
  print(table(K=paste0("K T",tert_fs(d$n_family)), p_hat=paste0("p T",tert_fs(d$p_hat_H)))) }

cat("\n\n#### 5. MEDIAN FIELD LOWER BOUND vs MEDIAN TRUE theta(Hhat), HR SCALE\n")
cat("####    theta(Hhat) = betaHhat_H, the realized conditional target on the super-population.\n\n")
L <- do.call(rbind, lapply(cells, function(x){ d <- ev(readRDS(dfile(x))$results)
  data.frame(cell=lab(x), n_eval=nrow(d),
             med_fld_H_lo1s=stats::median(d$fld_H_lo1s),
             med_theta_Hhat=stats::median(d$betaHhat_H),
             med_fld_H_est2=stats::median(d$fld_H_est2),
             med_naive_H=stats::median(d$nv_H_est),
             bound_minus_theta=stats::median(d$fld_H_lo1s)-stats::median(d$betaHhat_H),
             bound_over_theta=stats::median(d$fld_H_lo1s)/stats::median(d$betaHhat_H),
             med_paired_ratio=stats::median(d$fld_H_lo1s/d$betaHhat_H),
             planted_marg_H=unname(readRDS(dfile(x))$truth$marg_H),
             share_bound_ge_1.00=mean(d$fld_H_lo1s>=1.00),
             share_bound_ge_1.25=mean(d$fld_H_lo1s>=1.25)) }))
print(L, row.names=FALSE, digits=4)
