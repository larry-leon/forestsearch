# Per-cell headline numbers for REPORT_grfmr_2026-09-11.md and, extended to the
# full 18-cell grid, REPORT_grfmr_completion_2026-09-12.md.
# Reads whatever grfmr combined bundles are on disk and emits the report's
# tables in markdown.  Re-runnable: cells absent from disk are skipped, never
# rendered empty.
#
# Coverage figures are computed over the detected replicates of the cell each
# row names.  Nothing here certifies GRF and no acceptance criterion is applied.
#
# HR 1.00 cells (grfmr completion): the rate is a SELECTION RATE and nothing
# more -- the planted region is differentially null against a benefiting
# complement and clears the sub-null log(0.90) floor, so returning it is an
# admissible selection.  Bound location carries the question, so the null
# cells get their own location table with FS beside GRF.
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
R <- file.path(QMD_DIR, "results/")
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||all(is.na(a))) b else a
wil <- function(p, n, conf = .95) { if (!is.finite(n)||n<=0) return(c(NA,NA))
  z <- stats::qnorm(1-(1-conf)/2)
  c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/(1+z^2/n),
    (p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/(1+z^2/n)) }
wf <- function(p,n) sprintf("%.4f [%.4f, %.4f]", p, wil(p,n)[1], wil(p,n)[2])

cells <- do.call(rbind, lapply(c("12.4%","31%"), function(pv)
  do.call(rbind, lapply(c(1.50,1.75,1.00), function(hr)
    data.frame(prev=pv, hr=hr, n=c(500L,1000L,1500L), stringsAsFactors=FALSE)))))
gf <- function(hr,n,pv) sprintf("%sgrf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfmr_combined_1_2000.rds",
                                R, round(100*hr), n, if (pv=="31%") "_z1q60" else "")
# The DESIGNATED comparator, as gate2G.R resolves it.  HR 1.00: tier2 at 12.4%
# n 500, p12ext at 12.4% n 1000/1500; cert20 at 31% for all n (no e1stud HR
# 1.00 bundle).  map1 / s7 bundles are never substituted.
fsf <- function(hr,n,pv) { if (pv=="12.4%") {
    camp <- if (abs(hr-1.75)<1e-9) "tier2" else if (abs(hr-1.00)<1e-9 && n==500L) "tier2" else "p12ext"
    sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
  } else { camp <- if (n==500L && abs(hr-1.00)>1e-9) "e1stud" else "cert20"
    sprintf("%sfs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_z1q60_nb20_%s_combined_1_2000.rds", R, round(100*hr), n, camp) } }
# Location, on the same row set fs_extraction.R uses: detected, bound and theta finite.
loc <- function(r) { d <- r[r$detected %in% 1L & is.finite(r$fld_H_lo1s) & is.finite(r$betaHhat_H), , drop=FALSE]
  list(n=nrow(d), med_bound=stats::median(d$fld_H_lo1s), med_theta=stats::median(d$betaHhat_H),
       paired=stats::median(d$fld_H_lo1s/d$betaHhat_H),
       s100=mean(d$fld_H_lo1s>=1.00), s125=mean(d$fld_H_lo1s>=1.25)) }

rows <- list()
for (i in seq_len(nrow(cells))) {
  f <- gf(cells$hr[i], cells$n[i], cells$prev[i]); if (!file.exists(f)) next
  b <- readRDS(f); r <- b$results; m <- b$meta
  d <- r[r$detected %in% 1L & is.finite(r$betaHhat_H) & is.finite(r$betaHhat_Hc) &
         is.finite(r$fld_H_lo1s) & is.finite(r$fld_Hc_up1s_s), , drop=FALSE]
  A <- r$admitted_n[is.finite(r$admitted_n)]; K <- r$n_family[is.finite(r$n_family)]
  q <- stats::quantile(A, c(.25,.50,.75,.90), names=FALSE)
  fb <- if (file.exists(fsf(cells$hr[i],cells$n[i],cells$prev[i])))
          readRDS(fsf(cells$hr[i],cells$n[i],cells$prev[i])) else NULL
  fd <- if (!is.null(fb)) { fr <- fb$results
    fr[fr$detected %in% 1L & is.finite(fr$betaHhat_H) & is.finite(fr$fld_H_lo1s) &
       is.finite(fr$fld_Hc_up1s_s), , drop=FALSE] } else NULL
  gl <- loc(r); fl <- if (!is.null(fb)) loc(fb$results) else NULL
  ND <- r[!(r$detected %in% 1L), , drop=FALSE]
  cov <- function(x) mean(x)
  rows[[length(rows)+1]] <- data.frame(
    prev=cells$prev[i], hr=cells$hr[i], n=cells$n[i],
    detection=mean(r$detected %in% 1L), n_eval=nrow(d),
    n_nondet=nrow(ND), nd_adm_na=sum(is.na(ND$admitted_n)), nd_adm_0=sum(ND$admitted_n %in% 0L),
    adm_min=min(A), adm_q25=q[1], adm_med=q[2], adm_q75=q[3], adm_p90=q[4], adm_max=max(A),
    adm_cv=stats::sd(A)/mean(A), adm_eq0=sum(A==0L), adm_na=sum(!is.finite(r$admitted_n)),
    K_min=min(K), K_med=stats::median(K), K_max=max(K), K_cv=stats::sd(K)/mean(K),
    share_adm=stats::median(r$admitted_n[is.finite(r$admitted_n)&is.finite(r$n_family)] /
                            r$n_family[is.finite(r$admitted_n)&is.finite(r$n_family)]),
    prev_trial=mean(r$n_true)/m$n_sample,
    nv_H=cov(d$betaHhat_H >= d$nv_H_lo & d$betaHhat_H <= d$nv_H_hi),
    fld_H=cov(d$betaHhat_H >= d$fld_H_lo1s),
    fld_Hc_s=cov(d$betaHhat_Hc <= d$fld_Hc_up1s_s),
    ij2_H=cov(d$betaHhat_H >= d$mr_H_lo & d$betaHhat_H <= d$mr_H_hi),
    ij2_Hc=cov(d$betaHhat_Hc >= d$mr_Hc_lo & d$betaHhat_Hc <= d$mr_Hc_hi),
    joint_s=cov(d$betaHhat_H >= d$fld_joint_s_bonf_loH & d$betaHhat_Hc <= d$fld_joint_s_bonf_upHc),
    sens=mean(d$sens), spec=mean(d$spec), ppv=mean(d$ppv), npv=mean(d$npv),
    mean_Hhat=mean(d$n_sel),
    med_bound=stats::median(d$fld_H_lo1s), med_theta=stats::median(d$betaHhat_H),
    share_ge100=mean(d$fld_H_lo1s>=1.00), share_ge125=mean(d$fld_H_lo1s>=1.25),
    L_n=gl$n, L_med_bound=gl$med_bound, L_med_theta=gl$med_theta, L_paired=gl$paired,
    L_s100=gl$s100, L_s125=gl$s125, planted=unname(b$truth$marg_H),
    fs_camp=if (!is.null(fb)) fb$meta$campaign_tag %||% NA else NA,
    fs_focus=if (!is.null(fb)) fb$meta$sg_focus %||% NA else NA,
    fs_eps=if (!is.null(fb)) fb$meta$effect_neighborhood %||% NA else NA,
    fs_det=if (!is.null(fd)) mean(fb$results$detected %in% 1L) else NA,
    fs_fld_H=if (!is.null(fd)) cov(fd$betaHhat_H >= fd$fld_H_lo1s) else NA,
    fs_ij2_H=if (!is.null(fd)) cov(fd$betaHhat_H >= fd$mr_H_lo & fd$betaHhat_H <= fd$mr_H_hi) else NA,
    fs_K_med=if (!is.null(fb)) stats::median(fb$results$n_family, na.rm=TRUE) else NA,
    fsL_n=if (!is.null(fl)) fl$n else NA, fsL_med_bound=if (!is.null(fl)) fl$med_bound else NA,
    fsL_med_theta=if (!is.null(fl)) fl$med_theta else NA, fsL_paired=if (!is.null(fl)) fl$paired else NA,
    fsL_s100=if (!is.null(fl)) fl$s100 else NA, fsL_s125=if (!is.null(fl)) fl$s125 else NA,
    stringsAsFactors=FALSE)
}
if (!length(rows)) { cat("No grfmr combined bundle on disk yet.\n"); quit(save="no") }
X <- do.call(rbind, rows)
cl <- function(i) sprintf("%s HR %.2f n %d", X$prev[i], X$hr[i], X$n[i])
isnull <- abs(X$hr - 1.00) < 1e-9

cat("\n### Detection (harm cells) / selection rate (HR 1.00 cells), admitted_n and the enumerated pool\n\n")
cat("| cell | detection or **selection rate** [Wilson] | n_eval | **admitted_n** min / q25 / **med** / q75 / p90 / max | adm CV | adm = 0 | adm NA | non-det (adm NA / adm 0) | n_family (enumerated pool) min / **med** / max | K CV | med share |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(X))) cat(sprintf(
 "| %s | %s | %d | %g / %g / **%g** / %g / %g / %g | %.3f | %d | %d | %d (%d / %d) | %g / **%g** / %g | %.4f | %.4f |\n",
 cl(i), wf(X$detection[i], 2000L), X$n_eval[i], X$adm_min[i], X$adm_q25[i], X$adm_med[i],
 X$adm_q75[i], X$adm_p90[i], X$adm_max[i], X$adm_cv[i], X$adm_eq0[i], X$adm_na[i],
 X$n_nondet[i], X$nd_adm_na[i], X$nd_adm_0[i],
 X$K_min[i], X$K_med[i], X$K_max[i], X$K_cv[i], X$share_adm[i]))

cat("\n### Coverage of every product, absolute levels [Wilson]\n\n")
cat("Over detected replicates throughout. At HR 1.00 the H block is the identified region, differentially null against a benefiting complement.\n\n")
cat("| cell | naive H | field H (1-sided) | field-s Hc (1-sided) | IJ 2-sided H | IJ 2-sided Hc | joint_s Bonferroni |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(X))) cat(sprintf("| %s | %s | %s | %s | %s | %s | %s |\n", cl(i),
 wf(X$nv_H[i],X$n_eval[i]), wf(X$fld_H[i],X$n_eval[i]), wf(X$fld_Hc_s[i],X$n_eval[i]),
 wf(X$ij2_H[i],X$n_eval[i]), wf(X$ij2_Hc[i],X$n_eval[i]), wf(X$joint_s[i],X$n_eval[i])))

cat("\n### Classification and bound location\n\n")
cat("| cell | sens | spec | PPV | NPV | mean \\|Hhat\\| | med bound | med theta(Hhat) | share >= 1.00 | share >= 1.25 |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(X))) cat(sprintf("| %s | %.4f | %.4f | %.4f | %.4f | %.1f | %.4f | %.4f | %s | %s |\n",
 cl(i), X$sens[i], X$spec[i], X$ppv[i], X$npv[i], X$mean_Hhat[i], X$med_bound[i], X$med_theta[i],
 wf(X$share_ge100[i],X$n_eval[i]), wf(X$share_ge125[i],X$n_eval[i])))

if (any(isnull)) {
cat("\n### HR 1.00 cells: selection rate and bound location, FS beside GRF, criterion named\n\n")
cat("Two separate quantities: the selection rate, and where the field lower bound sits. Rows over detected replicates with bound and theta(Hhat) finite (fs_extraction.R's row set).\n\n")
cat("| cell | engine | criterion | matched? | selection rate [Wilson] | n_eval | med lower bound | med theta(Hhat) | bound - theta | bound / theta | paired ratio | planted | **share >= 1.00** | **share >= 1.25** |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in which(isnull)) {
  mt <- identical(X$fs_focus[i],"effMaxSG") && isTRUE(all.equal(X$fs_eps[i],0.20))
  cat(sprintf("| %s | GRF | grfmr / effMaxSG / eps 0.2 | -- | %s | %d | %.4f | %.4f | %+.4f | %.4f | %.4f | %.4f | %s | %s |\n",
    cl(i), wf(X$detection[i],2000L), X$L_n[i], X$L_med_bound[i], X$L_med_theta[i],
    X$L_med_bound[i]-X$L_med_theta[i], X$L_med_bound[i]/X$L_med_theta[i], X$L_paired[i], X$planted[i],
    wf(X$L_s100[i],X$L_n[i]), wf(X$L_s125[i],X$L_n[i])))
  cat(sprintf("| %s | FS | %s / %s / eps %s | %s | %s | %d | %.4f | %.4f | %+.4f | %.4f | %.4f | %.4f | %s | %s |\n",
    cl(i), X$fs_camp[i], X$fs_focus[i], format(X$fs_eps[i]), if (mt) "**yes**" else "no",
    wf(X$fs_det[i],2000L), X$fsL_n[i], X$fsL_med_bound[i], X$fsL_med_theta[i],
    X$fsL_med_bound[i]-X$fsL_med_theta[i], X$fsL_med_bound[i]/X$fsL_med_theta[i], X$fsL_paired[i], X$planted[i],
    wf(X$fsL_s100[i],X$fsL_n[i]), wf(X$fsL_s125[i],X$fsL_n[i])))
}
cat("\nConfound, with every comparison: FS and GRF differ in identifier, in family construction and\n")
cat("in detection set; at 12.4% they differ in the selection criterion as well (maxeffCons eps 0.10\n")
cat("against effMaxSG eps 0.20), so a 12.4% gap cannot be read as engine behaviour even in part.\n")
}

cat("\n### Beside the FS comparator, with its criterion named\n\n")
cat("| cell | matched? | FS comparator | GRF / FS detection or selection | GRF field / FS field | GRF IJ2 / FS IJ2 | GRF pool / FS family (med) |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(X))) {
  mt <- identical(X$fs_focus[i],"effMaxSG") && isTRUE(all.equal(X$fs_eps[i],0.20))
  cat(sprintf("| %s | %s | %s / %s / eps %s | %.4f / %.4f | %.4f / %.4f | %.4f / %.4f | %g / %g |\n",
    cl(i), if (mt) "**yes**" else "no", X$fs_camp[i], X$fs_focus[i], format(X$fs_eps[i]),
    X$detection[i], X$fs_det[i], X$fld_H[i], X$fs_fld_H[i], X$ij2_H[i], X$fs_ij2_H[i],
    X$K_med[i], X$fs_K_med[i])) }
cat("\nConfound, with every comparison: FS and GRF differ in identifier, in family construction and\n")
cat("in detection set; at 12.4% they differ in the selection criterion as well (maxeffCons eps 0.10\n")
cat("against effMaxSG eps 0.20), so a 12.4% gap cannot be read as engine behaviour even in part.\n")
saveRDS(X, "grfmr_numbers.rds")
