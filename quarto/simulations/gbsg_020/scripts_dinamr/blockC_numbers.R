# Block C numbers for REPORT_dinamr_blockC_2026-09-11, extracted with
# summary_dinamr.qmd's OWN definitions (cov-fns, strat-fns, strat-xtab -- verbatim in
# form) from the committed bundles.  Reads committed artifacts only; runs no simulation.
# Block C numbers, extracted with summary_dinamr.qmd's OWN definitions (verbatim
# in form) from the committed bundles.  Reading committed artifacts only.
suppressPackageStartupMessages(library(forestsearch))
setwd("/Users/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020")
R <- "results/"
wil <- function(p, n, z = stats::qnorm(0.975)) {
  ctr <- (p + z^2/(2*n))/(1 + z^2/n)
  hw  <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1 + z^2/n); c(ctr-hw, ctr+hw) }
subst_s <- function(r) { for (s in c("est2","up1s","lo1s","lo2s","hi2s","lo_se","hi_se","se","lam_mean"))
  r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
err_sd <- function(r, block, est) { d <- r[r$detected %in% 1L, ]
  e <- switch(est, naive = d[[paste0("nv_", block, "_est")]],
                   mr    = d[[paste0("mr_", block, "_est")]],
                   fld   = d[[paste0("fld_", block, "_est2")]])
  ok <- is.finite(e) & is.finite(d[[paste0("betaHhat_", block)]])
  stats::sd(log(e[ok]) - log(d[[paste0("betaHhat_", block)]][ok])) }
cov_block <- function(r, block, side) {
  t1 <- fs_sim_bias_coverage(r, block = block, estimators = c("naive","mr","fld"), side = side)
  t1$sd_err <- vapply(as.character(t1$estimator), function(e) err_sd(r, block, e), 0)
  t1$construction <- c(naive="naive", mr="IJ two-term", fld="field")[as.character(t1$estimator)]
  if (block == "Hc") { rs <- subst_s(r)
    t2 <- fs_sim_bias_coverage(rs, block = "Hc", estimators = "fld", side = side)
    t2$sd_err <- err_sd(rs, "Hc", "fld"); t2$construction <- "field-s"; t1 <- rbind(t1, t2) }
  t1$se_over_sd_err <- t1$se_mean / t1$sd_err
  t1$b_err <- t1$bias_log / t1$sd_err
  z95_ <- stats::qnorm(0.95)
  t1$cov1_ref_err <- if (identical(side,"lower")) stats::pnorm(z95_*t1$se_over_sd_err - t1$b_err) else
                                                  stats::pnorm(z95_*t1$se_over_sd_err + t1$b_err)
  t1 }
tert_fs <- function(v) { if (!any(is.finite(v))) return(rep(NA_integer_, length(v)))
  br <- stats::quantile(v, c(0,1/3,2/3,1), na.rm=TRUE, names=FALSE)
  if (length(unique(br)) < 2L) return(rep(1L, length(v)))
  cut(v, breaks=unique(br), include.lowest=TRUE, labels=FALSE) }

cells <- list(list("A",500L),list("A",1000L),list("A",1500L),
              list("B",500L),list("B",1000L),list("B",1500L))
lab <- function(b,n) sprintf("%s n%d", if (b=="A") "12.4%" else "31%", n)
dfile <- function(b,n) sprintf("%sdina_effMaxSG_fb_mr_field_m1_h100_knoise0_n%d%s_nb20_dinamr_combined_1_2000.rds",
                               R, n, if (b=="A") "" else "_z1q60")
ffile <- function(b,n) {
  if (b == "A") sprintf("%sfs_maxeffCons_fb_mr_field_m1_h100_knoise0_n%d_%s_combined_1_2000.rds",
                        R, n, if (n == 500L) "tier2" else "p12ext")
  else sprintf("%sfs_effMaxSG_fb_mr_field_m1_h100_knoise0_n%d_z1q60_nb20_cert20_combined_1_2000.rds", R, n) }

cat("################ 1. STANDARD CONSTRUCTIONS, BLOCK C ################\n")
cat("bias_log = retained bias on the log scale; sd_emp = marginal SD; sd_err = error SD;\n")
cat("se_mean = mean SE; cov1 = one-sided 95% on the exposed side with Wilson limits;\n")
cat("cov2 = two-sided; cov1_ref / cov1_ref_err = Gaussian reference on marginal / error SD.\n\n")
for (x in cells) { b <- x[[1]]; n <- x[[2]]; r <- readRDS(dfile(b,n))$results
  t <- rbind(cbind(block="Hhat (lower)",   cov_block(r,"H","lower")),
             cbind(block="Hhat^c (upper)", cov_block(r,"Hc","upper")))
  t$construction <- factor(t$construction, levels=c("naive","field","field-s","IJ two-term"))
  t <- t[order(t$block, t$construction), ]
  cat(sprintf("---- %s ----\n", lab(b,n)))
  print(t[,c("block","construction","n","bias_log","sd_emp","sd_err","b","b_err","se_mean","r",
             "se_over_sd_err","cov1","cov1_wilson_lo","cov1_wilson_hi","cov2","cov1_ref","cov1_ref_err")],
        row.names=FALSE, digits=4); cat("\n") }

det_rows2 <- function(r) { d <- r[r$detected %in% 1L, , drop=FALSE]
  keep <- is.finite(d$betaHhat_H) & is.finite(d$betaHhat_Hc) & is.finite(d$n_family) &
          is.finite(d$p_hat_H) & is.finite(d$fld_H_lo1s) & is.finite(d$fld_Hc_up1s_s)
  d[keep, , drop=FALSE] }
cvw <- function(p,n) { w <- wil(p,n); sprintf("%.4f [%.4f, %.4f]", p, w[1], w[2]) }

cat("\n\n################ 2. JOINT, IJ MISS SPLIT, NAIVE -- BLOCK C ################\n\n")
J <- do.call(rbind, lapply(cells, function(x){ b<-x[[1]]; n<-x[[2]]; d<-det_rows2(readRDS(dfile(b,n))$results)
  N <- nrow(d)
  data.frame(cell=lab(b,n), n_eval=N,
    joint_bonf   = cvw(mean(d$betaHhat_H>=d$fld_joint_bonf_loH   & d$betaHhat_Hc<=d$fld_joint_bonf_upHc), N),
    joint_bonf_s = cvw(mean(d$betaHhat_H>=d$fld_joint_s_bonf_loH & d$betaHhat_Hc<=d$fld_joint_s_bonf_upHc), N),
    naive_H_2s   = cvw(mean(d$betaHhat_H>=d$nv_H_lo & d$betaHhat_H<=d$nv_H_hi), N),
    naive_Hc_2s  = cvw(mean(d$betaHhat_Hc>=d$nv_Hc_lo & d$betaHhat_Hc<=d$nv_Hc_hi), N)) }))
print(J, row.names=FALSE)
cat("\nIJ two-term two-sided, miss split by side:\n")
M <- do.call(rbind, lapply(cells, function(x){ b<-x[[1]]; n<-x[[2]]; d<-det_rows2(readRDS(dfile(b,n))$results)
  N<-nrow(d)
  data.frame(cell=lab(b,n), n_eval=N,
    ij2_H=cvw(mean(d$betaHhat_H>=d$mr_H_lo & d$betaHhat_H<=d$mr_H_hi),N),
    miss_below_H=mean(d$betaHhat_H<d$mr_H_lo), miss_above_H=mean(d$betaHhat_H>d$mr_H_hi),
    ij2_Hc=cvw(mean(d$betaHhat_Hc>=d$mr_Hc_lo & d$betaHhat_Hc<=d$mr_Hc_hi),N),
    miss_below_Hc=mean(d$betaHhat_Hc<d$mr_Hc_lo), miss_above_Hc=mean(d$betaHhat_Hc>d$mr_Hc_hi)) }))
print(M, row.names=FALSE, digits=4)

cat("\n\n################ 3. FIELD LOWER BY n_family TERTILE AND p-hat BIN ################\n")
cat("coverage of beta(Hhat) by the field lower bound, Wilson; retained bias = mean(log(fld_H_est2) - log(betaHhat_H))\n\n")
for (by in c("K","p")) {
  cat(sprintf("=== stratified by %s ===\n", if (by=="K") "n_family (K tertiles, plus overlapping K=1 and K<=5)" else "p-hat tertiles"))
  S <- do.call(rbind, lapply(cells, function(x){ b<-x[[1]]; n<-x[[2]]; d<-det_rows2(readRDS(dfile(b,n))$results)
    v <- if (by=="K") d$n_family else d$p_hat_H; g <- tert_fs(v)
    rows <- lapply(sort(unique(g[!is.na(g)])), function(t){ dk <- d[!is.na(g)&g==t,,drop=FALSE]
      data.frame(cell=lab(b,n), stratum=sprintf("%s T%d [%s, %s]", if(by=="K")"K" else "p-hat", t,
                 format(min(v[!is.na(g)&g==t])), format(max(v[!is.na(g)&g==t]))), overlapping=FALSE,
                 n=nrow(dk), field_cov=cvw(mean(dk$betaHhat_H>=dk$fld_H_lo1s), nrow(dk)),
                 retained_bias=mean(log(dk$fld_H_est2)-log(dk$betaHhat_H))) })
    if (by=="K") for (cut in c(1L,5L)) { dk <- d[if(cut==1L) d$n_family==1L else d$n_family<=5L,,drop=FALSE]
      if (nrow(dk)) rows[[length(rows)+1]] <- data.frame(cell=lab(b,n),
        stratum=sprintf("K %s %d (overlapping)", if(cut==1L) "=" else "<=", cut), overlapping=TRUE,
        n=nrow(dk), field_cov=cvw(mean(dk$betaHhat_H>=dk$fld_H_lo1s), nrow(dk)),
        retained_bias=mean(log(dk$fld_H_est2)-log(dk$betaHhat_H))) }
    rows[[length(rows)+1]] <- data.frame(cell=lab(b,n), stratum="all detected (overlapping)", overlapping=TRUE,
      n=nrow(d), field_cov=cvw(mean(d$betaHhat_H>=d$fld_H_lo1s), nrow(d)),
      retained_bias=mean(log(d$fld_H_est2)-log(d$betaHhat_H)))
    do.call(rbind, rows) }))
  print(S, row.names=FALSE, digits=4); cat("\n") }

cat("\n################ 4. JOINT COUNT TABLE: K stratum x p-hat bin ################\n\n")
for (x in cells) { b<-x[[1]]; n<-x[[2]]; d<-det_rows2(readRDS(dfile(b,n))$results)
  gk<-tert_fs(d$n_family); gp<-tert_fs(d$p_hat_H)
  tb <- table(K=paste0("K T",gk), p_hat=paste0("p T",gp))
  cat(sprintf("---- %s (n_eval %d) ----\n", lab(b,n), nrow(d))); print(tb)
  cat(sprintf("  overlapping: K = 1 -> %d ; K <= 5 -> %d\n\n", sum(d$n_family==1L), sum(d$n_family<=5L))) }

cat("\n################ 5. FS COMPARATOR BESIDE DINA, BLOCK C ################\n\n")
FS <- do.call(rbind, lapply(cells, function(x){ b<-x[[1]]; n<-x[[2]]
  fb <- readRDS(ffile(b,n)); fr <- fb$results; d <- det_rows2(readRDS(dfile(b,n))$results); f <- det_rows2(fr)
  data.frame(cell=lab(b,n),
    criterion_matched = identical(fb$meta$sg_focus,"effMaxSG") && isTRUE(all.equal(fb$meta$effect_neighborhood,0.20)),
    fs_criterion = sprintf("%s / %s / eps %s", fb$meta$campaign_tag, fb$meta$sg_focus, format(fb$meta$effect_neighborhood)),
    dina_field = cvw(mean(d$betaHhat_H>=d$fld_H_lo1s), nrow(d)),
    fs_field   = cvw(mean(f$betaHhat_H>=f$fld_H_lo1s), nrow(f)),
    dina_fields_up = cvw(mean(d$betaHhat_Hc<=d$fld_Hc_up1s_s), nrow(d)),
    fs_fields_up   = cvw(mean(f$betaHhat_Hc<=f$fld_Hc_up1s_s), nrow(f)),
    dina_joint_s = cvw(mean(d$betaHhat_H>=d$fld_joint_s_bonf_loH & d$betaHhat_Hc<=d$fld_joint_s_bonf_upHc), nrow(d)),
    fs_joint_s   = cvw(mean(f$betaHhat_H>=f$fld_joint_s_bonf_loH & f$betaHhat_Hc<=f$fld_joint_s_bonf_upHc), nrow(f)),
    dina_ij2 = cvw(mean(d$betaHhat_H>=d$mr_H_lo & d$betaHhat_H<=d$mr_H_hi), nrow(d)),
    fs_ij2   = cvw(mean(f$betaHhat_H>=f$mr_H_lo & f$betaHhat_H<=f$mr_H_hi), nrow(f))) }))
print(FS[,c("cell","criterion_matched","fs_criterion","dina_field","fs_field")], row.names=FALSE)
cat("\n"); print(FS[,c("cell","dina_fields_up","fs_fields_up","dina_joint_s","fs_joint_s")], row.names=FALSE)
cat("\n"); print(FS[,c("cell","dina_ij2","fs_ij2")], row.names=FALSE)
