suppressMessages(library(forestsearch))
options(width = 250)
R <- "results/"
`%||%` <- function(a,b) if (is.null(a)||length(a)==0||all(is.na(a))) b else a
cells <- list(c(1.50,500),c(1.50,1000),c(1.50,1500),c(1.75,500),c(1.75,1000),c(1.75,1500))
lab <- function(x) sprintf("A: HR %.2f, n %d", x[1], as.integer(x[2]))
file <- function(x) sprintf("%sdina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_nb20_dinamr_combined_1_2000.rds",
                            R, round(100*x[1]), as.integer(x[2]))
fsfile <- function(x) { camp <- if (abs(x[1]-1.75)<1e-9) "tier2" else "p12ext"
  sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*x[1]), as.integer(x[2]), camp) }
B <- lapply(cells, function(x) readRDS(file(x))); names(B) <- vapply(cells, lab, "")
labs <- names(B); NN <- vapply(cells, function(x) as.integer(x[2]), 0L); HR <- vapply(cells, function(x) x[1], 0)

# ---- the document's own cov machinery (summary_dinamr.qmd chunk `cov-fns`) ----
subst_s <- function(r) { for (s in c("est2","up1s","lo1s","lo2s","hi2s","lo_se","hi_se","se","lam_mean"))
  r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
err_sd <- function(r, block, est) { d <- r[r$detected %in% 1L, ]
  e <- switch(est, naive=d[[paste0("nv_",block,"_est")]], mr=d[[paste0("mr_",block,"_est")]], fld=d[[paste0("fld_",block,"_est2")]])
  ok <- is.finite(e) & is.finite(d[[paste0("betaHhat_",block)]])
  stats::sd(log(e[ok]) - log(d[[paste0("betaHhat_",block)]][ok])) }
cov_block <- function(r, block, side) {
  t1 <- fs_sim_bias_coverage(r, block=block, estimators=c("naive","mr","fld"), side=side)
  t1$sd_err <- vapply(as.character(t1$estimator), function(e) err_sd(r, block, e), 0)
  t1$construction <- c(naive="naive", mr="IJ two-term", fld="field")[as.character(t1$estimator)]
  if (block=="Hc") { rs <- subst_s(r); t2 <- fs_sim_bias_coverage(rs, block="Hc", estimators="fld", side=side)
    t2$sd_err <- err_sd(rs,"Hc","fld"); t2$construction <- "field-s"; t1 <- rbind(t1,t2) }
  t1$se_over_sd_err <- t1$se_mean/t1$sd_err; t1$b_err <- t1$bias_log/t1$sd_err; t1 }

cat("\n##################################################################################\n")
cat("### 1. STANDARD TABLE, BLOCK A (12.4%), per cell, both subgroups, every product\n")
cat("###    Coverage = CONDITIONAL-ON-PROPOSED-FAMILY estimand, detected replicates only.\n")
cat("##################################################################################\n")
for (i in seq_along(labs)) { k <- labs[i]; r <- B[[k]]$results
  det <- mean(r$detected %in% 1L)
  cat(sprintf("\n--- %s   [detection %.4f = %d/2000] ---\n", k, det, sum(r$detected %in% 1L)))
  tab <- rbind(cbind(subgroup="Hhat (lower)",   cov_block(r,"H","lower")),
               cbind(subgroup="Hhat^c (upper)", cov_block(r,"Hc","upper")))
  tab$construction <- factor(tab$construction, levels=c("naive","field","field-s","IJ two-term"))
  tab <- tab[order(tab$subgroup, tab$construction), ]
  print(tab[, c("subgroup","construction","n","bias_log","sd_emp","sd_err","b","b_err","se_mean","r",
                "se_over_sd_err","cov1","cov1_wilson_lo","cov1_wilson_hi","cov2","cov1_ref")],
        row.names=FALSE, digits=4) }
