# Part B readout: the GRF cost surface from the five 36-replicate probes.
# TASK_dinamr_blockC_grfprobe_2026-09-11.  COST AND MECHANISM ONLY.
#
# NOT reported here, deliberately: any coverage table, any FS or DINA
# comparison, any acceptance criterion, any recommendation.  GRF is NOT
# "FS-analogous" -- a GRF-to-FS or GRF-to-DINA comparison differs in
# identifier, family construction, detection set, selection criterion AND the
# scale of the selection criterion (Larry, 2026-09-11).
#
# Gate 3 (alignment) is asserted here from the resolved template values.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")
TPL     <- file.path(QMD_DIR, "sim_fs_maxeffCons_fb_mr_field_m1_template.qmd")

## ---- GATE 3: alignment, STOP on failure ------------------------------------
# grf_selection and grf_select_statistic are template literals with no
# environment override anywhere in the document, so they resolve to whatever
# the source says.  Read them back from the source rather than asserting them.
L   <- readLines(TPL)
gv  <- function(nm) {
  i <- grep(sprintf("^%s\\s*<-", nm), L)
  if (!length(i)) return(NA_character_)
  v <- sub(sprintf("^%s\\s*<-\\s*", nm), "", L[i[1]])
  list(value = trimws(sub("#.*$", "", v)), line = i[1])
}
gs <- gv("grf_selection"); gt <- gv("grf_select_statistic"); gd <- gv("dmin.grf")
cat("===== GATE 3 -- ALIGNMENT =====\n")
cat(sprintf("  grf_selection        resolves to %s   (template line %d)\n", gs$value, gs$line))
cat(sprintf("  grf_select_statistic resolves to %s   (template line %d)\n", gt$value, gt$line))
cat(sprintf("  dmin.grf             resolves to %s   (template line %d)\n", gd$value, gd$line))
ov <- grep("FS_S7_GRF|env_chr\\(\"FS_S7_.*GRF|grf_selection\\s*<-\\s*\\.env|grf_select_statistic\\s*<-\\s*\\.env", L)
cat(sprintf("  environment overrides for either knob: %s\n",
            if (length(ov)) paste("FOUND at lines", paste(ov, collapse = ",")) else "none (both are literals)"))
g3 <- identical(gs$value, "\"frontier\"") && identical(gt$value, "\"effect\"")
cat(sprintf("  GATE 3: %s\n", if (g3) "PASS" else "**FAIL -- STOP**"))
if (!g3) stop("Gate 3 alignment failed; Part B does not proceed.")

## ---- the five probes -------------------------------------------------------
wilson <- function(x, n, z = 1.959964) {
  if (!n) return(c(NA, NA)); p <- x/n; d <- 1 + z^2/n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d,
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d)
}
corners <- list(
  list(hr=1.50, n= 500L, z=FALSE, lab="g_p124_h150_n500"),
  list(hr=1.50, n=1500L, z=FALSE, lab="g_p124_h150_n1500"),
  list(hr=1.50, n= 500L, z=TRUE,  lab="g_p31_h150_n500"),
  list(hr=1.50, n=1500L, z=TRUE,  lab="g_p31_h150_n1500"),
  list(hr=1.00, n= 500L, z=FALSE, lab="g_p124_h100_n500"))
pf <- function(cc) file.path(RES, sprintf(
  "grf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfprobe_res_1_36.rds",
  round(100*cc$hr), cc$n, if (cc$z) "_z1q60" else ""))

recov <- c("fld_recov_sens_H","fld_recov_ppv_H","fld_recov_sens_Hc","fld_recov_npv_Hc",
           "fld_recov_q10","fld_recov_q50","fld_recov_q90","fld_recov_share1","fld_recov_n_used")
peak_mb <- function(lab) {
  f <- file.path(SCRATCH, "logs", sprintf("rss_%s.txt", lab))
  g <- file.path(SCRATCH, "logs", sprintf("time_%s.txt", lab))
  s <- if (file.exists(f)) { v <- read.table(f, header = FALSE); max(v$V2)/1024 } else NA_real_
  t <- if (file.exists(g)) {
    ln <- grep("maximum resident set size", readLines(g, warn = FALSE), value = TRUE)
    if (length(ln)) max(as.numeric(sub("^\\s*([0-9]+).*", "\\1", ln)))/1024/1024 else NA_real_
  } else NA_real_
  c(sampler_tree_MB = s, time_l_maxrss_MB = t)
}

cat("\n\n===== PER-PROBE COST RECORD (36 replicates, 12 workers) =====\n")
ALL <- list()
for (cc in corners) {
  f <- pf(cc)
  cat(sprintf("\n---------- %s : HR %.2f, n %d, prevalence %s ----------\n",
              cc$lab, cc$hr, cc$n, if (cc$z) "31%" else "12.4%"))
  if (!file.exists(f)) { cat("  NOT ON DISK:", basename(f), "\n"); next }
  b <- readRDS(f); r <- b$results
  s <- r$fit_mr_secs; K <- r$n_family
  det <- sum(r$detected %in% 1L); N <- nrow(r); w <- wilson(det, N)
  cat(sprintf("  wall per replicate (s)      : median %.2f   p90 %.2f   max %.2f   (mean %.2f, total %.0f)\n",
      median(s, na.rm=TRUE), quantile(s, .9, names=FALSE, na.rm=TRUE),
      max(s, na.rm=TRUE), mean(s, na.rm=TRUE), sum(s, na.rm=TRUE)))
  Kf <- K[is.finite(K)]
  if (length(Kf)) {
    cat(sprintf("  proposed-family size        : min %g  q10 %g  q25 %g  med %g  q75 %g  q90 %g  max %g\n",
        min(Kf), quantile(Kf,.1,names=FALSE), quantile(Kf,.25,names=FALSE), median(Kf),
        quantile(Kf,.75,names=FALSE), quantile(Kf,.9,names=FALSE), max(Kf)))
    ok <- is.finite(K) & is.finite(s)
    cat(sprintf("  wall against family size    : Spearman rho %.3f over %d replicates; ",
        suppressWarnings(cor(K[ok], s[ok], method = "spearman")), sum(ok)))
    if (sum(ok) >= 4) {
      q <- cut(K[ok], breaks = unique(quantile(K[ok], c(0,.5,1))), include.lowest = TRUE)
      cat(sprintf("median wall by family half: %s\n",
          paste(sprintf("%s -> %.2fs", levels(q), tapply(s[ok], q, median)), collapse = " | ")))
    } else cat("\n")
  } else cat("  proposed-family size        : no finite n_family (all replicates undetected)\n")
  cat(sprintf("  detection                   : %d / %d = %.4f  [Wilson 95%%: %.4f, %.4f]\n",
      det, N, det/N, w[1], w[2]))
  pm <- peak_mb(cc$lab)
  cat(sprintf("  peak memory                 : process-tree sampler %.0f MB ; /usr/bin/time -l max RSS %.0f MB\n",
      pm[1], pm[2]))
  D <- r[r$detected %in% 1L, , drop = FALSE]
  pres <- function(k) k %in% names(r)
  popd <- function(k) pres(k) && nrow(D) && !all(is.na(D[[k]]))
  cat(sprintf("  p-hat   present/populated   : %s / %s\n",
      all(sapply(c("p_hat_H","p_hat_sum","p_hat_top1"), pres)),
      all(sapply(c("p_hat_H","p_hat_sum","p_hat_top1"), popd))))
  cat(sprintf("  rho-c   present/populated   : %s / %s\n",
      pres("fld_Hc_scale_ratio"), popd("fld_Hc_scale_ratio")))
  cat(sprintf("  nine recovery columns       : %d of 9 present, %d of 9 populated\n",
      sum(sapply(recov, pres)), sum(sapply(recov, popd))))
  miss <- recov[!sapply(recov, popd)]
  if (length(miss)) cat("     not populated:", paste(miss, collapse = ", "), "\n")
  ALL[[cc$lab]] <- data.frame(lab = cc$lab, hr = cc$hr, n = cc$n,
    prev = if (cc$z) "31%" else "12.4%", reps = N, detection = det/N,
    wilson_lo = w[1], wilson_hi = w[2],
    s_med = median(s, na.rm=TRUE), s_q90 = quantile(s,.9,names=FALSE,na.rm=TRUE),
    s_max = max(s, na.rm=TRUE), s_total = sum(s, na.rm=TRUE),
    K_med = if (length(Kf)) median(Kf) else NA, K_q90 = if (length(Kf)) quantile(Kf,.9,names=FALSE) else NA,
    K_max = if (length(Kf)) max(Kf) else NA,
    peak_tree_MB = pm[1], peak_time_MB = pm[2])
}
if (length(ALL)) {
  A <- do.call(rbind, ALL)
  cat("\n\n===== COST SURFACE =====\n")
  print(A[, c("prev","hr","n","reps","detection","s_med","s_q90","s_max","K_med","K_q90","K_max","peak_tree_MB")],
        row.names = FALSE, digits = 4)
  cat("\nNo coverage table, no FS/DINA comparison, no acceptance criterion, no recommendation.\n")
  saveRDS(A, file.path(SCRATCH, "grfprobe_surface.rds"))
}
