# O1 readout -- columns matched to REPORT_stopB_md_harm_grid.md section 1 so the
# maxeffCons grid can be read directly beside the hr grid.
# All quantities ORIENTED (positive = harm) via the 1830ca92 bridge.
o1_row <- function(path) {
  b <- readRDS(path); r <- b$results; m <- b$meta
  det <- r$detected %in% 1L
  rd  <- r[which(det), , drop = FALSE]
  .orient <- if (isTRUE(m$adverse_outcome)) 1 else -1
  bH <- .orient * rd$betaHhat_H
  nb <- rd$nv_H_est - bH
  mb <- rd$mr_H_est - bH
  se <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
  lo <- rd$mr_H_lo; hi <- rd$mr_H_hi
  ok <- is.finite(bH) & is.finite(lo) & is.finite(hi)
  cov <- mean(bH[ok] >= lo[ok] & bH[ok] <= hi[ok])
  nsel <- rd$n_harm
  qs <- stats::quantile(nsel, c(.1,.25,.5,.75,.9), na.rm = TRUE)
  data.frame(
    n = m$cell_n, s = nrow(r),
    detected = sum(r$detected %in% 1L), mr_ok = sum(r$mr_ok %in% 1L),
    naive_bias = mean(nb, na.rm = TRUE),
    mr_bias    = mean(mb, na.rm = TRUE),
    mr_se      = se(mb),
    mr_bias_over_se = mean(mb, na.rm = TRUE) / se(mb),
    coverage = cov, mc_se = sqrt(cov * (1 - cov) / sum(ok)),
    width = mean(hi - lo, na.rm = TRUE),
    sens = mean(rd$sens, na.rm = TRUE), spec = mean(rd$spec, na.rm = TRUE),
    ppv  = mean(rd$ppv,  na.rm = TRUE), npv  = mean(rd$npv,  na.rm = TRUE),
    nsel_q10 = qs[[1]], nsel_q25 = qs[[2]], nsel_med = qs[[3]],
    nsel_q75 = qs[[4]], nsel_q90 = qs[[5]],
    secs_per_rep_med = stats::median(r$fit_mr_secs, na.rm = TRUE),
    secs_per_rep_mean = mean(r$fit_mr_secs, na.rm = TRUE),
    cell_elapsed_s = m$elapsed_sec,
    row.names = NULL)
}
files <- Sys.glob("o1_maxeffCons_mr_grid/o1_maxeffCons_mr_n*_s1000.rds")
files <- files[order(as.integer(sub(".*_n([0-9]+)_.*", "\\1", files)))]
tab <- do.call(rbind, lapply(files, o1_row))
options(width = 250)
cat("=== O1: maxeffCons MR-only grid ===\n"); print(tab, digits = 8, row.names = FALSE)
cat("\n=== trend columns, hr-grid layout ===\n")
print(tab[, c("n","s","detected","mr_ok","naive_bias","mr_bias","mr_se",
              "mr_bias_over_se","coverage","mc_se","width")], digits = 8, row.names = FALSE)
