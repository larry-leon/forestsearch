# Stage 0d: single-replicate timing of the continuous twin at committed settings (ci_method = "ij").
# Evaluates the twin's setup-knobs / build-dgm / machinery chunks verbatim, then times record_replicate(1).
qmd <- "sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd"
txt <- readLines(qmd)
chunk <- function(label) {
  s <- grep(sprintf("^```\\{r %s[,}]", label), txt)[1]
  e <- s + which(grepl("^```\\s*$", txt[(s+1):length(txt)]))[1]
  txt[(s+1):(e-1)]
}
t_setup <- system.time(eval(parse(text = chunk("setup-knobs")), envir = globalenv()))
t_dgm   <- system.time(eval(parse(text = chunk("build-dgm")),   envir = globalenv()))
eval(parse(text = chunk("machinery")), envir = globalenv())
cat(sprintf("setup %.1f s | DGM build %.1f s | cores=%d n_workers(default)=%d\n",
            t_setup[["elapsed"]], t_dgm[["elapsed"]], parallel::detectCores(logical = FALSE), n_workers))
t1 <- system.time(r1 <- record_replicate(1L))
cat(sprintf("record_replicate(1): wall %.1f s | fit_mr_secs %.1f | status %s | n_harm %d | nv_H_est %.10f | mr_H_est %.10f | mr_H_se_ij %.10f | sg_def %s\n",
            t1[["elapsed"]], r1$fit_mr_secs, r1$status, r1$n_harm, r1$nv_H_est, r1$mr_H_est, r1$mr_H_se_ij, r1$sg_def))
b <- readRDS("mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds")$results
b1 <- b[b$sim_id == 1L, ]
num <- names(r1)[vapply(r1, is.numeric, logical(1))]
num <- setdiff(num, c("fit_mr_secs", "fb_secs"))
rel <- sapply(num, function(k) { a <- r1[[k]]; z <- b1[[k]]
  if (is.na(a) && is.na(z)) 0 else if (is.na(a) || is.na(z)) NA else abs(a - z) / max(abs(z), 1e-300) })
cat("committed row 1 sg_def:", b1$sg_def, "\n")
cat("max rel diff over numeric pre-existing columns (excl. timing):", format(max(rel, na.rm = TRUE), digits = 3), "\n")
print(round(rel[rel > 1e-8 | is.na(rel)], 12))
saveRDS(list(r1 = r1, t1 = t1, t_setup = t_setup, t_dgm = t_dgm), file.path(Sys.getenv("SP"), "stage0_rep1_ij.rds"))
