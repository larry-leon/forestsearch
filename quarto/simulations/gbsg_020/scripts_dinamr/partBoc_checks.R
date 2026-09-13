# Part B OC smoke (TASK_partB_enabling_2026-09-12): two readouts beside the OC
# table.  Reads committed bundles only.
#  1. The smoke (MR OFF) against the committed MR-ON bundle at the SAME cell and
#     seeds, sim_id 1-30: identification and classification should be identical
#     (MR cannot change the subgroup; forestsearch_main.R:1041-1044).  DINA and
#     GRF effMaxSG eps 0.20 against dinamr / grfmr (Mac); FS maxeffCons against
#     p12ext (pop-os, a cross-machine comparison).
#  2. A cell-profile adjustment of the flat one-cell projection, using each
#     engine's committed 18-cell profile of the field-excluded bound
#     (partB_stage0_readout.rds): multiplier = mean over 18 cells / value at
#     this cell.  An approximation, labelled as such -- the bound still carries
#     the non-field MR cost, and FS's grid is mostly pop-os.
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
cmp <- list(
  c("dina", "dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_nomr_pBoc_res_1_30.rds",
    "dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_dinamr_res_1_1000.rds"),
  c("grf", "grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_nomr_pBoc_res_1_30.rds",
    "grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_grfmr_res_1_1000.rds"),
  c("consistency", "fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n500_nomr_pBoc_res_1_30.rds",
    "fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n500_p12ext_res_1_1000.rds"))
cols <- c("detected", "status", "sg_def", "n_sel", "n_harm", "n_true", "label",
          "sens", "spec", "ppv", "npv", "betaHhat_H", "betaHhat_Hc", "admitted_n")
cat("=== 1. SMOKE (MR off) vs COMMITTED (MR on), sim_id 1-30, same cell and seeds ===\n")
for (x in cmp) {
  a <- readRDS(file.path(RES, x[2])); b <- readRDS(file.path(RES, x[3]))
  ra <- a$results; rb <- b$results[match(ra$sim_id, b$results$sim_id), ]
  cc <- intersect(cols, names(rb))
  ok <- vapply(cc, function(k) identical(ra[[k]], rb[[k]]), logical(1))
  eqv <- vapply(cc, function(k) isTRUE(all.equal(ra[[k]], rb[[k]], tolerance = 1e-12, check.attributes = FALSE)), logical(1))
  cat(sprintf("\n%s vs %s (%s, %s, R %s):\n", x[1], sub("_res_1_1000.rds", "", x[3]), b$meta$hostname, b$meta$forestsearch_version, b$meta$r_version))
  cat(sprintf("  identical(): %s\n", paste(sprintf("%s=%s", cc, ok), collapse = " ")))
  if (any(!eqv)) cat(sprintf("  NOT equal within 1e-12: %s\n", paste(cc[!eqv], collapse = ", ")))
  cat(sprintf("  truth identical: %s ; all.equal 1e-12: %s\n", identical(a$truth, b$truth),
              isTRUE(all.equal(a$truth, b$truth, tolerance = 1e-12))))
  for (k in cc[!ok]) { d <- which(!mapply(identical, ra[[k]], rb[[k]]))
    cat(sprintf("  %s differs on sim_id %s (max abs diff %s)\n", k, paste(ra$sim_id[d], collapse = ","),
        if (is.numeric(ra[[k]])) format(max(abs(ra[[k]][d] - rb[[k]][d]), na.rm = TRUE)) else "n/a")) }
}
cat("\n=== the replicate only maxeff detects (consistency) ===\n")
me <- readRDS(file.path(RES, "fs_maxeff_fb_mr_field_m1_h150_knoise0_n500_nomr_pBoc_res_1_30.rds"))$results
mc <- readRDS(file.path(RES, "fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n500_nomr_pBoc_res_1_30.rds"))$results
w <- which(me$detected %in% 1L & !(mc$detected %in% 1L))
print(me[w, c("sim_id", "sg_def", "n_sel", "n_true", "sens", "ppv", "betaHhat_H")])
d2 <- which(me$detected %in% 1L & mc$detected %in% 1L & me$sg_def != mc$sg_def)
cat("maxeff vs maxeffCons, jointly detected but different rules:\n")
print(data.frame(sim_id = me$sim_id[d2], maxeff = me$sg_def[d2], maxeffCons = mc$sg_def[d2]))

cat("\n=== 2. CELL-PROFILE-ADJUSTED PROJECTION (approximation) ===\n")
pc <- readRDS("partB_stage0_readout.rds")$per_cell
oc <- readRDS("partBoc_table.rds")
mult <- sapply(c("consistency", "dina", "grf"), function(e) {
  s <- pc[pc$engine == e, ]
  here <- s$bnd_worker_h[s$prev == "12.4%" & s$hr == 1.50 & s$n == 500]
  mean(s$bnd_worker_h) / here })
print(round(mult, 3))
PJ <- oc$projection
PJ$mult <- unname(mult[PJ$engine])
PJ$compute_adj_h <- round(PJ$compute_h * PJ$mult, 2)
PJ$wall_adj_h <- round(PJ$compute_adj_h + PJ$overhead_h, 2)
print(PJ, row.names = FALSE)
TOT <- aggregate(cbind(wall_h, wall_adj_h) ~ reps, PJ, sum)
print(TOT[order(-TOT$reps), ], row.names = FALSE)
saveRDS(list(multiplier = mult, projection_adj = PJ, total_adj = TOT), "partBoc_checks.rds")
