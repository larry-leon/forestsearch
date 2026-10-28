# Post-process, post-selection vignette stage 0 (2026-10-28). Run from this directory
# after run_ladder.sh. Reads rungs/*.rds, fits/*.rds, ladder.csv and the Table S2
# payload copy in manuscript_source/; writes comparison_table_s2.csv and
# budget.csv and prints everything the report quotes.
options(width = 140)
P <- readRDS("manuscript_source/analysis_gbsg_survival_multimethod_payload.rds")
L <- read.csv("ladder.csv", stringsAsFactors = FALSE)
rd <- function(r) { f <- file.path("rungs", paste0("rung_", r, ".rds")); if (file.exists(f)) readRDS(f) else NULL }

# ---- 1. M-FS vs M-FS-ij: g$debiased (and complement) all.equal -----------------
a <- readRDS("fits/fit_M-FS.rds")$mr_inference; b <- readRDS("fits/fit_M-FS-ij.rds")$mr_inference
pick <- function(g) list(H = g$debiased[c("est", "lower", "upper", "se_ij")],
                         Hc = g$complement$debiased[c("est", "lower", "upper", "se_ij")],
                         naive_H = g$naive[c("est", "lower", "upper")])
cat("==== all.equal(M-FS debiased, M-FS-ij debiased), tolerance 0 ====\n")
print(all.equal(pick(a), pick(b), tolerance = 0))
cat("identical():", identical(pick(a), pick(b)), "\n")
cat("ci_method:", a$ci_method, "vs", b$ci_method, "; field block present:", !is.null(a$field), "vs", !is.null(b$field), "\n")
cat(sprintf("MR multiplier-loop timing_seconds: %.3f (field run) vs %.3f (ij run)\n", a$timing_seconds, b$timing_seconds))

# ---- 2. Table S2 comparison --------------------------------------------------
ids <- c(FS = "M-FS", DINA = "M-DINA", GRF = "M-GRF"); bids <- c(FS = "B-FS", DINA = "B-DINA", GRF = "B-GRF")
rows <- list()
for (m in names(ids)) {
  q <- rd(ids[[m]])$mr; rr <- rd(ids[[m]])$row; bb <- rd(bids[[m]])
  for (reg in c("H", "Hc")) {
    pt <- P$table[P$table$method == m & P$table$region == reg, ]
    s <- if (reg == "H") "_H" else "_Hc"
    n_now <- if (reg == "H") rr$N_H else rr$N_Hc
    fbH <- if (!is.null(bb)) (if (reg == "H") bb$H_estimates else bb$Hc_estimates) else NULL
    rows[[length(rows) + 1L]] <- data.frame(
      method = m, region = reg, n_payload = pt$n, n_now = n_now,
      naive_payload = sprintf("%.6f (%.6f, %.6f)", pt$naive_est, pt$naive_lo, pt$naive_hi),
      naive_now = sprintf("%.6f (%.6f, %.6f)", q[[paste0("naive", s)]], q[[paste0("naive_lo", s)]], q[[paste0("naive_hi", s)]]),
      naive_2dp_match = identical(sprintf("%.2f", c(pt$naive_est, pt$naive_lo, pt$naive_hi)),
                                  sprintf("%.2f", c(q[[paste0("naive", s)]], q[[paste0("naive_lo", s)]], q[[paste0("naive_hi", s)]]))),
      naive_maxabs_diff = max(abs(c(pt$naive_est, pt$naive_lo, pt$naive_hi) -
                                  c(q[[paste0("naive", s)]], q[[paste0("naive_lo", s)]], q[[paste0("naive_hi", s)]]))),
      mr_payload = sprintf("%.6f (%.6f, %.6f)", pt$mr_est, pt$mr_lo, pt$mr_hi),
      mr_now = sprintf("%.6f (%.6f, %.6f)", q[[paste0("mr", s)]], q[[paste0("mr_lo", s)]], q[[paste0("mr_hi", s)]]),
      mr_2dp_match = identical(sprintf("%.2f", c(pt$mr_est, pt$mr_lo, pt$mr_hi)),
                               sprintf("%.2f", c(q[[paste0("mr", s)]], q[[paste0("mr_lo", s)]], q[[paste0("mr_hi", s)]]))),
      mr_maxabs_diff = max(abs(c(pt$mr_est, pt$mr_lo, pt$mr_hi) -
                               c(q[[paste0("mr", s)]], q[[paste0("mr_lo", s)]], q[[paste0("mr_hi", s)]]))),
      fb_payload_B1000 = sprintf("%.6f (%.6f, %.6f)", pt$fb_est, pt$fb_lo, pt$fb_hi),
      fb_now_B5 = if (is.null(fbH)) NA_character_ else sprintf("%.6f (%.6f, %.6f)", fbH$H2, fbH$H2_lower, fbH$H2_upper),
      stringsAsFactors = FALSE)
  }
}
cmp <- do.call(rbind, rows)
write.csv(cmp, "comparison_table_s2.csv", row.names = FALSE)
cat("\n==== Table S2 comparison (payload = León & Anderson 2026a; now = this ladder) ====\n"); print(cmp, row.names = FALSE)

# payload meta timings for reference
cat("\n==== payload meta$timings (0.2.0, Mac Studio M4 Max, 11 workers for FB) ====\n")
print(do.call(rbind, lapply(P$meta$timings, function(x) as.data.frame(x[c("n_family", "boot_sec", "mr_sec", "speedup")]))))

# ---- 3. bootstrap rungs ------------------------------------------------------
cat("\n==== bootstrap rungs (nb_boots = 5) ====\n")
for (r in bids) { x <- rd(r); if (is.null(x)) { cat(r, ": no record\n"); next }
  print(x$row); cat("  SG_CIs:\n"); str(x$SG_CIs) }

# ---- 4. budget arithmetic -------------------------------------------------------
cat("\n==== budget arithmetic (single core, plan = sequential) ====\n")
bud <- do.call(rbind, lapply(names(ids), function(m) {
  mr <- rd(ids[[m]])$row; bb <- rd(bids[[m]])$row
  per <- if (is.null(bb)) NA_real_ else bb$per_rep_sec_total_over_5
  data.frame(method = m, M_fit_sec = mr$fit_sec, M_rung_wall_s = L$wall_s[L$rung == ids[[m]]],
             B5_call_sec = if (is.null(bb)) NA else bb$call_sec, per_rep_sec = per,
             rep1_sec = if (is.null(bb)) NA else bb$rep1_sec,
             B_in_30s = floor(30 / per), B_in_60s = floor(60 / per), B_in_120s = floor(120 / per))
}))
print(bud, row.names = FALSE)
write.csv(bud, "budget.csv", row.names = FALSE)
fsr <- rd("M-FS")$row; fij <- rd("M-FS-ij")$row
cat(sprintf("\nFS-ij saving: fit %.2f s (field) - %.2f s (ij) = %.2f s; rung wall %.2f - %.2f = %.2f s\n",
            fsr$fit_sec, fij$fit_sec, fsr$fit_sec - fij$fit_sec,
            L$wall_s[L$rung == "M-FS"], L$wall_s[L$rung == "M-FS-ij"],
            L$wall_s[L$rung == "M-FS"] - L$wall_s[L$rung == "M-FS-ij"]))
cat(sprintf("Three M fits (as the document ran them): %.2f s fit time total (%s)\n",
            sum(bud$M_fit_sec), paste(sprintf("%s %.2f", bud$method, bud$M_fit_sec), collapse = ", ")))

# ---- 5. overhead-adjusted budget ---------------------------------------------
# call_sec = fixed per-call overhead + sum over replicates of tmins_iteration.
# The overhead (index matrix, H/Hc observed fits, get_dfRes, table build) does
# not scale with B, so B affordable in T s = floor((T - overhead) / mean_iter).
cat("\n==== overhead-adjusted budget (B affordable in T seconds of the bootstrap call alone; M fit excluded) ====\n")
bud2 <- do.call(rbind, lapply(names(bids), function(m) {
  x <- rd(bids[[m]]); it <- x$results$tmins_iteration * 60
  ov <- x$row$call_sec - sum(it); mi <- mean(it)
  data.frame(method = m, call_sec = x$row$call_sec, sum_iter_sec = round(sum(it), 2), overhead_sec = round(ov, 2),
             mean_iter_sec = round(mi, 3), max_iter_sec = round(max(it), 2),
             B_in_30s = floor((30 - ov) / mi), B_in_60s = floor((60 - ov) / mi), B_in_120s = floor((120 - ov) / mi))
}))
print(bud2, row.names = FALSE)
write.csv(bud2, "budget_overhead_adjusted.csv", row.names = FALSE)
