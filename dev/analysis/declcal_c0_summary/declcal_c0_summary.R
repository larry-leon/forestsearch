# declcal_c0_summary.R -- S1.8 payload: every declaration-calibration table value,
# computed from the committed declcal / declcalc0 payloads
# (TASK_declcal_c0_summary_2026-09-25).
#
# Read-only on the payloads.  Types no number: every value is a mean, count,
# maximum or median of a payload column.  Column semantics: scripts_dinamr/
# declcalc0_run.R@a46bf9b7 and REPORT_declcal_c0_campaign_reads_2026-09-24.md, Q6.
#
# usage (from the repository root): Rscript dev/analysis/declcal_c0_summary/declcal_c0_summary.R
# writes quarto/simulations/gbsg_020/results/declcalc0_summary_tables.csv

res_dir <- "quarto/simulations/gbsg_020/results"
out_csv <- file.path(res_dir, "declcalc0_summary_tables.csv")
if (!dir.exists(res_dir)) stop("run from the repository root", call. = FALSE)

c0s    <- c("c070", "c075", "c080", "c085")
alphas <- c("05" = 0.05, "10" = 0.10)
k_fix  <- 2.0   # fixed cutoff on the executed family (the record's "tuned fixed p* = 0.9545")

B_ids <- paste0("B", 1:6)
C_ids <- paste0("C", 1:4)
load_set <- function(prefix) {
  files <- c(sprintf("%s_inull_%s_res_1_2000.rds", prefix, B_ids),
             sprintf("%s_power_%s_res_1_2000.rds", prefix, C_ids))
  lapply(setNames(files, c(B_ids, C_ids)), function(f) {
    x <- readRDS(file.path(res_dir, f))
    stopifnot(identical(x$meta$cell_status, "complete"), nrow(x$results) == 2000L,
              all(x$results$status == "ok"))
    list(file = f, r = x$results, design_hr = x$meta$target_hr, n = x$meta$n,
         cell = x$meta$cell_id, p_star = x$meta$p_star)
  })
}
c0p <- load_set("declcalc0")   # protected-level (shifted) rule
clp <- load_set("declcal")     # claim-threshold (unshifted) rule

rows <- list()
add <- function(table, row, col, alpha, c0, cell, design_hr, n, value, count,
                denominator, source_column, source_payload) {
  rows[[length(rows) + 1L]] <<- data.frame(
    table = table, row = row, col = col, alpha = alpha, c0 = c0, cell = cell,
    design_hr = design_hr, n = n, value = value, count = count,
    denominator = denominator, source_column = source_column,
    source_payload = source_payload, stringsAsFactors = FALSE)
}
# indicator of a column; NA counts as not declared
ind <- function(p, col) {
  v <- p$r[[col]]
  if (is.null(v)) stop("column ", col, " not in ", p$file, call. = FALSE)
  as.integer(!is.na(v) & v == 1L)
}
ind_fix <- function(p) as.integer(!is.na(p$r$max_T_post) & p$r$max_T_post >= k_fix)
rate_row <- function(table, row, alpha, c0, p, x, src_col) {
  add(table, row, p$cell, alpha, c0, p$cell, p$design_hr, p$n,
      mean(x), sum(x), length(x), src_col, p$file)
}
# the B maximum, ties reported as "B2;B5"
max_row <- function(table, row, alpha, c0, set, xfun, src_col) {
  cnt <- vapply(B_ids, function(b) sum(xfun(set[[b]])), numeric(1))
  den <- vapply(B_ids, function(b) nrow(set[[b]]$r), numeric(1))
  rt  <- cnt / den
  at  <- B_ids[rt == max(rt)]
  add(table, row, "max_B", alpha, c0, paste(at, collapse = ";"),
      paste(unique(vapply(set[at], `[[`, 0, "design_hr")), collapse = ";"),
      paste(unique(vapply(set[at], `[[`, 0L, "n")), collapse = ";"),
      max(rt), cnt[at[1]], den[at[1]], src_col,
      paste(vapply(set[at], `[[`, "", "file"), collapse = ";"))
}

# ---- (a) Table S3: conventional screen, B1-B6 -------------------------------
for (b in B_ids) rate_row("S3", "conventional", NA, NA, c0p[[b]],
                          ind(c0p[[b]], "declared_conv"), "declared_conv")

# ---- (b) Tables S4 (alpha 0.10) and S5 (alpha 0.05) --------------------------
for (a in names(alphas)) {
  tb <- if (a == "10") "S4" else "S5"
  al <- alphas[[a]]
  for (c0 in c0s) {
    col <- sprintf("declared_cal%s_%s", a, c0)
    row <- sprintf("calibrated_%s", c0)
    for (b in B_ids) rate_row(tb, row, al, c0, c0p[[b]], ind(c0p[[b]], col), col)
    max_row(tb, row, al, c0, c0p, function(p) ind(p, col), col)
    for (cc in C_ids) rate_row(tb, row, al, c0, c0p[[cc]], ind(c0p[[cc]], col), col)
  }
  # comparators
  max_row(tb, "conventional", al, NA, c0p, function(p) ind(p, "declared_conv"), "declared_conv")
  for (cc in C_ids) rate_row(tb, "conventional", al, NA, c0p[[cc]],
                             ind(c0p[[cc]], "declared_conv"), "declared_conv")
  fcol <- sprintf("max_T_post>=%.1f", k_fix)
  for (b in B_ids) rate_row(tb, "fixed_k2.0", al, NA, c0p[[b]], ind_fix(c0p[[b]]), fcol)
  max_row(tb, "fixed_k2.0", al, NA, c0p, ind_fix, fcol)
  for (cc in C_ids) rate_row(tb, "fixed_k2.0", al, NA, c0p[[cc]], ind_fix(c0p[[cc]]), fcol)
  ccol <- sprintf("declared_cal%s", a)
  for (b in B_ids) rate_row(tb, "claim_threshold", al, NA, clp[[b]], ind(clp[[b]], ccol), ccol)
  max_row(tb, "claim_threshold", al, NA, clp, function(p) ind(p, ccol), ccol)
  for (cc in C_ids) rate_row(tb, "claim_threshold", al, NA, clp[[cc]], ind(clp[[cc]], ccol), ccol)

  # implied p* of the comparators: conventional from meta$p_star (must agree across
  # cells); fixed cutoff 2*pnorm(k_fix)-1; claim threshold pooled over every cell and
  # every n, under both kappa summaries (the per-n values are in S6)
  all_ids <- c(B_ids, C_ids)
  allf <- paste(vapply(c0p[all_ids], `[[`, "", "file"), collapse = ";")
  ps <- unique(vapply(c0p[all_ids], `[[`, 0, "p_star"))
  if (length(ps) != 1L) stop("meta$p_star differs across cells", call. = FALSE)
  add(tb, "conventional", "implied_pstar", al, NA, paste(all_ids, collapse = ";"), NA, NA,
      ps, NA, length(all_ids), "meta$p_star", allf)
  add(tb, "fixed_k2.0", "implied_pstar", al, NA, NA, NA, NA,
      2 * pnorm(k_fix) - 1, NA, NA, sprintf("2*pnorm(%.1f)-1", k_fix), NA)
  kcol <- sprintf("kappa_hat_%s", a)
  kl <- lapply(clp[all_ids], function(p) p$r[[kcol]])
  kp <- unlist(kl, use.names = FALSE)
  cands <- list(pooled_median = median(kp, na.rm = TRUE),
                median_of_cell_medians = median(vapply(kl, median, 0, na.rm = TRUE)))
  for (def in names(cands))
    add(tb, "claim_threshold", paste0("implied_pstar_", def, "_allcells_alln"), al, NA,
        paste(all_ids, collapse = ";"), NA, NA, 2 * pnorm(cands[[def]]) - 1, NA,
        if (def == "pooled_median") sum(!is.na(kp)) else length(all_ids),
        sprintf("2*pnorm(%s)-1", kcol),
        paste(vapply(clp[all_ids], `[[`, "", "file"), collapse = ";"))
}

# ---- (c) Table S6: kappa_hat by sample size, four candidate definitions ------
# two summaries (pooled per-replicate median; median of per-cell medians) x two
# cell sets at each n (every cell, B and C; B cells only)
all_ids  <- c(B_ids, C_ids)
cellsets <- list(allcells = all_ids, Bcells = B_ids)
ns <- sort(unique(vapply(c0p, `[[`, 0L, "n")))
kappa_rows <- function(set, kcol, row, al, c0) {
  for (nn in ns) for (cs in names(cellsets)) {
    ids <- intersect(cellsets[[cs]], all_ids[vapply(set, `[[`, 0L, "n") == nn])
    k   <- lapply(set[ids], function(p) {
      v <- p$r[[kcol]]
      if (is.null(v)) stop("column ", kcol, " not in ", p$file, call. = FALSE)
      v })
    pooled <- unlist(k, use.names = FALSE)
    permed <- vapply(k, median, 0, na.rm = TRUE)
    cands  <- list(pooled_median = median(pooled, na.rm = TRUE),
                   median_of_cell_medians = median(permed))
    dens   <- c(pooled_median = sum(!is.na(pooled)), median_of_cell_medians = length(ids))
    for (def in names(cands)) {
      kk <- cands[[def]]
      common <- list(alpha = al, c0 = c0, cell = paste(ids, collapse = ";"),
                     design_hr = paste(unique(vapply(set[ids], `[[`, 0, "design_hr")), collapse = ";"),
                     n = nn, source_payload = paste(vapply(set[ids], `[[`, "", "file"), collapse = ";"))
      add("S6", row, paste0("kappa_", def, "_", cs), al, c0, common$cell, common$design_hr, nn,
          kk, NA, dens[[def]], kcol, common$source_payload)
      add("S6", row, paste0("level_", def, "_", cs), al, c0, common$cell, common$design_hr, nn,
          2 * pnorm(kk) - 1, NA, dens[[def]], sprintf("2*pnorm(%s)-1", kcol), common$source_payload)
    }
  }
}
for (a in names(alphas)) {
  al <- alphas[[a]]
  for (c0 in c0s) kappa_rows(c0p, sprintf("kappa_hat_%s_%s", a, c0),
                             sprintf("calibrated_%s", c0), al, c0)
  kappa_rows(clp, sprintf("kappa_hat_%s", a), "claim_threshold", al, NA)
}

# ---- (d) re-selection footprint: C cells, alpha 0.05 -------------------------
pool_num <- 0; pool_den <- 0
for (c0 in c0s) {
  dcol <- sprintf("declared_cal05_%s", c0); acol <- sprintf("n_admitted_cal05_%s", c0)
  num <- 0; den <- 0
  for (cc in C_ids) {
    p <- c0p[[cc]]; d <- ind(p, dcol) == 1L
    num <- num + sum(!is.na(p$r[[acol]][d]) & p$r[[acol]][d] > 1L); den <- den + sum(d)
  }
  add("footprint", sprintf("calibrated_%s", c0), "share_n_admitted_gt1", 0.05, c0,
      paste(C_ids, collapse = ";"), NA, NA, num / den, num, den,
      sprintf("%s>1 | %s==1", acol, dcol),
      paste(vapply(c0p[C_ids], `[[`, "", "file"), collapse = ";"))
  pool_num <- pool_num + num; pool_den <- pool_den + den
}
add("footprint", "pooled_c0", "share_n_admitted_gt1", 0.05, "c070;c075;c080;c085",
    paste(C_ids, collapse = ";"), NA, NA, pool_num / pool_den, pool_num, pool_den,
    "n_admitted_cal05_<c0>>1 | declared_cal05_<c0>==1",
    paste(vapply(c0p[C_ids], `[[`, "", "file"), collapse = ";"))

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE, na = "")
cat("wrote", out_csv, ":", nrow(out), "rows\n")
