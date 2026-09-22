#!/usr/bin/env Rscript
# nullmr findings (TASK_null_gbsg_mr_2026-09-21, Step 4).
# Reads the 18 nullmr bundles and, READ ONLY, the 18 committed nullid bundles,
# and writes the report's tables to stdout as markdown.  Nothing is re-run.
#
# COVERAGE IS THE COMMITTED ARITHMETIC (Step 1.2, logs/nullmr_step1_record.md):
# summary_grfmr.qmd's `cov-fns` chunk and the `wil` / `covs_cell` definitions
# of its `wilson-fn` chunk are extracted from the committed .qmd and eval'd
# verbatim (the mechanism of grfmr_tables.R:12-21), not re-implemented.
# usage: nullmr_findings.R [campaign-tag]  (default nullmr)
args <- commandArgs(trailingOnly = TRUE)
CAMP <- if (length(args)) args[1] else "nullmr"
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(forestsearch))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

# ---- the committed definitions, verbatim --------------------------------------
src <- readLines("summary_grfmr.qmd", warn = FALSE)
chunk <- function(nm) {
  st <- grep(sprintf("^```\\{r %s[,}]", nm), src); stopifnot(length(st) == 1L)
  en <- st + which(src[(st + 1):length(src)] == "```")[1]
  paste(src[(st + 1):(en - 1)], collapse = "\n") }
labs <- character(0)                       # the chunk's own VSN line then builds nothing
have <- function(x) !is.null(x) && NROW(x) > 0L
eval(parse(text = chunk("cov-fns")), envir = globalenv())      # subst_s, err_sd, cov_block
eval(parse(text = chunk("wilson-fn")), envir = globalenv())    # wil, z95, z975, covs_cell
stopifnot(exists("cov_block"), exists("covs_cell"), exists("wil"))

CELLS <- data.frame(
  cell = c("null0657_n500","null0721_n500","null0657_n1000","null0721_n1000",
           "null0657_n1500","null0721_n1500"),
  hr   = c(0.657, 0.721, 0.657, 0.721, 0.657, 0.721),
  n    = c(500L, 500L, 1000L, 1000L, 1500L, 1500L))
ENG <- c(consistency = "fs", dina = "dina", grf = "grf")
ENGLAB <- c(consistency = "FS", dina = "DINA", grf = "GRF")
REPS <- 2000L
mpath <- function(tag, hr, n)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_%s_res_1_2000.rds",
          tag, round(100 * hr), n, round(1000 * hr), CAMP)
ipath <- function(tag, hr, n)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr_nullid_res_1_2000.rds",
          tag, round(100 * hr), n, round(1000 * hr))

# covs_cell(k) reads B[[k]]$results and BLK / NN / HR / ARM[[k]]
B <- list(); BLK <- list(); NN <- list(); HR <- list(); ARM <- list(); NID <- list(); KEYS <- character(0)
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  k <- paste(CELLS$cell[i], ENGLAB[[e]])
  p <- mpath(ENG[[e]], CELLS$hr[i], CELLS$n[i])
  if (!file.exists(p)) { cat(sprintf("**Bundle absent:** `%s`\n\n", p)); next }
  B[[k]] <- readRDS(p); BLK[[k]] <- CELLS$cell[i]; NN[[k]] <- CELLS$n[i]; HR[[k]] <- CELLS$hr[i]
  ARM[[k]] <- "null"; KEYS <- c(KEYS, k)
  q <- ipath(ENG[[e]], CELLS$hr[i], CELLS$n[i]); NID[[k]] <- if (file.exists(q)) readRDS(q) else NULL
}
f4 <- function(v) ifelse(is.finite(v), formatC(v, format = "f", digits = 4), "-")
f3 <- function(v) ifelse(is.finite(v), formatC(v, format = "f", digits = 3), "-")
cw <- function(p, lo, hi, k) if (!is.finite(p)) sprintf("- [k=%d]", k) else sprintf("%s [%s, %s] (k=%d)", f4(p), f4(lo), f4(hi), k)
fw <- function(x, nn) { if (!nn) return("-"); w <- wil(x / nn, nn); sprintf("%d/%d = %.4f [%.4f, %.4f]", x, nn, x / nn, w[1], w[2]) }

# ---------------------------------------------------------------------------
cat("### Table 1 — coverage of β(Ĥ) and β(Ĥᶜ), over declaring replicates\n\n")
cat("Rate [Wilson 95%] (k = declaring replicates the rate rests on). field / field-s / Bonferroni / IJ: `covs_cell`",
    "(summary_grfmr.qmd:332-357) on its row set (declared, both targets and both one-sided field bounds finite).",
    "Naive: `cov_block()` → `fs_sim_bias_coverage()`, one-sided 95% on the exposed side (cov1) and two-sided (cov2),",
    "each on its own finite rows. IJ two-term is two-sided, as the committed reports define it. DINA and GRF rows are",
    "**conditional on the proposed family**.\n\n")
cat("| cell | identifier | declared | field lower β(Ĥ) | field-s upper β(Ĥᶜ) | Bonferroni pair (joint) | IJ two-term β(Ĥ) 2-sided | IJ two-term β(Ĥᶜ) 2-sided | naive β(Ĥ) 1-sided lower | naive β(Ĥᶜ) 1-sided upper | naive β(Ĥ) 2-sided | naive β(Ĥᶜ) 2-sided |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
COVT <- list()
for (k in KEYS) {
  r <- B[[k]]$results; nd <- sum(r$detected == 1L)
  v <- if (nd) covs_cell(k) else NULL
  nH <- if (nd) cov_block(r, "H", "lower") else NULL
  nC <- if (nd) cov_block(r, "Hc", "upper") else NULL
  gnv <- function(t, col) { x <- t[t$construction == "naive", ]; x[[col]] }
  COVT[[k]] <- list(v = v, nH = nH, nC = nC)
  if (!nd || is.null(v)) { cat(sprintf("| %s | %s | %d | - | - | - | - | - | - | - | - | - |\n", BLK[[k]], sub(".* ", "", k), nd)); next }
  ke <- v$n_eval
  # the naive rows' own n: fs_sim_bias_coverage reports max(n1, n2); cov1 rests on the rows with the bound finite
  n1H <- sum(r$detected == 1L & is.finite(r$betaHhat_H) & is.finite(r$nv_H_est) & is.finite(r$nv_H_se) & r$nv_H_est > 0)
  n1C <- sum(r$detected == 1L & is.finite(r$betaHhat_Hc) & is.finite(r$nv_Hc_est) & is.finite(r$nv_Hc_se) & r$nv_Hc_est > 0)
  n2H <- sum(r$detected == 1L & is.finite(r$betaHhat_H) & is.finite(r$nv_H_lo) & is.finite(r$nv_H_hi))
  n2C <- sum(r$detected == 1L & is.finite(r$betaHhat_Hc) & is.finite(r$nv_Hc_lo) & is.finite(r$nv_Hc_hi))
  cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n",
      BLK[[k]], sub(".* ", "", k), nd,
      cw(v$harm_field, v$harm_field_lo, v$harm_field_hi, ke),
      cw(v$comp_field_s, v$comp_field_s_lo, v$comp_field_s_hi, ke),
      cw(v$joint_s_bonf, v$joint_s_bonf_lo, v$joint_s_bonf_hi, ke),
      cw(v$ij2_H, v$ij2_H_lo, v$ij2_H_hi, ke),
      cw(v$ij2_Hc, v$ij2_Hc_lo, v$ij2_Hc_hi, ke),
      cw(gnv(nH, "cov1"), gnv(nH, "cov1_wilson_lo"), gnv(nH, "cov1_wilson_hi"), n1H),
      cw(gnv(nC, "cov1"), gnv(nC, "cov1_wilson_lo"), gnv(nC, "cov1_wilson_hi"), n1C),
      cw(gnv(nH, "cov2"), gnv(nH, "cov2_wilson_lo"), gnv(nH, "cov2_wilson_hi"), n2H),
      cw(gnv(nC, "cov2"), gnv(nC, "cov2_wilson_lo"), gnv(nC, "cov2_wilson_hi"), n2C)))
}

# ---------------------------------------------------------------------------
cat("\n### Table 2 — where the lower bounds on β(Ĥ) sit (HR scale)\n\n")
cat("Locations against 1.00 and 1.25, not tests. Naive and IJ two-term: the one-sided 95% Gaussian bound",
    "`exp(log est − z_0.95 · SE)` that `fs_sim_bias_coverage()` scores (SE `nv_H_se` / `mr_H_se_ij`); field:",
    "`fld_H_lo1s`; Bonferroni: its lower member `fld_joint_s_bonf_loH` (97.5%). *cond.* over declaring replicates",
    "with the bound finite (k); *uncond.* over all 2,000 (no declaration scores 0). The `nullid` row is the committed",
    "MR-off bundle's unadjusted one-sided bound `nv_H_lo1s` (its `nv_H` refit), printed for reference beside naive.\n\n")
cat("| cell | identifier | product | k | median LB | LB ≥ 1.00 cond. | LB ≥ 1.25 cond. | LB ≥ 1.00 uncond. | LB ≥ 1.25 uncond. |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
LOCT <- list()
lb_of <- function(r, p) {
  d <- r[r$detected == 1L, , drop = FALSE]
  switch(p,
    naive = ifelse(is.finite(d$nv_H_est) & is.finite(d$nv_H_se) & d$nv_H_est > 0, exp(log(d$nv_H_est) - z95 * d$nv_H_se), NA_real_),
    "IJ two-term" = ifelse(is.finite(d$mr_H_est) & is.finite(d$mr_H_se_ij) & d$mr_H_est > 0, exp(log(d$mr_H_est) - z95 * d$mr_H_se_ij), NA_real_),
    field = d$fld_H_lo1s,
    "Bonferroni lower" = d$fld_joint_s_bonf_loH,
    "nullid nv_H (MR off, ref.)" = d$nv_H_lo1s) }
for (k in KEYS) {
  prods <- c("naive", "nullid nv_H (MR off, ref.)", "IJ two-term", "field", "Bonferroni lower")
  for (p in prods) {
    rr <- if (grepl("^nullid", p)) NID[[k]]$results else B[[k]]$results
    if (is.null(rr)) next
    lb <- lb_of(rr, p); kf <- sum(is.finite(lb))
    a1 <- sum(lb >= 1.00, na.rm = TRUE); a2 <- sum(lb >= 1.25, na.rm = TRUE)
    LOCT[[paste(k, p)]] <- data.frame(key = k, product = p, k = kf, med = if (kf) stats::median(lb, na.rm = TRUE) else NA,
                                      a1 = a1, a2 = a2, N = nrow(rr))
    cat(sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s |\n", BLK[[k]], sub(".* ", "", k), p, kf,
        if (kf) f3(stats::median(lb, na.rm = TRUE)) else "-",
        if (kf) sprintf("%d/%d = %.4f", a1, kf, a1 / kf) else "-",
        if (kf) sprintf("%d/%d = %.4f", a2, kf, a2 / kf) else "-",
        fw(a1, nrow(rr)), fw(a2, nrow(rr))))
  }
}

# ---------------------------------------------------------------------------
cat("\n### Table 3 — context per cell: the realized trial against the target, the truths in Ĥ and Ĥᶜ, and mr_ok\n\n")
cat("`itt_est` is the whole-trial Cox HR (engine-independent; Gate C). The truths β(Ĥ), β(Ĥᶜ) are",
    "super-population marginal Cox HRs on uncensored potential outcomes of the selected region and its complement;",
    "the trial fits censored data (`nullc125` §5), so coverage is read beside that gap. Medians over declaring replicates.\n\n")
cat("| cell | identifier | target HR | itt_est median [Q1, Q3] | declared | median β(Ĥ) | median β(Ĥᶜ) | median β(Ĥᶜ) − target | mr_ok among declaring |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (k in KEYS) {
  r <- B[[k]]$results; d <- r[r$detected == 1L, , drop = FALSE]; nd <- nrow(d)
  q <- stats::quantile(r$itt_est, c(.25, .5, .75), names = FALSE)
  nmr <- sum(d$mr_ok == 1L)
  cat(sprintf("| %s | %s | %.3f | %s [%s, %s] | %d | %s | %s | %s | %s |\n", BLK[[k]], sub(".* ", "", k), HR[[k]],
      f4(q[2]), f4(q[1]), f4(q[3]), nd,
      if (nd) f4(stats::median(d$betaHhat_H)) else "-", if (nd) f4(stats::median(d$betaHhat_Hc)) else "-",
      if (nd) sprintf("%+.4f", stats::median(d$betaHhat_Hc) - HR[[k]]) else "-",
      fw(nmr, nd)))
}

# ---------------------------------------------------------------------------
cat("\n### Table 4 — identification against `nullid`, recorded, not gated\n\n")
cat("Per replicate, same `sim_id`, against the committed `nullid` bundle (Mac-Studio, R 4.5.2, MR off). Label and",
    "|Ĥ| are compared on replicates where both declared.\n\n")
cat("| cell | identifier | declared nullmr / nullid | declaration differs | label differs (both declared) | \\|Ĥ\\| differs (both declared) | sg_def differs (both declared) |\n")
cat("|---|---|---|---|---|---|---|\n")
IDT <- list()
for (k in KEYS) {
  if (is.null(NID[[k]])) next
  a <- B[[k]]$results; b <- NID[[k]]$results
  a <- a[order(a$sim_id), ]; b <- b[match(a$sim_id, b$sim_id), ]; stopifnot(identical(a$sim_id, b$sim_id))
  both <- a$detected == 1L & b$detected == 1L
  dd <- sum(a$detected != b$detected)
  dl <- sum(both & !(a$label %in% NA & b$label %in% NA) & (is.na(a$label) | is.na(b$label) | a$label != b$label))
  dn <- sum(both & (is.na(a$n_sel) | is.na(b$n_sel) | a$n_sel != b$n_sel))
  ds <- sum(both & (is.na(a$sg_def) | is.na(b$sg_def) | a$sg_def != b$sg_def))
  IDT[[k]] <- c(dd = dd, dl = dl, dn = dn, ds = ds, both = sum(both))
  cat(sprintf("| %s | %s | %d / %d | %d | %d / %d | %d / %d | %d / %d |\n", BLK[[k]], sub(".* ", "", k),
      sum(a$detected == 1L), sum(b$detected == 1L), dd, dl, sum(both), dn, sum(both), ds, sum(both)))
}

# ---------------------------------------------------------------------------
cat("\n### Table 5 — wall per cell, and the machine, R version and build\n\n")
wall <- function(cell, e) {
  f <- sprintf("scripts_dinamr/logs/%s_%s_%s_effMaxSG.log", CAMP, cell, e)
  if (!file.exists(f)) return(NA_real_)
  w <- grep("^WALL_SECONDS=", readLines(f, warn = FALSE), value = TRUE)
  if (!length(w)) NA_real_ else as.numeric(sub("^WALL_SECONDS=([0-9]+).*", "\\1", w[length(w)])) }
cat("| cell | wall s FS / DINA / GRF | cell s |\n|---|---|---|\n")
tot <- 0
for (c in CELLS$cell) {
  w <- vapply(names(ENG), function(e) wall(c, e), 0); tot <- tot + sum(w, na.rm = TRUE)
  cat(sprintf("| %s | %s | %s |\n", c, paste(ifelse(is.finite(w), sprintf("%.0f", w), "-"), collapse = " / "),
              if (all(is.finite(w))) sprintf("%.0f", sum(w)) else "-"))
}
cat(sprintf("| **total (%d renders)** | | **%.0f s = %.3f h** |\n\n", length(KEYS), tot, tot / 3600))
M <- do.call(rbind, lapply(KEYS, function(k) { m <- B[[k]]$meta
  data.frame(host = m$hostname %||% NA, r = m$r_version %||% NA, fs = m$forestsearch_version %||% NA,
             workers = m$n_workers %||% NA) }))
M <- unique(M)
cat("| hostnames | R | forestsearch | workers | bundles |\n|---|---|---|---|---|\n")
for (i in seq_len(nrow(M))) cat(sprintf("| %s | %s | %s | %s | %d |\n", M$host[i], M$r[i], M$fs[i], M$workers[i], length(KEYS)))

saveRDS(list(COVT = COVT, LOCT = do.call(rbind, LOCT), IDT = IDT), file.path("scripts_dinamr", "nullmr_findings.rds"))
