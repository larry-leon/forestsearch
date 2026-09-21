#!/usr/bin/env Rscript
# nullid findings (TASK_null_gbsg_identification_2026-09-21, Step 4).
# Reads the 18 committed nullid bundles and writes the report's tables and
# findings to stdout as markdown.  Reading only; nothing is re-run.
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

CELLS <- data.frame(
  cell = c("null0657_n500","null0721_n500","null0657_n1000","null0721_n1000",
           "null0657_n1500","null0721_n1500"),
  hr   = c(0.657, 0.721, 0.657, 0.721, 0.657, 0.721),
  n    = c(500L, 500L, 1000L, 1000L, 1500L, 1500L))
ENG <- c(consistency = "fs", dina = "dina", grf = "grf")
ENGLAB <- c(consistency = "FS", dina = "DINA", grf = "GRF")

bpath <- function(tag, hr, n)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr_nullid_res_1_2000.rds",
          tag, round(100 * hr), n, round(1000 * hr))

wilson <- function(x, nn, conf = 0.95) {
  if (nn == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2); p <- x / nn
  d <- 1 + z^2 / nn; c1 <- p + z^2 / (2 * nn)
  h <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2))
  c((c1 - h) / d, (c1 + h) / d)
}
mcse  <- function(v) { v <- v[is.finite(v)]; if (length(v) < 2L) NA_real_ else stats::sd(v) / sqrt(length(v)) }
q3    <- function(v) stats::quantile(v[is.finite(v)], c(0.25, 0.5, 0.75), names = FALSE)
f2    <- function(v, d = 3) formatC(v, format = "f", digits = d)
fq    <- function(v, d = 3) { q <- q3(v); sprintf("%s / %s / %s", f2(q[1], d), f2(q[2], d), f2(q[3], d)) }

# wall per cell, from the render logs the driver wrote
cellwall <- function(cell) {
  ws <- vapply(names(ENG), function(e) {
    f <- sprintf("scripts_dinamr/logs/nullid_%s_%s_effMaxSG.log", cell, e)
    if (!file.exists(f)) return(NA_real_)
    l <- grep("WALL_SECONDS", readLines(f, warn = FALSE), value = TRUE)
    if (!length(l)) NA_real_ else as.numeric(sub(".*WALL_SECONDS=([0-9]+).*", "\\1", l[1]))
  }, numeric(1))
  sum(ws, na.rm = TRUE)
}

B <- list()
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  p <- bpath(ENG[[e]], CELLS$hr[i], CELLS$n[i])
  B[[paste(CELLS$cell[i], e)]] <- if (file.exists(p)) readRDS(p) else NULL
}

cat("### Table 1 — declaration, size and specificity (cells x identifiers)\n\n")
cat("| cell | HR | n | identifier | declarations / 2000 | rate [Wilson 95%] |",
    "mean \\|H\\| (MC SE) | \\|H\\| Q1/med/Q3 | mean \\|H\\|/n (MC SE) |",
    "spec uncond. (MC SE) | spec cond. (MC SE) |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  b <- B[[paste(CELLS$cell[i], e)]]; if (is.null(b)) next
  r <- b$results; n <- CELLS$n[i]; N <- nrow(r)
  d <- r[r$detected == 1L, , drop = FALSE]; nd <- nrow(d)
  w <- wilson(nd, N)
  su <- ifelse(r$detected == 1L, r$spec, 1)          # unconditional: no declaration scores 1
  sc <- d$spec                                        # conditional: declaring replicates only
  cat(sprintf("| %s | %.3f | %d | %s | %d | %.4f [%.4f, %.4f] | %s (%s) | %s | %s (%s) | %s (%s) | %s (%s) |\n",
      CELLS$cell[i], CELLS$hr[i], n, ENGLAB[[e]], nd, nd / N, w[1], w[2],
      if (nd) f2(mean(d$n_sel), 1) else "-", if (nd) f2(mcse(d$n_sel), 2) else "-",
      if (nd) sprintf("%.0f / %.0f / %.0f", q3(d$n_sel)[1], q3(d$n_sel)[2], q3(d$n_sel)[3]) else "-",
      if (nd) f2(mean(d$n_sel) / n, 4) else "-", if (nd) f2(mcse(d$n_sel / n), 4) else "-",
      f2(mean(su, na.rm = TRUE), 4), f2(mcse(su), 4),
      if (nd) f2(mean(sc, na.rm = TRUE), 4) else "-", if (nd) f2(mcse(sc), 4) else "-"))
}

cat("\n### Table 2 — the unadjusted within-region estimate, and where its one-sided bound lands\n\n")
cat("| cell | identifier | HR(H) Q1/med/Q3 | model-based SE med | lower 1s Q1/med/Q3 |",
    "share lower >= 1.00 [Wilson] | share lower >= 1.25 [Wilson] |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  b <- B[[paste(CELLS$cell[i], e)]]; if (is.null(b)) next
  d <- b$results[b$results$detected == 1L, , drop = FALSE]; nd <- nrow(d)
  if (!nd) { cat(sprintf("| %s | %s | - | - | - | - | - |\n", CELLS$cell[i], ENGLAB[[e]])); next }
  lo <- d$nv_H_lo1s; k <- sum(is.finite(lo)); a <- sum(lo >= 1.00, na.rm = TRUE); bb <- sum(lo >= 1.25, na.rm = TRUE)
  wa <- wilson(a, k); wb <- wilson(bb, k)
  cat(sprintf("| %s | %s | %s | %s | %s | %d/%d = %.4f [%.4f, %.4f] | %d/%d = %.4f [%.4f, %.4f] |\n",
      CELLS$cell[i], ENGLAB[[e]], fq(d$nv_H_est), f2(stats::median(d$nv_H_se, na.rm = TRUE)),
      fq(lo), a, k, a / k, wa[1], wa[2], bb, k, bb / k, wb[1], wb[2]))
}

cat("\n### Table 3 — the candidate family and the screen statistic (FS only)\n\n")
cat("| cell | enumerated Q1/med/Q3 | clearing the floor Q1/med/Q3 | consistency-qualifying Q1/med/Q3 |",
    "p_sel med | p_max_qual med | floor>0 but no declaration |\n")
cat("|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) {
  b <- B[[paste(CELLS$cell[i], "consistency")]]; if (is.null(b)) next
  r <- b$results; d <- r[r$detected == 1L, , drop = FALSE]
  declined <- sum(r$n_cand_floor > 0L & r$detected == 0L, na.rm = TRUE)
  dn <- sum(r$n_cand_floor > 0L, na.rm = TRUE); wd <- wilson(declined, dn)
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %d/%d = %.4f [%.4f, %.4f] |\n",
      CELLS$cell[i],
      sprintf("%.0f / %.0f / %.0f", q3(r$n_cand_enum)[1], q3(r$n_cand_enum)[2], q3(r$n_cand_enum)[3]),
      sprintf("%.0f / %.0f / %.0f", q3(r$n_cand_floor)[1], q3(r$n_cand_floor)[2], q3(r$n_cand_floor)[3]),
      if (nrow(d)) sprintf("%.0f / %.0f / %.0f", q3(d$n_cons_qual)[1], q3(d$n_cons_qual)[2], q3(d$n_cons_qual)[3]) else "-",
      if (nrow(d)) f2(stats::median(d$p_sel, na.rm = TRUE)) else "-",
      if (nrow(d)) f2(stats::median(d$p_max_qual, na.rm = TRUE)) else "-",
      declined, dn, declined / dn, wd[1], wd[2]))
}

cat("\n### Table 4 — max_g T_g over the screened family, against z_0.95 = 1.645 (FS only)\n\n")
cat("| cell | n with max_g T_g | Q1 | median | Q3 | 90% | 95% | 99% | share > 1.645 [Wilson] |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) {
  b <- B[[paste(CELLS$cell[i], "consistency")]]; if (is.null(b)) next
  v <- b$results$maxT; v <- v[is.finite(v)]; k <- length(v)
  qq <- stats::quantile(v, c(0.25, 0.5, 0.75, 0.90, 0.95, 0.99), names = FALSE)
  a <- sum(v > 1.6448536); w <- wilson(a, k)
  cat(sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %d/%d = %.4f [%.4f, %.4f] |\n",
      CELLS$cell[i], k, f2(qq[1]), f2(qq[2]), f2(qq[3]), f2(qq[4]), f2(qq[5]), f2(qq[6]),
      a, k, a / k, w[1], w[2]))
}

cat("\n### Table 5 — the composition of H-hat: covariates and cut directions\n\n")
parse_terms <- function(sg) {
  if (!length(sg) || is.na(sg)) return(character(0))
  tt <- regmatches(sg, gregexpr("!?\\{[^}]*\\}", sg))[[1]]
  vapply(tt, function(s) {
    neg <- startsWith(s, "!")
    body <- gsub("^!?\\{|\\}$", "", s)
    v  <- sub("^\\s*([A-Za-z._][A-Za-z0-9._]*).*$", "\\1", body)
    op <- if (grepl("<=", body)) "<=" else if (grepl(">=", body)) ">=" else
          if (grepl("==", body)) "==" else if (grepl("<", body)) "<" else
          if (grepl(">", body)) ">" else "?"
    sprintf("%s %s%s", v, if (neg) "NOT " else "", op)
  }, character(1), USE.NAMES = FALSE)
}
for (e in names(ENG)) {
  cat(sprintf("\n**%s** — share of declaring replicates whose rule contains the term, pooled over the six cells:\n\n", ENGLAB[[e]]))
  all_terms <- list(); tot <- 0L
  for (i in seq_len(nrow(CELLS))) {
    b <- B[[paste(CELLS$cell[i], e)]]; if (is.null(b)) next
    d <- b$results[b$results$detected == 1L, , drop = FALSE]
    tot <- tot + nrow(d)
    for (s in d$sg_def) all_terms[[length(all_terms) + 1L]] <- unique(parse_terms(s))
  }
  tb <- sort(table(unlist(all_terms)), decreasing = TRUE)
  cat("| term | replicates | share of declaring |\n|---|---|---|\n")
  for (k in names(tb)) cat(sprintf("| `%s` | %d | %.4f |\n", k, tb[[k]], tb[[k]] / tot))
  cat(sprintf("\n(declaring replicates pooled: %d)\n", tot))
}

cat("\n### Table 6 — measured wall per cell\n\n")
cat("| cell | wall (s) | wall (h) |\n|---|---|---|\n")
tw <- 0
for (i in seq_len(nrow(CELLS))) {
  w <- cellwall(CELLS$cell[i]); tw <- tw + w
  cat(sprintf("| %s | %.0f | %.3f |\n", CELLS$cell[i], w, w / 3600))
}
cat(sprintf("| **total (18 renders)** | **%.0f** | **%.3f** |\n", tw, tw / 3600))

b1 <- B[[paste(CELLS$cell[1], "consistency")]]
cat(sprintf("\nBuild: forestsearch %s ; host %s ; workers %s ; R %s.%s\n",
            b1$meta$forestsearch_version, b1$meta$host %||% Sys.info()[["nodename"]],
            b1$meta$n_workers, R.version$major, R.version$minor))
