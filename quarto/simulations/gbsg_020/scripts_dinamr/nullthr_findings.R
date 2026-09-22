#!/usr/bin/env Rscript
# nullc125 findings (TASK_null_gbsg_thresholds_2026-09-21_v2, Step 4).
# Reads the 18 nullc125 bundles and, READ ONLY, the 18 committed nullid
# bundles (the 0.90 / 0.80 screen), and writes the report's tables to stdout
# as markdown.  Nothing is re-run.
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
SCR <- data.frame(key = c("nullid", "nullc125"), lab = c("0.90/0.80", "1.25/1.00"),
                  thr = c("", "_c125c100"))
REPS <- 2000L

bpath <- function(tag, hr, n, s)
  sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr%s_%s_res_1_2000.rds",
          tag, round(100 * hr), n, round(1000 * hr), SCR$thr[SCR$key == s], s)

wilson <- function(x, nn, conf = 0.95) {
  if (nn == 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf) / 2); p <- x / nn
  d <- 1 + z^2 / nn; c1 <- p + z^2 / (2 * nn)
  h <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2))
  c((c1 - h) / d, (c1 + h) / d)
}
mcse <- function(v) { v <- v[is.finite(v)]; if (length(v) < 2L) NA_real_ else stats::sd(v) / sqrt(length(v)) }
q3   <- function(v) { v <- v[is.finite(v)]; if (!length(v)) rep(NA_real_, 3) else stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE) }
f2   <- function(v, d = 3) ifelse(is.finite(v), formatC(v, format = "f", digits = d), "-")
fq0  <- function(v) { q <- q3(v); if (all(is.na(q))) "-" else sprintf("%.0f / %.0f / %.0f", q[1], q[2], q[3]) }
fq   <- function(v, d = 3) { q <- q3(v); if (all(is.na(q))) "-" else sprintf("%s / %s / %s", f2(q[1], d), f2(q[2], d), f2(q[3], d)) }
fw   <- function(x, nn) { if (!nn) return("-"); w <- wilson(x, nn); sprintf("%d/%d = %.4f [%.4f, %.4f]", x, nn, x / nn, w[1], w[2]) }

B <- list()
for (s in SCR$key) for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  p <- bpath(ENG[[e]], CELLS$hr[i], CELLS$n[i], s)
  B[[paste(s, CELLS$cell[i], e)]] <- if (file.exists(p)) readRDS(p) else NULL
}
missing <- names(B)[vapply(B, is.null, logical(1))]
all_keys <- as.vector(outer(outer(SCR$key, CELLS$cell, paste), names(ENG), paste))
missing <- setdiff(all_keys, names(B)[!vapply(B, is.null, logical(1))])
if (length(missing)) cat(sprintf("**Bundles absent:** %s\n\n", paste(missing, collapse = "; ")))
get <- function(s, c, e) B[[paste(s, c, e)]]

# ---------------------------------------------------------------------------
cat("### Table 1 — screen x cell x identifier\n\n")
cat("Every declaration is false here. `k` is the number of declaring replicates a conditional summary rests on.",
    "Lower-bound shares: *cond.* over declaring replicates with a finite bound, *uncond.* over all 2,000 (a replicate",
    "that declares nothing scores 0).\n\n")
cat("| cell | identifier | screen c1/c2 | declared / 2000 | rate [Wilson 95%] | mean \\|H\\|/n (MC SE) | spec uncond. (MC SE) |",
    "spec cond. (MC SE) | med HR(H) unadj. | med true β(Ĥ) (HR scale) [populated] | LB ≥ 1.00 cond. | LB ≥ 1.25 cond. |",
    "LB ≥ 1.00 uncond. | LB ≥ 1.25 uncond. |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) for (s in SCR$key) {
  b <- get(s, CELLS$cell[i], e); if (is.null(b)) next
  r <- b$results; n <- CELLS$n[i]; N <- nrow(r)
  d <- r[r$detected == 1L, , drop = FALSE]; nd <- nrow(d)
  w <- wilson(nd, N)
  su <- ifelse(r$detected == 1L, r$spec, 1)
  lo <- d$nv_H_lo1s; kf <- sum(is.finite(lo))
  a1 <- sum(lo >= 1.00, na.rm = TRUE); a2 <- sum(lo >= 1.25, na.rm = TRUE)
  bh <- d$betaHhat_H; kb <- sum(is.finite(bh))
  cat(sprintf("| %s | %s | %s | %d | %.4f [%.4f, %.4f] | %s | %s (%s) | %s | %s | %s | %s | %s | %s | %s |\n",
      CELLS$cell[i], ENGLAB[[e]], SCR$lab[SCR$key == s], nd, nd / N, w[1], w[2],
      if (nd) sprintf("%s (%s) [k=%d]", f2(mean(d$n_sel) / n, 4), f2(mcse(d$n_sel / n), 4), nd) else "- [k=0]",
      f2(mean(su, na.rm = TRUE), 4), f2(mcse(su), 4),
      if (nd) sprintf("%s (%s) [k=%d]", f2(mean(d$spec, na.rm = TRUE), 4), f2(mcse(d$spec), 4), nd) else "- [k=0]",
      if (nd) sprintf("%s [k=%d]", f2(stats::median(d$nv_H_est, na.rm = TRUE)), sum(is.finite(d$nv_H_est))) else "-",
      if (kb) sprintf("%s [%d/%d]", f2(stats::median(bh, na.rm = TRUE)), kb, nd) else sprintf("- [0/%d]", nd),
      if (kf) sprintf("%d/%d = %.4f", a1, kf, a1 / kf) else "-",
      if (kf) sprintf("%d/%d = %.4f", a2, kf, a2 / kf) else "-",
      fw(a1, N), fw(a2, N)))
}

# ---------------------------------------------------------------------------
cat("\n### Table 2 — FS: the candidate family, the consistency screen, and max_g T_g against 1.645\n\n")
cat("`clearing the floor` is at the screen's own c1, so the family `max_g T_g` is taken over differs by screen.\n\n")
cat("| cell | screen | enumerated Q1/med/Q3 | clearing the floor Q1/med/Q3 | empty floor family | consistency-qualifying Q1/med/Q3 [k] |",
    "floor > 0 but no declaration [Wilson] | max_g T_g Q1/med/Q3 | 90% / 95% / 99% | share > 1.645 [Wilson] |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) for (s in SCR$key) {
  b <- get(s, CELLS$cell[i], "consistency"); if (is.null(b)) next
  r <- b$results; d <- r[r$detected == 1L, , drop = FALSE]
  declined <- sum(r$n_cand_floor > 0L & r$detected == 0L, na.rm = TRUE)
  dn <- sum(r$n_cand_floor > 0L, na.rm = TRUE)
  v <- r$maxT[is.finite(r$maxT)]; k <- length(v)
  qq <- if (k) stats::quantile(v, c(0.90, 0.95, 0.99), names = FALSE) else rep(NA_real_, 3)
  cat(sprintf("| %s | %s | %s | %s | %d/2000 | %s | %s | %s | %s / %s / %s | %s |\n",
      CELLS$cell[i], SCR$lab[SCR$key == s], fq0(r$n_cand_enum), fq0(r$n_cand_floor),
      sum(r$n_cand_floor == 0L, na.rm = TRUE),
      if (nrow(d)) sprintf("%s [k=%d]", fq0(d$n_cons_qual), nrow(d)) else "- [k=0]",
      fw(declined, dn), fq(v), f2(qq[1]), f2(qq[2]), f2(qq[3]), fw(sum(v > 1.6448536), k)))
}

# ---------------------------------------------------------------------------
cat("\n### Table 3 — per cell: the realized whole-trial ITT HR (`itt_est`), and the cross-machine draw check\n\n")
cat("`itt_est` exists on the nullc125 bundles only (Gate C: identical across the three identifiers on all 2,000 rows).",
    "nullid has no `itt_est`, but its `or_Hc_est` is the same Cox fit of treatment alone on the whole trial, written on",
    "its declaring rows; the last column compares the two on the rows where nullid's FS declared.\n\n")
cat("| cell | target marginal HR | itt_est Q1 / median / Q3 | mean (MC SE) | nullid or_Hc_est vs nullc125 itt_est (FS rows) |\n")
cat("|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) {
  b <- get("nullc125", CELLS$cell[i], "consistency"); if (is.null(b)) next
  r <- b$results[order(b$results$sim_id), ]
  x <- get("nullid", CELLS$cell[i], "consistency")
  cmp <- "-"
  if (!is.null(x)) {
    rx <- x$results[order(x$results$sim_id), ]
    ok <- is.finite(rx$or_Hc_est)
    dd <- abs(unname(rx$or_Hc_est[ok]) - r$itt_est[ok])
    cmp <- sprintf("%d rows; max |diff| %.3e; identical %s", sum(ok), max(dd), identical(unname(rx$or_Hc_est[ok]), r$itt_est[ok]))
  }
  cat(sprintf("| %s | %.3f | %s | %s (%s) | %s |\n", CELLS$cell[i], CELLS$hr[i], fq(r$itt_est, 4),
              f2(mean(r$itt_est), 4), f2(mcse(r$itt_est), 4), cmp))
}

# ---------------------------------------------------------------------------
cat("\n### Table 4 — the composition of Ĥ, per identifier and screen\n\n")
parse_terms <- function(sg) {
  if (!length(sg) || is.na(sg)) return(character(0))
  tt <- regmatches(sg, gregexpr("!?\\{[^}]*\\}", sg))[[1]]
  vapply(tt, function(s) {
    neg <- startsWith(s, "!")
    body <- gsub("^!?\\{|\\}$", "", s)
    v  <- sub("^\\s*([A-Za-z._][A-Za-z0-9._]*).*$", "\\1", body)
    op <- if (grepl("<=", body)) "<=" else if (grepl(">=", body)) ">=" else
          if (grepl("==", body)) "==" else if (grepl("<", body)) "<" else
          if (grepl(">", body)) ">" else ""
    if (!nzchar(op)) sprintf("%s%s (indicator)", if (neg) "NOT " else "", v)
    else sprintf("%s%s %s", if (neg) "NOT " else "", v, op)
  }, character(1), USE.NAMES = FALSE)
}
comp <- function(s, e) {
  terms <- list(); tot <- 0L; nf <- integer(0)
  for (i in seq_len(nrow(CELLS))) {
    b <- get(s, CELLS$cell[i], e); if (is.null(b)) next
    d <- b$results[b$results$detected == 1L, , drop = FALSE]
    tot <- tot + nrow(d)
    for (sg in d$sg_def) { u <- unique(parse_terms(sg)); terms[[length(terms) + 1L]] <- u; nf <- c(nf, length(u)) }
  }
  list(tb = table(unlist(terms)), tot = tot, nf = nf)
}
for (e in names(ENG)) {
  cs <- lapply(SCR$key, comp, e = e); names(cs) <- SCR$key
  allt <- union(names(cs$nullid$tb), names(cs$nullc125$tb))
  sh <- function(s, t) { x <- cs[[s]]$tb; if (t %in% names(x) && cs[[s]]$tot) x[[t]] / cs[[s]]$tot else 0 }
  ord <- allt[order(-vapply(allt, sh, numeric(1), s = "nullc125"), -vapply(allt, sh, numeric(1), s = "nullid"))]
  cat(sprintf("\n**%s** — share of declaring replicates whose rule contains the term, pooled over the six cells",
              ENGLAB[[e]]))
  cat(sprintf(" (declaring replicates pooled: 0.90/0.80 %d, 1.25/1.00 %d; one-factor rules: %s / %s):\n\n",
              cs$nullid$tot, cs$nullc125$tot,
              fw(sum(cs$nullid$nf == 1L), length(cs$nullid$nf)), fw(sum(cs$nullc125$nf == 1L), length(cs$nullc125$nf))))
  cat("| term | 0.90/0.80 replicates (share) | 1.25/1.00 replicates (share) |\n|---|---|---|\n")
  for (t in ord) {
    g <- function(s) { x <- cs[[s]]$tb; k <- if (t %in% names(x)) x[[t]] else 0L
                       sprintf("%d (%.4f)", k, if (cs[[s]]$tot) k / cs[[s]]$tot else NA) }
    cat(sprintf("| `%s` | %s | %s |\n", t, g("nullid"), g("nullc125")))
  }
}

# ---------------------------------------------------------------------------
cat("\n### Table 5 — wall per cell, and the machine, R version and build of each campaign\n\n")
cellwall <- function(s, cell) {
  vapply(names(ENG), function(e) {
    f <- sprintf("scripts_dinamr/logs/%s_%s_%s_effMaxSG.log", s, cell, e)
    if (!file.exists(f)) return(NA_real_)
    l <- grep("WALL_SECONDS", readLines(f, warn = FALSE), value = TRUE)
    if (!length(l)) NA_real_ else as.numeric(sub(".*WALL_SECONDS=([0-9]+).*", "\\1", l[length(l)]))
  }, numeric(1))
}
cat("| cell | 0.90/0.80 wall s (FS / DINA / GRF) | 0.90/0.80 cell s | 1.25/1.00 wall s (FS / DINA / GRF) | 1.25/1.00 cell s |\n")
cat("|---|---|---|---|---|\n")
tot <- c(nullid = 0, nullc125 = 0)
for (i in seq_len(nrow(CELLS))) {
  w <- lapply(SCR$key, cellwall, cell = CELLS$cell[i]); names(w) <- SCR$key
  for (s in SCR$key) tot[[s]] <- tot[[s]] + sum(w[[s]], na.rm = TRUE)
  cat(sprintf("| %s | %s | %.0f | %s | %.0f |\n", CELLS$cell[i],
              paste(f2(w$nullid, 0), collapse = " / "), sum(w$nullid, na.rm = TRUE),
              paste(f2(w$nullc125, 0), collapse = " / "), sum(w$nullc125, na.rm = TRUE)))
}
cat(sprintf("| **total (18 renders)** | | **%.0f s = %.3f h** | | **%.0f s = %.3f h** |\n",
            tot[["nullid"]], tot[["nullid"]] / 3600, tot[["nullc125"]], tot[["nullc125"]] / 3600))
cat("\n| campaign | screen | hostnames | R | forestsearch | workers | bundles |\n|---|---|---|---|---|---|---|\n")
for (s in SCR$key) {
  ms <- lapply(B[grepl(paste0("^", s, " "), names(B))], `[[`, "meta"); ms <- ms[!vapply(ms, is.null, logical(1))]
  u <- function(k) paste(unique(vapply(ms, function(m) as.character(m[[k]] %||% NA), character(1))), collapse = ", ")
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %d |\n", s, SCR$lab[SCR$key == s], u("hostname"),
              u("r_version"), u("forestsearch_version"), u("n_workers"), length(ms)))
}
