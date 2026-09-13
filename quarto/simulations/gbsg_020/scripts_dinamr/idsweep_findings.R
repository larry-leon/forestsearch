# idsweep findings readout (TASK_idsweep_2026-09-12): the compact tables
# REPORT_idsweep quotes.  Reads the idsweep bundles and idsweep_gateI.rds only;
# no earlier bundle is used as data.  summary_idsweep.qmd is authoritative for
# the full per-cell tables; this prints the across-cell headlines.
#
# usage (from scripts_dinamr/): Rscript idsweep_findings.R > logs/idsweep_findings.txt
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results"); REPS <- 500L; TAG <- "idsweep"
op <- options(width = 220)
cells <- read.table(file.path(SCRATCH, "idsweep.cells"), col.names = c("z", "n", "hr", "tag"),
                    colClasses = c("character", "integer", "numeric", "character"))
cells$prev <- ifelse(cells$z == "-", "12.4%", "31%")
ENG <- c("consistency", "dina", "grf")
runs <- rbind(data.frame(engine = "consistency", focus = c("effMaxSG", "effMinSG", "maxeffCons", "maxeff", "maxSG", "minSG")),
              data.frame(engine = "dina", focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")),
              data.frame(engine = "grf",  focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")))
runs$crit <- ifelse(runs$engine != "consistency" & runs$focus == "maxeffCons", "eff", runs$focus)
CRIT <- c("effMaxSG", "effMinSG", "maxeffCons", "maxeff", "eff", "maxSG", "minSG")
path <- function(e, f, ci) file.path(RES, sprintf("%s_%s_fb_mr_field_m1_h%03d_knoise0_n%d%s%s_nomr_%s_res_1_%d.rds",
  if (e == "consistency") "fs" else e, forestsearch::fs_focus_tag(e, f), round(100 * cells$hr[ci]), cells$n[ci],
  if (cells$z[ci] == "-") "" else "_z1q60", if (f %in% c("effMaxSG", "effMinSG")) "_nb20" else "", TAG, REPS))
wilson <- function(x, n) { z <- qnorm(.975); p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }

B <- list(); S <- list()
for (ci in seq_len(nrow(cells))) for (j in seq_len(nrow(runs))) {
  p <- path(runs$engine[j], runs$focus[j], ci); if (!file.exists(p)) next
  r <- readRDS(p)$results; B[[paste(ci, runs$engine[j], runs$crit[j])]] <- r
  d <- r$status %in% "DETECTED"; nd <- sum(d); w <- wilson(nd, nrow(r))
  m <- function(cc) mean(r[[cc]][d], na.rm = TRUE)
  S[[length(S) + 1L]] <- data.frame(ci = ci, cell = cells$tag[ci], prev = cells$prev[ci], hr = cells$hr[ci], n = cells$n[ci],
    engine = runs$engine[j], crit = runs$crit[j], det = nd, rate = nd / nrow(r), lo = w[1], hi = w[2],
    sens = m("sens"), spec = m("spec"), ppv = m("ppv"), npv = m("npv"),
    ratio = m("n_sel") / m("n_true"), und_npv = sum(d & is.na(r$npv)), whole = sum(d & r$n_sel == cells$n[ci], na.rm = TRUE))
}
S <- do.call(rbind, S)
done <- sort(unique(S$ci[ave(S$ci, S$ci, FUN = length) == 16]))
cat(sprintf("=== idsweep findings: %d complete cells of 18 (%s) ===\n\n", length(done), paste(cells$tag[done], collapse = ", ")))
S <- S[S$ci %in% done, ]
S$crit <- factor(S$crit, levels = CRIT); S$engine <- factor(S$engine, levels = ENG)
S$hrlab <- ifelse(S$hr == 1, "HR 1.00 (selection)", sprintf("HR %.2f (detection)", S$hr))

f3 <- function(x) sprintf("%.3f", x)
rng <- function(x) if (all(is.na(x))) "NA" else sprintf("%.3f-%.3f", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
cat("--- 1. RATE by engine x criterion, by HR block: range over cells, and by n (mean over prevalence x HR in the block) ---\n")
for (hb in unique(S$hrlab)) {
  s <- S[S$hrlab == hb, ]
  a <- aggregate(rate ~ engine + crit, s, rng); names(a)[3] <- "range_over_cells"
  for (nn in c(500, 1000, 1500)) { x <- s[s$n == nn, ]; if (nrow(x)) {
    b <- aggregate(rate ~ engine + crit, x, function(v) f3(mean(v))); a[[paste0("mean_n", nn)]] <- b$rate[match(paste(a$engine, a$crit), paste(b$engine, b$crit))] } }
  cat("\n", hb, "\n"); print(a[order(a$engine, a$crit), ], row.names = FALSE)
}
cat("\n--- 1b. Is detection flat across criteria within engine?  max - min rate over criteria, per engine x cell ---\n")
fl <- aggregate(rate ~ engine + cell, S, function(v) max(v) - min(v)); names(fl)[3] <- "spread"
fl$which_max <- vapply(seq_len(nrow(fl)), function(i) { x <- S[S$engine == fl$engine[i] & S$cell == fl$cell[i], ]
  paste(as.character(x$crit[x$rate == max(x$rate)]), collapse = "/") }, "")
print(fl[order(fl$engine, match(fl$cell, cells$tag)), ], row.names = FALSE)

cat("\n--- 1c. RATE by engine x prevalence x HR against n (one cell per entry; flat across criteria except consistency maxeff, shown separately) ---\n")
S$crit_grp <- ifelse(S$engine == "consistency" & S$crit == "maxeff", "maxeff", "other criteria")
r1 <- aggregate(rate ~ engine + crit_grp + prev + hr + n, S, function(v) sprintf("%.3f", mean(v)))
r1w <- reshape(r1, idvar = c("engine", "crit_grp", "prev", "hr"), timevar = "n", direction = "wide")
print(r1w[order(r1w$engine, r1w$crit_grp, r1w$prev, r1w$hr), ], row.names = FALSE)
cat("  (for 'other criteria' the value is identical across criteria wherever 1b shows spread 0)\n")

cat("\n--- 2b. SIZE ratio by engine x criterion x prevalence, harm cells (HR 1.50 and 1.75 averaged) and HR 1.00, against n ---\n")
S$block <- ifelse(S$hr == 1, "HR 1.00", "HR 1.50/1.75")
for (bl in c("HR 1.50/1.75", "HR 1.00")) {
  x <- S[S$block == bl, ]; if (!nrow(x)) next
  a2 <- aggregate(ratio ~ engine + crit + prev + n, x, function(v) sprintf("%.2f", mean(v)))
  w2 <- reshape(a2, idvar = c("engine", "crit", "prev"), timevar = "n", direction = "wide")
  cat("\n", bl, "\n"); print(w2[order(w2$engine, w2$crit, w2$prev), ], row.names = FALSE) }

cat("\n--- 3b. CLASSIFICATION (sens / spec / PPV / NPV) by engine x criterion x prevalence, harm cells (HR 1.50 and 1.75 averaged) and HR 1.00, against n ---\n")
for (bl in c("HR 1.50/1.75", "HR 1.00")) {
  x <- S[S$block == bl, ]; if (!nrow(x)) next
  a3 <- aggregate(cbind(sens, spec, ppv, npv) ~ engine + crit + prev + n, x, mean, na.action = na.pass)
  a3$v <- sprintf("%.2f/%.2f/%.2f/%.2f", a3$sens, a3$spec, a3$ppv, a3$npv)
  w3 <- reshape(a3[, c("engine", "crit", "prev", "n", "v")], idvar = c("engine", "crit", "prev"), timevar = "n", direction = "wide")
  cat("\n", bl, "\n"); print(w3[order(w3$engine, w3$crit, w3$prev), ], row.names = FALSE) }

cat("\n--- 2. SIZE: mean|Hhat|/mean|H| by engine x criterion; by n (mean over cells at that n) and range over cells ---\n")
a <- aggregate(ratio ~ engine + crit, S, function(v) sprintf("%.2f-%.2f", min(v), max(v))); names(a)[3] <- "range"
for (nn in c(500, 1000, 1500)) { x <- S[S$n == nn, ]; if (nrow(x)) {
  b <- aggregate(ratio ~ engine + crit, x, function(v) sprintf("%.2f", mean(v))); a[[paste0("n", nn)]] <- b$ratio[match(paste(a$engine, a$crit), paste(b$engine, b$crit))] } }
print(a[order(a$engine, a$crit), ], row.names = FALSE)
cat("size ordering per engine x cell (largest first):\n")
ord <- aggregate(ratio ~ engine + cell, S, length)
ord$order <- vapply(seq_len(nrow(ord)), function(i) { x <- S[S$engine == ord$engine[i] & S$cell == ord$cell[i], ]
  paste(as.character(x$crit[order(-x$ratio)]), collapse = " > ") }, "")
tab <- table(paste(ord$engine, ord$order)); print(as.data.frame(tab, responseName = "cells"), row.names = FALSE)

cat("\n--- 3. CLASSIFICATION by engine x criterion: mean over cells at each n (sens / spec / PPV / NPV) ---\n")
for (nn in c(500, 1000, 1500)) { x <- S[S$n == nn, ]; if (!nrow(x)) next
  b <- aggregate(cbind(sens, spec, ppv, npv) ~ engine + crit, x, mean, na.action = na.pass)
  cat(sprintf("\nn %d (over %d cells)\n", nn, length(unique(x$ci))))
  print(transform(b[order(b$engine, b$crit), ], sens = f3(sens), spec = f3(spec), ppv = f3(ppv), npv = f3(npv)), row.names = FALSE) }
cat("\nrange over cells:\n")
b <- aggregate(cbind(sens, spec, ppv, npv) ~ engine + crit, S, rng, na.action = na.pass)
print(b[order(b$engine, b$crit), ], row.names = FALSE)
cat("\nundefined NPV (selection = whole trial), detected replicates:\n")
u <- S[S$und_npv > 0 | S$whole > 0, c("cell", "engine", "crit", "det", "und_npv", "whole")]
if (nrow(u)) print(u, row.names = FALSE) else cat("  none\n")

cat("\n--- 4. CRITERION AGREEMENT: identical rule on k of m jointly detected replicates ---\n")
AG <- list()
for (ci in done) for (e in ENG) {
  ks <- runs$crit[runs$engine == e]; cmb <- combn(length(ks), 2)
  for (q in seq_len(ncol(cmb))) { a <- B[[paste(ci, e, ks[cmb[1, q]])]]; b <- B[[paste(ci, e, ks[cmb[2, q]])]]
    both <- a$status %in% "DETECTED" & b$status %in% "DETECTED"; m <- sum(both); k <- sum(both & a$sg_def == b$sg_def, na.rm = TRUE)
    w <- wilson(k, m)
    AG[[length(AG) + 1L]] <- data.frame(cell = cells$tag[ci], n = cells$n[ci], engine = e, pair = paste(ks[cmb[1, q]], ks[cmb[2, q]], sep = "/"),
                                        k = k, m = m, share = k / m, lo = w[1], hi = w[2]) } }
AG <- do.call(rbind, AG)
for (pr in c("maxeffCons/maxeff", "effMinSG/minSG", "effMinSG/eff", "effMaxSG/maxSG", "effMaxSG/eff", "minSG/eff")) {
  x <- AG[AG$pair == pr, ]; if (!nrow(x)) next
  cat(sprintf("\n%s\n", pr))
  print(transform(x[order(match(x$engine, ENG), match(x$cell, cells$tag)), c("engine", "cell", "k", "m", "share", "lo", "hi")],
                  share = f3(share), lo = f3(lo), hi = f3(hi)), row.names = FALSE) }
cat("\nshare identical, per engine x pair: median [min, max] over cells\n")
sm <- aggregate(share ~ engine + pair, AG, function(v) sprintf("%.3f [%.3f, %.3f]", median(v), min(v), max(v)))
print(sm[order(match(sm$engine, ENG), -as.numeric(sub(" .*", "", sm$share))), ], row.names = FALSE)

cat("\n--- 5. Gate I record: realized prevalence, same draws, undefined rates ---\n")
gi <- readRDS(file.path(SCRATCH, "idsweep_gateI.rds"))
for (k in cells$tag[cells$tag %in% names(gi)]) { g <- gi[[k]]
  sd <- g$same_draws
  cat(sprintf("%-16s gate %d/%d  prevalence realized %s (super %s)  within-cell n_true identical %s  same-draws vs committed: %s\n", k,
              g$pass, g$fail, paste(unique(sprintf("%.4f", g$per_run$prev_realized)), collapse = ","),
              sprintf("%.4f", g$per_run$prev_super[1]), g$within_n_true,
              if (is.null(sd)) "no committed bundle covering 1-500" else
                paste(sprintf("%s=%s", sd$campaign, sd$identical), collapse = " "))) }
saveRDS(list(per_run = S, agreement = AG), file.path(SCRATCH, "idsweep_findings.rds"))
options(op)
