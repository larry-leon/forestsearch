# declcal_c0approx_findings.R -- read-only post-processing for
# TASK_declcal_c0_approx_2026-09-22 (sections 2-4).  Run from the gbsg_020
# directory with C0A_CELLDIR pointing at the six per-cell capture payloads
# written by declcal_c0approx.sh.
#   1. combine the captures (reps 1..C0A_MAXREP per cell) -> results/declcal_c0approx_res.rds
#   2. identity gate vs the committed declcal_inull payloads (at 81752681):
#      max_T_pre and G_pre identical (STOP, exit 6, otherwise); kappa_hat_05
#      B 2000 - B 500 differences reported
#   3. per n and c0: median / min / IQR / max of kappa_hat_{05,10}_c0 over the
#      captures at that n; each median applied as a FIXED cutoff to max_T_pre
#      of all 2000 replicates of every committed B and C cell at that n
#   4. tables -> scripts_dinamr/logs/declcal_c0approx.txt
suppressPackageStartupMessages(library(data.table))
celldir <- Sys.getenv("C0A_CELLDIR")
stopifnot(nzchar(celldir), dir.exists(celldir))
pin <- "81752681"
maxrep  <- as.integer(Sys.getenv("C0A_MAXREP", "20"))        # captures used: reps 1..maxrep per cell
out_rds <- Sys.getenv("C0A_OUT_RDS", "results/declcal_c0approx_res.rds")
out_txt <- Sys.getenv("C0A_OUT_TXT", "scripts_dinamr/logs/declcal_c0approx.txt")

B_cells <- sprintf("B%d", 1:6); C_cells <- sprintf("C%d", 1:4)
c0_grid <- c(0.70, 0.75, 0.80, 0.85)
c0_sfx  <- sprintf("c%03d", as.integer(round(100 * c0_grid)))
ref_path <- function(cell) sprintf("results/declcal_%s_%s_res_1_2000.rds",
                                   if (substr(cell, 1, 1) == "B") "inull" else "power", cell)

# committed payloads must be the ones at the pin
ref_files <- vapply(c(B_cells, C_cells), ref_path, "")
dirty <- system2("git", c("diff", "--name-only", pin, "--", ref_files), stdout = TRUE)
if (length(dirty)) stop("committed declcal payloads differ from ", pin, ": ", paste(dirty, collapse = ", "))
ref <- lapply(setNames(ref_files, c(B_cells, C_cells)), function(p) readRDS(p)$results)

# ---- 1. combine ---------------------------------------------------------------
pl <- lapply(setNames(B_cells, B_cells), function(cl)
  readRDS(file.path(celldir, sprintf("declcal_c0approx_%s_res_1_20.rds", cl))))
res <- rbindlist(lapply(pl, function(p) as.data.table(p$results)), fill = TRUE)
aux <- rbindlist(lapply(names(pl), function(cl) cbind(cell_id = cl, as.data.table(pl[[cl]]$aux))), fill = TRUE)
res <- res[rep <= maxrep]; aux <- aux[rep <= maxrep]
setorder(res, n, cell_id, rep)
stopifnot(nrow(res) == 6L * maxrep, all(sort(unique(res$rep)) == seq_len(maxrep)), all(res$status == "ok"), all(res$B_cal == 2000L),
          all(vapply(pl, function(p) isTRUE(p$meta$final) && p$meta$cell_status == "complete", TRUE)))

L <- character(0); say <- function(...) { s <- sprintf(...); L <<- c(L, s); cat(s, "\n", sep = "") }
f4 <- function(x) sprintf("%.4f", x)
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d; hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(ctr - hw, ctr + hw)
}
rate_str <- function(ind) {
  ind[is.na(ind)] <- 0L; x <- sum(ind); n <- length(ind); w <- wilson(x, n)
  sprintf("%s [%s, %s]", f4(x / n), f4(w[1]), f4(w[2]))
}

say("declcal c0 approximate table -- plug-in fixed-cutoff rates from a per-n median kappa_hat(c0), %d field captures", nrow(res))
say("task: dev/tasks/TASK_declcal_c0_approx_2026-09-22.md ; generated %s", format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
say("captures: reps 1-%d of B1..B6 (%d in all),", maxrep, nrow(res))
say("   B = 2000, centred Poisson multipliers, c1 = c2 = 1.0, template floors, declaration_c0 = (0.70, 0.75, 0.80, 0.85)")
say("forestsearch %s built %s ; git %s", pl[[1]]$meta$forestsearch_version, pl[[1]]$meta$forestsearch_built,
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE))
say("committed comparators: results/declcal_{inull,power}_*_res_1_2000.rds at %s (B_cal 500), read-only", pin)
say("")
say("THESE ARE PLUG-IN FIXED-CUTOFF RATES: rate = mean(max_T_pre >= median kappa_hat(c0) at that n), not the")
say("per-replicate calibrated rule (which uses each replicate's own kappa_hat). The campaign task supplies the exact version.")
say("median over the %d captures at each n (%d from each of the two B cells at that n); per-cell medians in section 1b.", 2L * maxrep, maxrep)

# ---- 2. identity gate ------------------------------------------------------------
say(""); say("## 0. Identity gate (task section 2)")
say("")
say("| cell | n | reps | max_T_pre identical | G_pre identical | kappa_hat_05 diff (B2000 - B500): mean | min | max | max abs |")
say("|---|---|---|---|---|---|---|---|---|")
id_fail <- character(0); all_d <- numeric(0)
for (cl in B_cells) {
  r <- res[cell_id == cl]; m <- ref[[cl]][match(r$rep, ref[[cl]]$rep), ]
  okT <- identical(r$max_T_pre, m$max_T_pre) && all(m$status == "ok")
  okG <- identical(as.integer(r$G_pre), as.integer(m$G_pre))
  if (!okT) id_fail <- c(id_fail, sprintf("%s max_T_pre reps %s", cl, paste(r$rep[r$max_T_pre != m$max_T_pre], collapse = ",")))
  if (!okG) id_fail <- c(id_fail, sprintf("%s G_pre reps %s", cl, paste(r$rep[r$G_pre != m$G_pre], collapse = ",")))
  d <- r$kappa_hat_05 - m$kappa_hat_05; all_d <- c(all_d, d)
  say("| %s | %d | %d | %s | %s | %+.4f | %+.4f | %+.4f | %.4f |", cl, r$n[1], nrow(r), okT, okG,
      mean(d), min(d), max(d), max(abs(d)))
}
say("")
say("kappa_hat_05 B2000 - B500 over all %d: mean %+.4f, SD %.4f, max |diff| %.4f ; |diff| < 0.15 on %d of %d",
    length(all_d), mean(all_d), sd(all_d), max(abs(all_d)), sum(abs(all_d) < 0.15), length(all_d))
say("per-draw monotone checks not part of this task; fidelity gate (declared_conv vs search) per cell log: all passed (exit 0)")
if (length(id_fail)) {
  say("IDENTITY GATE: FAIL -- %s", paste(id_fail, collapse = " ; ")); say("STOP.")
  writeLines(L, out_txt); quit(status = 6L)
}
say("IDENTITY GATE: PASS (max_T_pre and G_pre identical on all %d replicates)", nrow(res))

# ---- 3. approximation --------------------------------------------------------------
lev <- data.table(lab = c(sprintf("%.2f", c0_grid), "c2 (1.00, unshifted)"),
                  sfx = c(paste0("_", c0_sfx), ""), c0 = c(c0_grid, 1.0))
alphas <- c("05" = 0.05, "10" = 0.10)
qs <- function(x) c(med = median(x), min = min(x), q1 = unname(quantile(x, 0.25)),
                    q3 = unname(quantile(x, 0.75)), max = max(x))
med <- list()   # key "n|lab|a" -> median kappa
say(""); say("## 1. Table 1 -- median kappa_hat(c0) per n over the captures (B = 2000)")
say("")
say("| n | c0 | median k05 | min | IQR | max | implied p* (median k05) | median k10 | min | IQR | max | mean fw_1645 |")
say("|---|---|---|---|---|---|---|---|---|---|---|---|")
for (nn in sort(unique(res$n))) for (i in seq_len(nrow(lev))) {
  r <- res[n == nn]; s <- lev$sfx[i]
  a <- qs(r[[paste0("kappa_hat_05", s)]]); b <- qs(r[[paste0("kappa_hat_10", s)]])
  fw <- if (s == "") r$alpha_FW_hat_1645 else r[[paste0("fw_1645", s)]]
  med[[sprintf("%d|%s|05", nn, lev$lab[i])]] <- a[["med"]]
  med[[sprintf("%d|%s|10", nn, lev$lab[i])]] <- b[["med"]]
  say("| %d | %s | %.4f | %.4f | %.4f-%.4f | %.4f | %.4f | %.4f | %.4f | %.4f-%.4f | %.4f | %.4f |",
      nn, lev$lab[i], a[["med"]], a[["min"]], a[["q1"]], a[["q3"]], a[["max"]], 2 * pnorm(a[["med"]]) - 1,
      b[["med"]], b[["min"]], b[["q1"]], b[["q3"]], b[["max"]], mean(fw))
}
say("")
say("c2 row: the unshifted kappa_hat at B = 2000 from the same captures (reference; not used below).")
say("fw_1645 = mean over the B draws of 1{M*(c0) > 1.6449}: the family-wise size of the p* = 0.90 rule if the true null sat at c0.")

say(""); say("## 1b. Per-cell medians (%d captures each) -- sensitivity of the pooled median", maxrep)
say("")
say("| cell | n | c0 | median k05 | median k10 |"); say("|---|---|---|---|---|")
for (cl in B_cells) for (i in seq_len(nrow(lev))) {
  r <- res[cell_id == cl]; s <- lev$sfx[i]
  say("| %s | %d | %s | %.4f | %.4f |", cl, r$n[1], lev$lab[i],
      median(r[[paste0("kappa_hat_05", s)]]), median(r[[paste0("kappa_hat_10", s)]]))
}

# ---- 4. plug-in rates ----------------------------------------------------------------
cells <- c(B_cells, C_cells)
n_of <- vapply(cells, function(cl) ref[[cl]]$n[1], 0)
na_pre <- vapply(cells, function(cl) sum(is.na(ref[[cl]]$max_T_pre)), 0)
stopifnot(all(vapply(cells, function(cl) nrow(ref[[cl]]) == 2000L && all(ref[[cl]]$status == "ok"), TRUE)))
rates <- list()
rate_x <- function(cl, cut) { m <- ref[[cl]]$max_T_pre; x <- !is.na(m) & m >= cut; x }
for (a in names(alphas)) {
  say(""); say("## 2. Table 2 -- declaration rate at alpha %.2f, Wilson 95%% (x / 2000)", alphas[[a]])
  say("")
  say("| cell | n | p* 0.90 as executed (committed) | c0 0.70 | c0 0.75 | c0 0.80 | c0 0.85 | c0 = c2 (committed, per-replicate) |")
  say("|---|---|---|---|---|---|---|---|")
  for (cl in cells) {
    R <- ref[[cl]]; nn <- n_of[[cl]]
    v <- c(rate_str(R$declared_conv),
           vapply(sprintf("%.2f", c0_grid), function(lb) {
             x <- rate_x(cl, med[[sprintf("%d|%s|%s", nn, lb, a)]])
             rates[[sprintf("%s|%s|%s", cl, lb, a)]] <<- mean(x); rate_str(x) }, ""),
           rate_str(R[[paste0("declared_cal", a)]]))
    rates[[sprintf("%s|conv|%s", cl, a)]] <- mean(R$declared_conv)
    rates[[sprintf("%s|c2|%s", cl, a)]] <- mean(R[[paste0("declared_cal", a)]])
    say("| %s | %d | %s |", cl, nn, paste(v, collapse = " | "))
  }
}
say("")
say("c0 columns: fixed cutoff = the per-n median kappa_hat(c0) of Table 1, applied to max_T_pre (pre-reduction family) of the committed replicates.")
say("p* 0.90 as executed = declared_conv (rounded rule, post-reduction family). c0 = c2 = declared_cal05 / declared_cal10 (per-replicate kappa_hat, B = 500).")
say("max_T_pre NA counts (counted as not declared): %s", paste(sprintf("%s %d", cells, na_pre), collapse = ", "))

for (cl in cells) rates[[sprintf("%s|k2", cl)]] <- mean(rate_x(cl, 2.0))
worst <- function(key) { v <- vapply(B_cells, function(cl) rates[[sprintf(key, cl)]], 0); sprintf("%s (%s)", f4(max(v)), B_cells[which.max(v)]) }
for (a in names(alphas)) {
  say(""); say("## 3. Table 3 -- worst uniform-benefit rate vs power, alpha %.2f", alphas[[a]])
  say("")
  say("| rule | worst B rate (cell) | HR 1.5 n 1000 (C1) | HR 1.5 n 1500 (C2) | HR 2.0 n 1000 (C3) | HR 2.0 n 1500 (C4) |")
  say("|---|---|---|---|---|---|")
  row <- function(lab, key) say("| %s | %s | %s |", lab, worst(key),
                                paste(vapply(C_cells, function(cl) f4(rates[[sprintf(key, cl)]]), ""), collapse = " | "))
  row("fixed p* 0.9545 (k 2.0, committed, max_T_pre)", "%s|k2")
  for (lb in sprintf("%.2f", c0_grid)) row(sprintf("c0 %s (plug-in median kappa)", lb), paste0("%s|", lb, "|", a))
  row("c0 = c2 (committed, per-replicate kappa_hat)", paste0("%s|c2|", a))
}
say("")
say("fixed k 2.0 row: max_T_pre >= 2.0 on the committed payloads (declcal_fixedk_practical.txt, max_T_pre column).")

# ---- payload -------------------------------------------------------------------------
payload <- list(
  results = as.data.frame(res), aux = as.data.frame(aux),
  approx = list(median_kappa = unlist(med), rates = unlist(rates), kappa05_diff_B2000_B500 = all_d),
  meta = list(task = "TASK_declcal_c0_approx_2026-09-22", pin = pin, maxrep = maxrep,
              declaration_c0 = c0_grid, c0_suffix = c0_sfx, B_cal = 2000L,
              cells = lapply(pl, `[[`, "meta"), written_at = Sys.time()))
saveRDS(payload, out_rds)
writeLines(L, out_txt)
cat(sprintf("wrote %s and %s\n", out_rds, out_txt))
