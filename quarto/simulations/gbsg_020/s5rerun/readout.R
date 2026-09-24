# Section 5 full re-run read-out (TASK_section5_full_rerun_2026-09-24 section 4).
# Per cell: s5rerun combined bundle (after) against the committed published bundle (before),
# paired by sim_id.  Shifts via s5rerun/pair.R (a copy of mrs5sweep/pair.R); coverage via the
# package's fs_sim_bias_coverage() with the summary documents' subst_s() for field-s.
# usage (from gbsg_020/): Rscript s5rerun/readout.R
suppressMessages(library(forestsearch))
G <- "."; D <- file.path(G, "s5rerun")
cells <- read.table(file.path(G, "mrs5sweep/cells.txt"), comment.char = "#",
  col.names = c("id","src","hr","n","z1q","cw","cwk","stem"), stringsAsFactors = FALSE)
subst_s <- function(r) { for (s in c("est2","up1s","lo1s","lo2s","hi2s","lo_se","hi_se","se","lam_mean"))
  r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
covv <- function(r, block, side) { rr <- if (block == "Hc") subst_s(r) else r
  t <- fs_sim_bias_coverage(rr, block = block, estimators = "fld", side = side); t$cov1[1] }
f4 <- function(v) sprintf("%+.5f", v)
rows <- list(); covrows <- list(); notrun <- character(0)
for (i in seq_len(nrow(cells))) {
  c1 <- cells[i, ]; zt <- if (c1$z1q != "-") "_z1q60" else ""
  af <- sprintf("results/fs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_s5rerun_combined_1_2000.rds",
                round(100 * c1$hr), c1$n, zt)
  if (!file.exists(af)) { notrun <- c(notrun, c1$id); next }
  pf <- file.path(D, "pair", paste0(c1$id, ".rds"))
  st <- system2("Rscript", c(file.path(D, "pair.R"), c1$id, af, c1$stem, pf), stdout = TRUE, stderr = TRUE)
  p <- readRDS(pf); s <- p$stats; g <- function(q, k) s[s$qty == q, k]
  A <- readRDS(af)$results; B <- readRDS(c1$stem)$results; B <- B[match(A$sim_id, B$sim_id), ]
  d <- A$detected == 1L
  ind <- list(lo = function(r) r$fld_H_lo1s <= r$betaHhat_H, up = function(r) r$fld_Hc_up1s_s >= r$betaHhat_Hc,
              bonf_lo = function(r) r$fld_joint_s_bonf_loH <= r$betaHhat_H,
              bonf_up = function(r) r$fld_joint_s_bonf_upHc >= r$betaHhat_Hc)
  flipcov <- vapply(ind, function(f) sum(d & (f(A) != f(B)), na.rm = TRUE), 0)
  rows[[c1$id]] <- data.frame(cell = c1$id, src = c1$src, hr = c1$hr, n = c1$n,
    decl = sprintf("%.4f / %.4f", p$det_rate_before, p$det_rate_after), decl_same = p$det_same,
    est = paste(f4(c(g("est","mean"), g("est","median"), g("est","p05"), g("est","p95"))), collapse = " / "),
    lo  = paste(f4(c(g("lo","mean"),  g("lo","median"),  g("lo","p05"),  g("lo","p95"))),  collapse = " / "),
    up  = paste(f4(c(g("up","mean"),  g("up","median"),  g("up","p05"),  g("up","p95"))),  collapse = " / "),
    region_chg = p$region_changed,
    cross = sum(s$flips075[s$qty %in% c("lo","up","bonf_lo","bonf_up")] + s$flips125[s$qty %in% c("lo","up","bonf_lo","bonf_up")]),
    cov_flips = sum(flipcov), errors = p$errors, stringsAsFactors = FALSE)
  covrows[[c1$id]] <- data.frame(cell = c1$id,
    lo_before = covv(B, "H", "lower"), lo_after = covv(A, "H", "lower"),
    up_before = covv(B, "Hc", "upper"), up_after = covv(A, "Hc", "upper"),
    flips_lo = flipcov["lo"], flips_up = flipcov["up"], flips_bonf_lo = flipcov["bonf_lo"], flips_bonf_up = flipcov["bonf_up"])
}
T <- do.call(rbind, rows); C <- do.call(rbind, covrows)
saveRDS(list(table = T, coverage = C, not_run = notrun), file.path(D, "readout.rds"))
hdr <- c("| cell | source | HR | n | decl. rate before / after | est shift mean / med / p05 / p95 | field lower shift mean / med / p05 / p95 | field-s upper shift mean / med / p05 / p95 | region chg | bound crossings 0.75/1.25 | coverage-indicator changes | errors |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|")
md <- c(hdr, sprintf("| %s | %s | %.2f | %d | %s | %s | %s | %s | %d | %d | %d | %d |", T$cell, T$src, T$hr, T$n, T$decl,
  T$est, T$lo, T$up, T$region_chg, T$cross, T$cov_flips, T$errors), "",
  "| cell | field lower cov. before | after | field-s upper cov. before | after |", "|---|---|---|---|---|",
  sprintf("| %s | %.4f | %.4f | %.4f | %.4f |", C$cell, C$lo_before, C$lo_after, C$up_before, C$up_after))
if (length(notrun)) md <- c(md, "", paste("Not run / not complete:", paste(notrun, collapse = ", ")))
writeLines(md, file.path(D, "readout.md")); cat(md, sep = "\n")
