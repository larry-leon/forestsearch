#!/usr/bin/env Rscript
# nullmr Step 3.4 -- grid projection from the two smoke cells, n-scaling MEASURED.
# Per identifier, the per-replicate cost fit_mr_secs (search + MR) is split by
# declaration: MR runs only on declaring replicates, so
#   cost(cell, engine) = (1 - d) * s0(n) + d * s1(n)
# with s0 / s1 the mean non-declaring / declaring cost measured at n 500
# (null0657_n500 smoke) and n 1500 (null0721_n1500 smoke), linear in n between
# (n 1000 = midpoint), and d the cell's declaration rate from the committed
# nullid bundle (identification is nullid's; the identity gate).  A class with
# no smoke replicate at an n borrows the other n's pooled cost, scaled by the
# measured ratio of all-replicate means.  Load factor at 64 workers: memory
# 'contention-factor-per-n' (63w 1.17 at n500, 3.27 at n2000, linear in n), as
# nullc125 used; render overhead = the smoke renders' wall minus their compute.
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)
W <- 64L
ENG <- c(consistency = "fs", dina = "dina", grf = "grf")
sm <- function(tag, hr, n) readRDS(sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nullmrsmoke_quickrun_res_1_20.rds",
                                           tag, round(100 * hr), n, round(1000 * hr)))$results
wall <- function(cell, e) { l <- readLines(sprintf("scripts_dinamr/logs/nullmrsmoke_%s_%s_effMaxSG.log", cell, e))
  as.numeric(sub("^WALL_SECONDS=([0-9]+).*", "\\1", grep("^WALL_SECONDS=", l, value = TRUE))) }
cat("nullmr Step 3.1 smoke and the grid projection\n")
cat(sprintf("TASK_null_gbsg_mr_2026-09-21 ; %s ; host %s ; R %s ; forestsearch %s ; %d workers\n\n",
            format(Sys.time()), Sys.info()[["nodename"]], getRversion(), utils::packageVersion("forestsearch"), W))
S <- list(); ovh <- c()
for (e in names(ENG)) {
  a <- sm(ENG[[e]], 0.657, 500L); b <- sm(ENG[[e]], 0.721, 1500L)
  f <- function(r, d) { x <- r$fit_mr_secs[r$detected == d]; if (length(x)) mean(x) else NA_real_ }
  S[[e]] <- data.frame(n = c(500, 1500), s0 = c(f(a, 0L), f(b, 0L)), s1 = c(f(a, 1L), f(b, 1L)),
                       all = c(mean(a$fit_mr_secs), mean(b$fit_mr_secs)),
                       nd = c(sum(a$detected), sum(b$detected)))
  wa <- wall("null0657_n500", e); wb <- wall("null0721_n1500", e)
  ovh <- c(ovh, wa - max(a$fit_mr_secs), wb - max(b$fit_mr_secs))
  cat(sprintf("== %s ==\n", e))
  cat(sprintf("  n 500  (null0657): declared %2d/20 ; mean s non-declaring %s ; declaring %s ; all %.3f ; render wall %d s\n",
              S[[e]]$nd[1], format(round(S[[e]]$s0[1], 3)), format(round(S[[e]]$s1[1], 3)), S[[e]]$all[1], wa))
  cat(sprintf("  n 1500 (null0721): declared %2d/20 ; mean s non-declaring %s ; declaring %s ; all %.3f ; render wall %d s\n",
              S[[e]]$nd[2], format(round(S[[e]]$s0[2], 3)), format(round(S[[e]]$s1[2], 3)), S[[e]]$all[2], wb))
  cat(sprintf("  measured n-scaling of the all-replicate mean, n1500 / n500: %.3f\n", S[[e]]$all[2] / S[[e]]$all[1]))
}
ovh_s <- max(ovh)
cat(sprintf("\n  render overhead (wall - slowest replicate), max over the 6 smoke renders: %.0f s\n\n", ovh_s))
fill <- function(v, all) { if (is.na(v[1])) v[1] <- v[2] * all[1] / all[2]; if (is.na(v[2])) v[2] <- v[1] * all[2] / all[1]; v }
at_n <- function(v, n) v[1] + (v[2] - v[1]) * (n - 500) / 1000
lf <- function(n) 1.17 + (3.27 - 1.17) * (n - 500) / 1500
CELLS <- data.frame(cell = c("null0657_n500","null0721_n500","null0657_n1000","null0721_n1000","null0657_n1500","null0721_n1500"),
                    hr = c(0.657, 0.721, 0.657, 0.721, 0.657, 0.721), n = c(500, 500, 1000, 1000, 1500, 1500))
tot <- 0; tot_nolf <- 0; tot_max <- 0
cat("| cell | engine | nullid decl. rate | s/rep | compute s @64w (x load) | render s |\n|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(CELLS))) for (e in names(ENG)) {
  p <- sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr_nullid_res_1_2000.rds",
               ENG[[e]], round(100 * CELLS$hr[i]), CELLS$n[i], round(1000 * CELLS$hr[i]))
  d <- mean(readRDS(p)$results$detected == 1L)
  s0 <- at_n(fill(S[[e]]$s0, S[[e]]$all), CELLS$n[i]); s1 <- at_n(fill(S[[e]]$s1, S[[e]]$all), CELLS$n[i])
  c1 <- ((1 - d) * s0 + d * s1) * 2000 / W
  r <- c1 * lf(CELLS$n[i]) + ovh_s
  tot <- tot + r; tot_nolf <- tot_nolf + c1 + ovh_s; tot_max <- tot_max + c1 * 3.27 + ovh_s
  cat(sprintf("| %s | %s | %.4f | %.3f | %.0f (x%.2f) | %.0f |\n", CELLS$cell[i], e, d, (1 - d) * s0 + d * s1, c1, lf(CELLS$n[i]), r))
}
cat(sprintf("\n  PROJECTED GRID: %.2f h  (no load factor %.2f h ; load factor 3.27 at every n %.2f h)\n", tot / 3600, tot_nolf / 3600, tot_max / 3600))
cat(sprintf("  bound 16 h -> %s\n", if (tot / 3600 < 16) "UNDER: the grid runs straight away (Step 3.4)" else "OVER: record and stop"))
