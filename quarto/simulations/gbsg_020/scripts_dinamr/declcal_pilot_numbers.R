# declcal_pilot_numbers.R -- Stage 1 read-out (TASK_declcal_CAMPAIGN_2026-09-22_v2,
# section 8): measured clock, B_cal rule, projection, per-replicate cap and the
# campaign hard cap.  usage: Rscript declcal_pilot_numbers.R <ceiling_hours> <block_C 0|1>
# Writes logs/declcal_pilot_numbers.txt and declcal_pilot_numbers.rds; prints
# "DECISION: CONTINUE" or "DECISION: STOP" and the Stage 2 settings.
args <- commandArgs(trailingOnly = TRUE)
ceiling_h <- as.numeric(args[1]); block_c <- as.integer(args[2])
p <- readRDS("../results/declcal_pilot_A2_res_1_200.rds")
r <- p$results; a <- p$aux[match(p$results$rep, p$aux$rep), ]
ok <- r$status == "ok"
out <- character(0)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
qs <- function(x) paste(sprintf("%.2f", stats::quantile(x, c(0, .25, .5, .75, 1), na.rm = TRUE)), collapse = " / ")
say("declcal pilot read-out -- cell A2 (complete null, n = 1000), %d replicates, B 2000 assembled once and sub-sampled", nrow(r))
say("status: %s", paste(names(table(r$status)), table(r$status), collapse = ", "))
say("workers %d ; cell elapsed %.0f s ; host %s ; forestsearch %s built %s",
    p$meta$n_workers, p$meta$elapsed_s, p$meta$hostname, p$meta$forestsearch_version, p$meta$forestsearch_built)
say("per-replicate wall (s), min / q25 / median / q75 / max: total %s", qs(r$wall_sec[ok]))
say("                                              search %s", qs(a$wall_search[ok]))
say("                                              field  %s", qs(a$wall_field[ok]))
say("G_pre min / q25 / median / q75 / max: %s", qs(r$G_pre[ok]))
say("G_post min / q25 / median / q75 / max: %s", qs(r$G_post[ok]))
say("fidelity: declared_conv vs search indicator: %d of %d agree", sum(r$declared_conv[ok] == a$search_declared[ok]), sum(ok))
say("realized conventional declaration fraction (as executed): %.4f ; exact z: %.4f", mean(r$declared_conv[ok]), mean(r$declared_conv_exact[ok]))
rates <- sapply(c(500L, 1000L, 2000L), function(B) c(
  cal05 = mean(a[[sprintf("cal05_B%d", B)]][ok]), cal10 = mean(a[[sprintf("cal10_B%d", B)]][ok]),
  k05 = median(a[[sprintf("kappa05_B%d", B)]][ok])))
colnames(rates) <- c("B500", "B1000", "B2000")
for (B in colnames(rates))
  say("  %-6s calibrated declaration rate alpha 0.05 %.4f ; alpha 0.10 %.4f ; median kappa_hat_05 %.4f",
      B, rates["cal05", B], rates["cal10", B], rates["k05", B])
# B_cal rule: smallest B whose alpha = 0.05 calibrated rate is within 0.005 of B = 2000's.
ref <- rates["cal05", "B2000"]
cand <- c(500L, 1000L)[abs(rates["cal05", c("B500", "B1000")] - ref) <= 0.005]
B_cal <- if (length(cand)) min(cand) else 2000L
say("B_cal rule (alpha 0.05 rate within 0.005 of B 2000's %.4f): |diff| B500 %.4f, B1000 %.4f -> B_cal = %d",
    ref, abs(rates["cal05", "B500"] - ref), abs(rates["cal05", "B1000"] - ref), B_cal)
med <- median(r$wall_sec[ok]); W <- p$meta$n_workers
N <- if (block_c == 1L) 26000L else 18000L
proj_s <- med * N / W
say("projection: median per-replicate wall %.3f s x %d replicates / %d workers = %.0f s = %.2f h (ceiling %.1f h)",
    med, N, W, proj_s, proj_s / 3600, ceiling_h)
say("  (assumption: the pilot ran with every worker busy on the same cell type, so contention is in the median; the pilot's B is 2000, so B_cal < 2000 makes this conservative on the field share)")
thr_obs <- p$meta$elapsed_s / nrow(r)
say("  cross-check from the pilot's own throughput: %.3f s of cell wall per replicate x %d = %.0f s = %.2f h (includes the ragged last round of 200 / %d)",
    thr_obs, N, thr_obs * N, thr_obs * N / 3600, W)
cap <- 10 * med; hard <- 1.5 * proj_s
say("per-replicate hard cap: 10 x %.3f = %.1f s ; campaign hard cap: 1.5 x %.0f = %.0f s (%.2f h)", med, cap, proj_s, hard, hard / 3600)
dec <- if (proj_s / 3600 <= ceiling_h) "CONTINUE" else "STOP"
say("DECISION: %s (projection %.2f h vs ceiling %.1f h)", dec, proj_s / 3600, ceiling_h)
writeLines(out, "logs/declcal_pilot_numbers.txt")
saveRDS(list(rates = rates, B_cal = B_cal, median_wall = med, workers = W, N = N,
             projection_s = proj_s, cap_s = cap, hardcap_s = hard, decision = dec,
             ceiling_h = ceiling_h, block_c = block_c,
             wall = list(total = r$wall_sec[ok], search = a$wall_search[ok], field = a$wall_field[ok]),
             G_pre = r$G_pre[ok], G_post = r$G_post[ok],
             conv_rate = mean(r$declared_conv[ok])), "declcal_pilot_numbers.rds")
