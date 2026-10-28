# Post-processing of the ladder, vignettes stage 0 (2026-10-27). No fit is run.
#   Rscript postprocess.R > postprocess_log.txt 2>&1   (from this directory)
suppressPackageStartupMessages(library(forestsearch))
options(width = 160)
rd <- function(r) readRDS(file.path("rungs", paste0("rung_", r, ".rds")))
fit_S1 <- readRDS("fits/fit_S1.rds")

cat("==== S1: top-level names of the forestsearch result ====\n")
print(names(fit_S1))
cat("\n==== S1: str(fit$mr_inference, max.level = 2) ====\n")
str(fit_S1$mr_inference, max.level = 2, list.len = 200, give.attr = FALSE)
cat("\n==== S1: names(fit$mr_inference$field), $complement, $joint_s ====\n")
print(names(fit_S1$mr_inference$field))
print(names(fit_S1$mr_inference$field$complement))
print(names(fit_S1$mr_inference$field$joint_s))
cat("\n==== S1: settings ====\n"); str(fit_S1$mr_inference$settings)
cat("\n==== S1: declaration_field present? ====\n")
dfld <- fit_S1$mr_inference$declaration_field
cat("is.null:", is.null(dfld), "\n"); if (!is.null(dfld)) str(dfld, max.level = 1, list.len = 50)
cat("\n==== fs_declaration_calibration(fit_S1, c0 = c(0.70, 0.75, 0.80, 0.85)) ====\n")
dc <- tryCatch(fs_declaration_calibration(fit_S1, c0 = c(0.70, 0.75, 0.80, 0.85)),
               error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(dc)) { str(dc, max.level = 1, list.len = 60); cat("-- print:\n"); print(dc) }

cat("\n==== print(fit_S1) ====\n"); print(fit_S1)
cat("\n==== summary(fit_S1) ====\n")
s <- tryCatch(summary(fit_S1), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(s)) print(s)
cat("\n==== print(fit_S1$mr_inference) ====\n")
tryCatch(print(fit_S1$mr_inference), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))

cat("\n==== S1 vs S1r: all.equal on region and every Table 1 quantity ====\n")
a <- rd("S1"); b <- rd("S1r")
keys <- c("naive_H", "naive_lo_H", "naive_hi_H", "corrected_H", "field_lower_H",
          "naive_Hc", "naive_lo_Hc", "naive_hi_Hc", "corrected_Hc", "field_s_upper_Hc",
          "bonf_lower_H", "bonf_upper_Hc", "n_family", "mr_draws", "field_R_out", "field_R_in")
print(all.equal(a$label, b$label)); print(all.equal(a[c("N_H","E_H","N_Hc","E_Hc")], b[c("N_H","E_H","N_Hc","E_Hc")]))
print(all.equal(a$mr[keys], b$mr[keys], tolerance = 0))
print(identical(a$mr[keys], b$mr[keys]))
cat("whole-object identical(fit_S1$mr_inference minus timings):\n")
strip <- function(x) { if (is.list(x)) { x <- x[!grepl("timing|_sec$", names(x))]; lapply(x, strip) } else x }
print(all.equal(strip(fit_S1$mr_inference), strip(readRDS("fits/fit_S1r.rds")$mr_inference), tolerance = 0))

cat("\n==== rung vs Table 1 (published, 3 dp) and vs committed payload (full precision) ====\n")
pub <- list(
  fs   = c(naive_H = 2.222, naive_lo_H = 1.169, naive_hi_H = 4.221, corrected_H = 1.276, field_lower_H = 0.634,
           naive_Hc = 0.605, naive_lo_Hc = 0.465, naive_hi_Hc = 0.788, corrected_Hc = 0.645, field_s_upper_Hc = 0.804,
           bonf_lower_H = 0.509, bonf_upper_Hc = 0.840),
  dina = c(naive_H = 1.315, naive_lo_H = 0.746, naive_hi_H = 2.317, corrected_H = 1.076, field_lower_H = 0.638,
           naive_Hc = 0.641, naive_lo_Hc = 0.489, naive_hi_Hc = 0.839, corrected_Hc = 0.658, field_s_upper_Hc = 0.842,
           bonf_lower_H = 0.572, bonf_upper_Hc = 0.887),
  grf  = c(naive_H = 1.951, naive_lo_H = 1.037, naive_hi_H = 3.673, corrected_H = 1.111, field_lower_H = 0.510,
           naive_Hc = 0.615, naive_lo_Hc = 0.472, naive_hi_Hc = 0.801, corrected_Hc = 0.663, field_s_upper_Hc = 0.846,
           bonf_lower_H = 0.404, bonf_upper_Hc = 0.874))
pl <- c(fs = "analysis_gbsg_mr", dina = "analysis_gbsg_dina_mr", grf = "analysis_gbsg_grf_mr")
payload_vals <- function(id) {
  p <- readRDS(file.path(path.expand("~/Documents/GitHub/fs-glms-interpretable/quarto/gbsg"), paste0(pl[[id]], "_payload.rds")))
  t <- p$table; H <- t[t$region == "H", ]; C <- t[t$region == "Hc", ]
  c(naive_H = H$naive_est, naive_lo_H = H$naive_lo, naive_hi_H = H$naive_hi, corrected_H = H$mr_est,
    field_lower_H = p$bounds$field_lower_H, naive_Hc = C$naive_est, naive_lo_Hc = C$naive_lo, naive_hi_Hc = C$naive_hi,
    corrected_Hc = C$mr_est, field_s_upper_Hc = p$bounds$field_s_upper_Hc,
    bonf_lower_H = p$bounds$bonf_lower_H, bonf_upper_Hc = p$bounds$bonf_upper_Hc)
}
cmp <- list()
for (r in c("S1", "S2", "S3", "D1", "D3", "G1", "G3")) {
  x <- rd(r); id <- x$identifier; v <- unlist(x$mr[names(pub[[id]])])
  pv <- payload_vals(id)
  cmp[[r]] <- data.frame(rung = r, quantity = names(v), rung_value = round(v, 6), published = pub[[id]],
                         log_diff_vs_published = round(log(v) - log(pub[[id]]), 4),
                         rung_3dp_eq_published = sprintf("%.3f", v) == sprintf("%.3f", pub[[id]]),
                         payload = round(pv, 6), abs_diff_vs_payload = signif(abs(v - pv), 3),
                         row.names = NULL)
}
cmp <- do.call(rbind, cmp); print(cmp, row.names = FALSE)
write.csv(cmp, "comparison_table1.csv", row.names = FALSE)

cat("\n==== FS search: counters vs Table S13, consistency stage ====\n")
s0 <- rd("S0"); str(s0$filter_counts); str(s0$n_candidates)
cat("n_hr_subgroups (find.grps$out.found$hr.subgroups rows):", s0$n_hr_subgroups, "\n")
cat("closed-form evaluations:", s0$n_closed_calls, " NA closed forms:", s0$n_closed_na, " fallbacks:", s0$n_fallback, "\n")
cat("p_cons >= 0.90:", sum(s0$pcons_closed$rate_closed >= 0.90), "\n")
cat("out_sg$result (candidates passing p* 0.90):\n"); print(s0$out_sg_result)
fit_S0 <- readRDS("fits/fit_S0.rds")
cat("\nS0 grp.consistency names:\n"); print(names(fit_S0$grp.consistency))
cat("S0 args_call_all$consistency_method:", fit_S0$args_call_all$consistency_method, "\n")
cat("S0 consistency_algorithm:", fit_S0$consistency_algorithm, "\n")
cat("S1 MR family size n_family:", fit_S1$mr_inference$n_family, "\n")
cat("\n==== D1 / G1 admission and family facts ====\n")
for (r in c("D1", "G1")) { x <- rd(r); cat(r, ": label", x$label, "| n_family", x$mr$n_family,
  "| closed-form calls in consistency stage", x$n_closed_calls, "| consistency_algorithm", format(x$consistency_algorithm), "\n") }
cat("\n==== per-rung MR seed / draws / field draws actually used ====\n")
for (r in c("S1", "S1r", "S2", "S3", "D1", "G1", "D3", "G3")) { m <- rd(r)$mr
  cat(sprintf("%-4s seed %s draws %s field R_out %s R_in %s | n_out_used H %s, n_in_used_mean H %.2f, n_out_used Hc %s | mr_sec %.2f field_sec H %.2f Hc %.2f\n",
      r, format(m$mr_seed), m$mr_draws, m$field_R_out, m$field_R_in, m$field_n_out_used_H, m$field_n_in_used_mean_H,
      m$field_n_out_used_Hc, m$mr_sec, m$field_sec_H, m$field_sec_Hc)) }
