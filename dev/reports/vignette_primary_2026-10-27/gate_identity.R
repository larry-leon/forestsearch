# Gate §2.4 / §2.5 of dev/tasks/cc_task_vignette_primary_2026-10-27.md.
# Run from the package root in a fresh Rscript; uses the installed forestsearch.
# Re-runs the vignette's code (knitr::purl + source), then compares the three fits
# with stage 0's fits S3 / D3 / G3 and ladder_quantities.csv, and records the
# funnel, candidates, partition and declaration values.

stage0 <- "dev/reports/vignettes_stage0_2026-10-27"
cat("forestsearch", as.character(packageVersion("forestsearch")),
    "| Built:", packageDescription("forestsearch")$Built, "\n")
cat("vignette md5:", unname(tools::md5sum("vignettes/forestsearch.qmd")), "\n\n")

# ---- 1. re-run the vignette ----------------------------------------------------
purled <- tempfile(fileext = ".R")
knitr::purl("vignettes/forestsearch.qmd", output = purled, quiet = TRUE)
t0 <- proc.time()[["elapsed"]]
sink(tempfile())                                    # the vignette's own printed output
source(purled, local = globalenv(), echo = FALSE)
sink()
cat(sprintf("purled vignette sourced in %.1f s\n\n", proc.time()[["elapsed"]] - t0))

# ---- 2. identity with stage 0 --------------------------------------------------
quantities <- function(fit) {
  g <- fit$mr_inference; cc <- g$complement; fb <- g$field
  in_H <- fit$df.est$treat.recommend[match(df$id, fit$df.est$id)] == 0
  list(region = fit$sg.harm, n_selected = g$n_selected, n_Hc = cc$n,
       naive_H = g$naive$est, naive_lo_H = g$naive$lower, naive_hi_H = g$naive$upper,
       naive_Hc = cc$naive$est, naive_lo_Hc = cc$naive$lower, naive_hi_Hc = cc$naive$upper,
       corrected_H = g$debiased$est, corrected_Hc = cc$debiased$est,
       field_lower_H = fb$lower_1s, field_s_upper_Hc = fb$complement$upper_1s_s,
       bonf_lower_H = fb$joint_s$bonf_lower_H, bonf_upper_Hc = fb$joint_s$bonf_upper_Hc,
       n_family = g$n_family,
       E_H = sum(df$status[in_H]), E_Hc = sum(df$status[!in_H]))
}
csv_col <- c(region = "region", n_selected = "N_H", n_Hc = "N_Hc",
             naive_H = "naive_H", naive_lo_H = "naive_lo_H", naive_hi_H = "naive_hi_H",
             naive_Hc = "naive_Hc", naive_lo_Hc = "naive_lo_Hc", naive_hi_Hc = "naive_hi_Hc",
             corrected_H = "corrected_H", corrected_Hc = "corrected_Hc",
             field_lower_H = "field_lower_H", field_s_upper_Hc = "field_s_upper_Hc",
             bonf_lower_H = "bonf_lower_H", bonf_upper_Hc = "bonf_upper_Hc",
             n_family = "n_family", E_H = "E_H", E_Hc = "E_Hc")
csv <- read.csv(file.path(stage0, "ladder_quantities.csv"), colClasses = "character")
fmt <- function(x) {
  if (is.character(x)) paste(x, collapse = " & ")
  else if (x == round(x)) sprintf("%d", as.integer(x)) else sprintf("%.6f", x)
}
fits <- list(FS = list(fit_fs, "S3"), DINA = list(fit_dina, "D3"), GRF = list(fit_grf, "G3"))
rows <- list()
for (nm in names(fits)) {
  v  <- quantities(fits[[nm]][[1]])
  rung <- fits[[nm]][[2]]
  s0 <- quantities(readRDS(file.path(stage0, "fits", paste0("fit_", rung, ".rds"))))
  cr <- csv[csv$rung == rung, ]
  for (q in names(v)) {
    csv_val <- cr[[csv_col[[q]]]]
    vig_fmt <- fmt(v[[q]])
    rows[[length(rows) + 1]] <- data.frame(
      quantity = q, fit = nm, rung = rung,
      vignette = vig_fmt, stage0_rds = fmt(s0[[q]]), csv = csv_val,
      equal_rds = identical(v[[q]], s0[[q]]),
      equal_csv = if (is.numeric(v[[q]])) identical(vig_fmt, sprintf(if (v[[q]] == round(v[[q]])) "%d" else "%.6f",
                                                                     as.numeric(csv_val)))
                  else identical(vig_fmt, csv_val))
  }
}
cmp <- do.call(rbind, rows)
old <- options(width = 200)
print(cmp, row.names = FALSE)
options(old)
cat(sprintf("\nGATE 2.4: %d comparisons vs rds, %d equal; vs csv, %d equal -> %s\n\n",
            nrow(cmp), sum(cmp$equal_rds), sum(cmp$equal_csv),
            if (all(cmp$equal_rds) && all(cmp$equal_csv)) "PASS" else "FAIL"))

# ---- 2.4 recorded checks (not gated) -------------------------------------------
cat("funnel:", funnel$candidates, "(expect 2211 1410 121 9 1744)\n")
cat(sprintf("candidates: %d rows; band >= %.2f; selected N %d; max N %d; max-beta N %d (expect 9, 2.03, 75, 82, 61)\n",
            nrow(cand), band, sel$N, max(cand$N), cand$N[which.max(cand$hr)]))
cat("three-way partition (expect N/events 35/19, 40/22, 78/42, 7/4, 526/212):\n")
print(partition[order(-partition$N), ], row.names = FALSE)

# ---- 2.5 declaration -----------------------------------------------------------
k75 <- tb$kappa_hat[tb$c0 == 0.75]
decl <- data.frame(
  quantity = c("kappa_hat (c0 0.75)", "implied level (c0 0.75)", "fw_size (c0 0.75)",
               "T_H", "T_max", "kappa_hat (c2)", "fw_size (c2)"),
  value    = sprintf("%.3f", c(k75, 2 * pnorm(k75) - 1, tb$fw_size[tb$c0 == 0.75],
                               T_H, T_max, cal_c2$kappa_hat, cal_c2$fw_size)),
  table_S17 = c("2.726", "0.994", "0.667", "2.438", "2.564", "3.503", "0.974"))
decl$equal <- decl$value == decl$table_S17
cat("\n")
print(decl, row.names = FALSE)
admits <- c(lengths(cal$c0$admitted_calibrated), length(cal_c2$admitted_calibrated))
cat("admits counts (c0 0.70, 0.75, 0.80, 0.85, c2):", admits, "\n")
cat("full declaration table:\n")
print(declaration, row.names = FALSE)
cat(sprintf("\nGATE 2.5: %s\n", if (all(decl$equal) && all(admits == 0)) "PASS" else "FAIL"))
