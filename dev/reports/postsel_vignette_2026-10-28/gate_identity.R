# Identity gate — post-selection vignette vs stage 0 (2026-10-28).
# Run from the package root:  Rscript dev/reports/postsel_vignette_2026-10-28/gate_identity.R
# Purls vignettes/survival-post-selection.qmd, sources it in this fresh process
# (reruns the three fits and the three bootstraps), then compares with stage 0's
# saved M fits. Writes identity.csv, table_s2_check.csv, bootstrap.csv, cost.csv,
# recovery.csv into the record directory and exits non-zero on a gated failure.

rec <- "dev/reports/postsel_vignette_2026-10-28"
s0  <- "dev/reports/postsel_stage0_2026-10-28/fits"
purled <- file.path(rec, "survival-post-selection_purled.R")
knitr::purl("vignettes/survival-post-selection.qmd", output = purled, quiet = TRUE)

t_all <- proc.time()[["elapsed"]]
source(purled, local = globalenv(), echo = FALSE)
cat(sprintf("Purled vignette sourced in %.1f s\n", proc.time()[["elapsed"]] - t_all))

M <- list(FS   = readRDS(file.path(s0, "fit_M-FS.rds")),
          DINA = readRDS(file.path(s0, "fit_M-DINA.rds")),
          GRF  = readRDS(file.path(s0, "fit_M-GRF.rds")))
V <- list(FS = fit_fs, DINA = fit_dina, GRF = fit_grf)

getters <- list(
  sg.harm            = function(f) f$sg.harm,
  n_selected         = function(f) f$mr_inference$n_selected,
  complement_n       = function(f) f$mr_inference$complement$n,
  naive_est          = function(f) f$mr_inference$naive$est,
  naive_lower        = function(f) f$mr_inference$naive$lower,
  naive_upper        = function(f) f$mr_inference$naive$upper,
  c_naive_est        = function(f) f$mr_inference$complement$naive$est,
  c_naive_lower      = function(f) f$mr_inference$complement$naive$lower,
  c_naive_upper      = function(f) f$mr_inference$complement$naive$upper,
  debiased_est       = function(f) f$mr_inference$debiased$est,
  debiased_lower     = function(f) f$mr_inference$debiased$lower,
  debiased_upper     = function(f) f$mr_inference$debiased$upper,
  c_debiased_est     = function(f) f$mr_inference$complement$debiased$est,
  c_debiased_lower   = function(f) f$mr_inference$complement$debiased$lower,
  c_debiased_upper   = function(f) f$mr_inference$complement$debiased$upper,
  n_family           = function(f) f$mr_inference$n_family,
  selection_bias     = function(f) f$mr_inference$selection_bias,
  ci_method          = function(f) f$mr_inference$ci_method)
fmt <- function(x) if (is.null(x)) "NULL" else paste(format(x, digits = 15), collapse = " & ")

rows <- list()
for (m in names(V)) for (q in names(getters)) {
  a <- getters[[q]](V[[m]]); b <- getters[[q]](M[[m]])
  eq <- identical(a, b) ||
    (is.numeric(a) && is.numeric(b) && length(a) == length(b) && all(abs(a - b) == 0))
  rows[[length(rows) + 1]] <- data.frame(quantity = q, identifier = m, vignette = fmt(a),
                                         stage0 = fmt(b), equal = eq,
                                         gated = q != "ci_method")
}
ident <- do.call(rbind, rows)
write.csv(ident, file.path(rec, "identity.csv"), row.names = FALSE)
print(ident, right = FALSE)

# Field quantities, FS only — recorded, not gated
fv <- fit_fs$mr_inference$field; fm <- M$FS$mr_inference$field
fq <- list(lower_1s = function(x) x$lower_1s,
           complement_upper_1s_s = function(x) x$complement$upper_1s_s,
           joint_s_bonf_lower_H = function(x) x$joint_s$bonf_lower_H,
           joint_s_bonf_upper_Hc = function(x) x$joint_s$bonf_upper_Hc)
field_tab <- do.call(rbind, lapply(names(fq), function(q) {
  a <- fq[[q]](fv); b <- fq[[q]](fm)
  data.frame(quantity = q, vignette = fmt(a), stage0 = fmt(b), identical = identical(a, b))
}))
write.csv(field_tab, file.path(rec, "field_fs.csv"), row.names = FALSE)
print(field_tab, right = FALSE)

rc <- fv$recovery
recov <- data.frame(sens_H = rc$sens_H, ppv_H = rc$ppv_H, q10 = rc$q10, q50 = rc$q50,
                    q90 = rc$q90, share_equal_1 = rc$share_equal_1, n_used = rc$n_used)
write.csv(recov, file.path(rec, "recovery.csv"), row.names = FALSE)
print(recov)

# Table S2: naive and multiplier gated; bootstrap recorded
s2 <- cbind(side_by_side[, c("method", "region", "n")],
            naive_vignette = side_by_side$naive, naive_s2 = table_s2$naive,
            naive_eq = side_by_side$naive == table_s2$naive,
            mult_vignette = side_by_side$multiplier, mult_s2 = table_s2$multiplier,
            mult_eq = side_by_side$multiplier == table_s2$multiplier,
            boot_vignette = side_by_side$bootstrap, boot_s2 = table_s2$bootstrap,
            n_s2 = table_s2$n)
write.csv(s2, file.path(rec, "table_s2_check.csv"), row.names = FALSE)
print(s2, right = FALSE)

boot_tab <- do.call(rbind, lapply(names(boots), function(m) data.frame(
  identifier = m, B = B[[m]],
  H2 = boots[[m]]$H_estimates$H2, H2_lower = boots[[m]]$H_estimates$H2_lower,
  H2_upper = boots[[m]]$H_estimates$H2_upper,
  Hc_H2 = boots[[m]]$Hc_estimates$H2, Hc_H2_lower = boots[[m]]$Hc_estimates$H2_lower,
  Hc_H2_upper = boots[[m]]$Hc_estimates$H2_upper, boot_secs = boot_secs[[m]])))
write.csv(boot_tab, file.path(rec, "bootstrap.csv"), row.names = FALSE)
print(boot_tab)
write.csv(cost, file.path(rec, "cost.csv"), row.names = FALSE)
print(cost)
write.csv(data.frame(identifier = names(secs), fit_secs = unname(secs)),
          file.path(rec, "fit_secs.csv"), row.names = FALSE)

fail_ident <- ident[ident$gated & !ident$equal, ]
fail_s2 <- !(all(s2$naive_eq) && all(s2$mult_eq))
cat(sprintf("\nGATE identity: %s (%d gated mismatches)\n",
            if (nrow(fail_ident) == 0) "PASS" else "FAIL", nrow(fail_ident)))
cat(sprintf("GATE Table S2 naive/multiplier: %s\n", if (fail_s2) "FAIL" else "PASS"))
cat(sprintf("Total elapsed %.1f s\n", proc.time()[["elapsed"]] - t_all))
if (nrow(fail_ident) > 0 || fail_s2) quit(status = 1)
