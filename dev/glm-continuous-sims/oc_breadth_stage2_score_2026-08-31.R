# =============================================================================
# oc_breadth_stage2_score_2026-08-31.R -- §4 scoring of the measured MD120 cell
# against the locked forecast (commit c9cb0ca2).  Every predicted number is
# READ from oc_breadth_ladder_2026-08-30_{forecast120,gate}.rds and
# oc_wrapper_grid_corrected_2026-08-30.rds; every measured number is COMPUTED
# from the payloads.  Task: dev/tasks/cc_task_oc_breadth_stage2_2026-08-31.md
# Machinery: sgdef_selection_2026-08-29.R (rule evaluation on df_super,
# signatures, TVD), oc_wrapper_confs_compare_2026-08-30.R (between-rule gap,
# multinomial floor), oc_wrapper_grid_corrected_2026-08-30.R (measured columns).
# Writes dev/glm-continuous-sims/oc_breadth_stage2_2026-08-31_score.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
DIR <- file.path("dev", "glm-continuous-sims")
PD  <- "quarto/simulations/actg175/continuous/mr_md_harm"
say <- function(...) cat(sprintf(...), "\n", sep = "")
f6 <- function(x) sprintf("%.6f", x)
read_pay <- function(p) { z <- readRDS(p); nm <- names(z); nm[is.na(nm)] <- ""; setNames(z, nm) }
set.seed(20260831L)   # the multinomial noise-floor simulation only

FC <- readRDS(file.path(DIR, "oc_breadth_ladder_2026-08-30_forecast120.rds"))
G0 <- readRDS(file.path(DIR, "oc_breadth_ladder_2026-08-30_gate.rds"))
CORR <- readRDS(file.path(DIR, "oc_wrapper_grid_corrected_2026-08-30.rds"))
GT <- readRDS(file.path(DIR, "oc_breadth_stage2_2026-08-31_gate.rds"))
c1_star <- FC$c1_star; k_inter <- FC$k_inter
stopifnot(identical(c1_star, GT$c1_star), identical(k_inter, GT$k_inter))
P1 <- read_pay(file.path(PD, "fs_maxeffCons_mr_md120_knoise0_n500_s1000_d5000", "fs_maxeffCons_mr_md120_knoise0_n500_res_1_1000.rds"))
P2 <- read_pay(file.path(PD, "fs_maxeffCons_mr_md120_knoise0_n500_c1star_s1000_d5000", "fs_maxeffCons_mr_md120_knoise0_n500_c1star_res_1_1000.rds"))
P0 <- read_pay(file.path(PD, "fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000", "fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds"))
P40 <- read_pay(CORR$alt$payloads[["500"]])
say("[0] payloads: run1 n_sims %d thr %s pkg %s | run2 n_sims %d thr %s pkg %s | null n_sims %d thr %s | MD40 n_sims %d",
    P1$meta$n_sims, format(P1$meta$effect_threshold, digits = 15), P1$meta$pkg_version,
    P2$meta$n_sims, format(P2$meta$effect_threshold, digits = 15), P2$meta$pkg_version, P0$meta$n_sims, P0$meta$effect_threshold, P40$meta$n_sims)
stopifnot(P1$meta$effect_threshold == 30, identical(P2$meta$effect_threshold, c1_star),
          P1$meta$seed_base == 8316951L, P2$meta$seed_base == 8316951L, P0$meta$seed_base == 8316951L,
          abs(P1$truth$effect_Q + 120) < 1e-9, abs(P2$truth$effect_Q + 120) < 1e-9, abs(P1$truth$beta_inter - k_inter) < 1e-12)
# the driver-recorded meta in full (same knobs as MD40 except the cell)
for (k in c("n_sample", "mr_draws", "sg_focus", "consistency_method", "consistency_threshold", "seed_base", "parallel_mode", "n_workers"))
  stopifnot(identical(P1$meta[[k]], P40$meta[[k]]), identical(P2$meta[[k]], P40$meta[[k]]))

# ---- DGMs and frames (Stage 1 / sgdef machinery) --------------------------------
driver_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L)); actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat_orig <- ifelse(actg_df$arms == 1L, 1L, 0L); actg_df$treat <- 1L - actg_df$treat_orig
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40; actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L)); actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L)); actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L)); actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo); actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  for (v in c("hemo","homo","drugs","race","gender","symptom","str2")) actg_df[[v]] <- as.factor(actg_df[[v]])
  actg_df
}
build <- function(data, model, k_inter) generate_glm_dgm(
  data = data, factor_vars = paste0("z", 1:12), outcome_var = "cd4_change", treatment_var = "treat",
  outcome_type = "continuous", effect_measure = "MD", subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = model, k_treat = 1, k_inter = k_inter, adverse_outcome = FALSE, n_super = 5000L, seed = 8316951L, verbose = FALSE)
dfd <- driver_frame()
dgm120 <- build(dfd, "alt", k_inter); dgm_null <- build(dfd, "null", 0)
frame120 <- fs_build_eval_frame(dgm120, outcome_type = "continuous"); frame_null <- fs_build_eval_frame(dgm_null, outcome_type = "continuous")
stopifnot(identical(dgm120$df_super$mu1, P1$truth$effect_Q * 0 + dgm120$df_super$mu1))
sc <- fs_dgm_scale(dgm120)$regions; tauQc <- abs(sc$m_tau[sc$region == "Qc"]); bint <- abs(sc$m_tau[sc$region == "Q"]) - tauQc
stopifnot(abs(abs(sc$m_tau[sc$region == "Q"]) - 120) < 1e-9)
fs_args <- CORR$fs_args
fam <- fs_oc_family_enumerate(dgm120, fs_args, n = 500, max_M = 5000L, verbose = FALSE)
stopifnot(identical(fam$lab, G0$family$lab), identical(fam$beta_g, FC$beta_g), fam$M == 1696L)
say("[0] DGM rebuilt: tauQc %.6f bint %.6f ; family M = %d identical to the forecast rung's", tauQc, bint, fam$M)

# signatures (sgdef, verbatim)
clause_sig <- function(cl) {
  cl <- trimws(cl); neg <- grepl("^!", cl)
  cl <- gsub("^!\\s*", "", cl); cl <- gsub("^\\{|\\}$", "", cl); cl <- trimws(cl)
  if (grepl("<=", cl, fixed = TRUE)) { v <- trimws(strsplit(cl, "<=", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) ">" else "<=")) }
  if (grepl(">", cl, fixed = TRUE) && !grepl(">=", cl, fixed = TRUE)) { v <- trimws(strsplit(cl, ">", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "<=" else ">")) }
  if (grepl("!=", cl, fixed = TRUE)) { v <- trimws(strsplit(cl, "!=", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "=1" else "=0")) }
  if (grepl("==", cl, fixed = TRUE)) { v <- trimws(strsplit(cl, "==", fixed = TRUE)[[1]][1]); return(paste0(v, if (neg) "=0" else "=1")) }
  paste0(cl, if (neg) "=0" else "=1")
}
rule_sig <- function(rule) { cls <- strsplit(rule, " & ", fixed = TRUE)[[1]]; paste(sort(vapply(cls, clause_sig, character(1))), collapse = " & ") }
tvd <- function(a, b) 0.5 * sum(abs(a - b))
align <- function(x, keys) { z <- setNames(rep(0, length(keys)), keys); z[names(x)] <- x; z }
se_mean <- function(x) { x <- x[is.finite(x)]; stats::sd(x) / sqrt(length(x)) }
bse <- function(p, n) sqrt(p * (1 - p) / n)

# population functionals of a realized rule on a frame (sgdef §3, verbatim)
eval_rules <- function(rules, frame, is_null) {
  inQ <- frame$flag_harm == 1; PQ <- mean(inQ)
  do.call(rbind, lapply(rules, function(r) {
    m <- forestsearch:::.fs_resolve_membership(frame, r)
    if (!identical(m$status, "ok") && !identical(m$status, "empty"))
      return(data.frame(sg_def = r, status = m$status, Pg_pop = NA_real_, nPg_pop = NA_real_, PQg_pop = NA_real_, sens_pop = NA_real_, spec_pop = NA_real_, beta_pop = NA_real_, stringsAsFactors = FALSE))
    g <- m$in_region; Pg <- mean(g); PgQ <- mean(g & inQ); PgQc <- mean(g & !inQ)
    PQg <- if (Pg > 0) PgQ / Pg else NA_real_
    data.frame(sg_def = r, status = m$status, Pg_pop = Pg, nPg_pop = 500 * Pg, PQg_pop = if (is_null) 0 else PQg,
               sens_pop = if (is_null) NA_real_ else PgQ / PQ, spec_pop = 1 - PgQc / (1 - PQ),
               beta_pop = if (is_null) tauQc else tauQc + bint * PQg, stringsAsFactors = FALSE)
  }))
}
prep <- function(P, frame, is_null) {
  r <- P$results; stopifnot(all(r$status %in% c("DETECTED", "NO-DETECTION")))
  r$declared_30 <- r$status == "DETECTED"
  r$That <- r$nv_H_est
  r$beta_or <- P$oc$targets$orient * r$betaHhat_H
  d <- r[r$declared_30, ]
  ev <- eval_rules(unique(d$sg_def), frame, is_null)
  stopifnot(all(ev$status == "ok"))
  r <- merge(r, ev[, c("sg_def", "Pg_pop", "nPg_pop", "PQg_pop", "sens_pop", "spec_pop", "beta_pop")], by = "sg_def", all.x = TRUE, sort = FALSE)
  r <- r[order(r$sim_id), ]
  if (!is_null) stopifnot(max(abs(r$beta_or[r$declared_30] - r$beta_pop[r$declared_30])) < 1e-8)   # the payload's exact beta(Hhat) == tauQc + bint*PQg
  r$sig <- NA_character_; r$sig[r$declared_30] <- vapply(r$sg_def[r$declared_30], rule_sig, character(1), USE.NAMES = FALSE)
  r
}
R1 <- prep(P1, frame120, FALSE); R2 <- prep(P2, frame120, FALSE); R0 <- prep(P0, frame_null, TRUE)
say("[0] run1: declared_30 %d / %d ; run2 (c1*): declared %d / %d ; null: declared_30 %d / %d ; T-hat range run1 [%.2f, %.2f]",
    sum(R1$declared_30), nrow(R1), sum(R2$declared_30), nrow(R2), sum(R0$declared_30), nrow(R0), min(R1$That, na.rm = TRUE), max(R1$That, na.rm = TRUE))
R1$declared_star <- R1$declared_30 & R1$That >= c1_star

# ---- 4.1 power at c1* -------------------------------------------------------------
n1 <- nrow(R1); pw <- mean(R1$declared_star); pw_se <- bse(pw, n1)
fc_pw <- FC$star$resample$det_rate; fc_pw_se <- FC$star$resample$det_rate_se
say("[4.1] power at c1* = %s: %d / %d = %.4f (SE %.4f) vs forecast %.4f (MC SE %.4f) ; gap %+.4f = %.2f combined SE",
    f6(c1_star), sum(R1$declared_star), n1, pw, pw_se, fc_pw, fc_pw_se, pw - fc_pw, (pw - fc_pw) / sqrt(pw_se^2 + fc_pw_se^2))
# strict vs non-strict (the search drops hr <= c1; the post-hoc rule is >=): ties
say("[4.1] replicates with T-hat exactly == c1*: %d", sum(R1$declared_30 & R1$That == c1_star))

# ---- 4.2 the measured declaration curve ---------------------------------------------
C1 <- 0:300
grid_alt <- FC$grid$resample$table; stopifnot(identical(as.numeric(grid_alt$c1), as.numeric(C1)))
curve <- data.frame(c1 = C1,
                    measured = sapply(C1, function(c) mean(R1$declared_30 & R1$That >= c)),
                    predicted = grid_alt$det_rate, predicted_se = grid_alt$det_rate_se)
curve$measured_se <- bse(curve$measured, n1); curve$gap <- curve$measured - curve$predicted
curve$z <- curve$gap / sqrt(curve$measured_se^2 + curve$predicted_se^2)
i_max <- which.max(abs(curve$gap))
say("[4.2] alternative curve: largest |gap| = %.4f at c1 = %d (measured %.4f, predicted %.4f, %.2f SE) ; |gap| > 0.02 at c1 in [%s] ; max |z| = %.2f at c1 = %d",
    abs(curve$gap[i_max]), curve$c1[i_max], curve$measured[i_max], curve$predicted[i_max], curve$z[i_max],
    paste(range(curve$c1[abs(curve$gap) > 0.02]), collapse = ", "), max(abs(curve$z)), curve$c1[which.max(abs(curve$z))])
print(curve[curve$c1 %in% c(0, 30, 60, 80, 100, 110, 120, 125, 130, 133, 135, 136, 140, 145, 150, 160, 170, 180, 200, 220, 250, 300), ], digits = 4, row.names = FALSE)

# ---- 4.3 the direct run as the gate on post-hoc scoring ----------------------------------
stopifnot(identical(R1$sim_id, R2$sim_id))
d2 <- mean(R2$declared_30); d2_se <- bse(d2, nrow(R2))
agree_decl <- R2$declared_30 == R1$declared_star
both <- R2$declared_30 & R1$declared_star
same_def <- rep(NA, n1); same_def[both] <- R1$sg_def[both] == R2$sg_def[both]
same_nharm <- rep(NA, n1); same_nharm[both] <- R1$n_harm[both] == R2$n_harm[both]
same_That <- rep(NA, n1); same_That[both] <- abs(R1$That[both] - R2$That[both]) < 1e-9
# membership: reconstruct the replicate's sample and resolve both rules on it wherever the strings differ
memb_same <- same_def
chk <- which(both & !same_def)
if (length(chk)) for (i in chk) {
  sd_i <- 8316951L + R1$sim_id[i]
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i); df <- simulate_from_glm_dgm(dgm120, n = 500, seed = sd_i)
  ma <- forestsearch:::.fs_resolve_membership(df, R1$sg_def[i]); mb <- forestsearch:::.fs_resolve_membership(df, R2$sg_def[i])
  memb_same[i] <- identical(ma$in_region, mb$in_region)
}
say("[4.3] run 2 (effect.threshold = c1*): declaration rate %d / %d = %.4f (SE %.4f) ; run1 post-hoc %.4f ; per replicate: declaration agrees %d / %d (disagree %d) ; among %d joint declarations: sg_def identical %d, membership identical %d, n_harm identical %d, T-hat identical %d",
    sum(R2$declared_30), nrow(R2), d2, d2_se, pw, sum(agree_decl), n1, sum(!agree_decl), sum(both), sum(same_def, na.rm = TRUE), sum(memb_same, na.rm = TRUE), sum(same_nharm, na.rm = TRUE), sum(same_That, na.rm = TRUE))
dis <- R1$sim_id[!agree_decl | (both & !memb_same)]
if (length(dis)) {
  say("[4.3] DISAGREEMENTS (sim_id | run1 declared_30 / T-hat / winner | run2 declared / T-hat / winner):")
  for (i in match(dis, R1$sim_id)) say("   %4d | %s / %.4f / %s | %s / %s / %s", R1$sim_id[i], R1$declared_30[i], R1$That[i], R1$sg_def[i], R2$declared_30[i], if (is.na(R2$That[i])) "NA" else sprintf("%.4f", R2$That[i]), R2$sg_def[i])
}
measured_cell <- if (length(dis)) "run2" else "run1_posthoc"
say("[4.3] post-hoc scoring gate: %s ; the measured cell for §4.5-4.6 is %s", if (length(dis)) "FAILED (disagreements above)" else "PASS (complete agreement)", measured_cell)
D <- if (measured_cell == "run2") R2[R2$declared_30, ] else R1[R1$declared_star, ]

# ---- 4.4 type-I at c1* ---------------------------------------------------------------
R0$declared_star <- R0$declared_30 & R0$That >= c1_star
n0 <- nrow(R0); t1 <- mean(R0$declared_star); t1_se <- bse(t1, n0)
fc_t1 <- FC$null_at_star$det_rate; fc_t1_se <- FC$null_at_star$det_rate_se
say("[4.4] null type-I at c1*: %d / %d = %.4f (SE %.4f) vs forecast %.5f (MC SE %.5f) ; gap %+.4f = %.2f combined SE",
    sum(R0$declared_star), n0, t1, t1_se, fc_t1, fc_t1_se, t1 - fc_t1, (t1 - fc_t1) / sqrt(t1_se^2 + fc_t1_se^2))
null_sw <- CORR$null$sweep$table[CORR$null$sweep$table$consistency_method == "resample", ]
null_iv <- CORR$null$invert$table[CORR$null$invert$table$consistency_method == "resample", ]
null_curve <- data.frame(c1 = C1, measured = sapply(C1, function(c) mean(R0$declared_30 & R0$That >= c)))
null_curve$measured_se <- bse(null_curve$measured, n0)
null_curve$predicted <- NA_real_; null_curve$predicted_se <- NA_real_
m <- match(null_sw$c1, null_curve$c1); null_curve$predicted[m] <- null_sw$det_rate; null_curve$predicted_se[m] <- null_sw$det_rate_se
null_curve$gap <- null_curve$measured - null_curve$predicted
# the inversion points (c1 at targets 0.05..0.95) and c1*
null_pts <- data.frame(c1 = c(null_iv$value, c1_star), predicted = c(null_iv$achieved, fc_t1), predicted_se = c(null_iv$achieved_se, fc_t1_se))
null_pts$measured <- sapply(null_pts$c1, function(c) mean(R0$declared_30 & R0$That >= c)); null_pts$measured_se <- bse(null_pts$measured, n0)
null_pts$gap <- null_pts$measured - null_pts$predicted; null_pts$z <- null_pts$gap / sqrt(null_pts$measured_se^2 + null_pts$predicted_se^2)
say("[4.4] null curve vs the corrected null sweep (c1 20..120 by 5): largest |gap| %.4f at c1 = %d ; at the null inversion points and c1*:",
    max(abs(null_curve$gap), na.rm = TRUE), null_curve$c1[which.max(abs(null_curve$gap))])
print(null_pts, digits = 4, row.names = FALSE)
print(null_curve[!is.na(null_curve$predicted), ], digits = 4, row.names = FALSE)

# ---- 4.5 the declared rule, conditional on declaration at c1* ------------------------------
row_meas <- function(D) {
  data.frame(
    quantity = c("EnH_sample", "EnH_pop", "Eppv_sample", "Eppv_pop", "Esens_sample", "Esens_pop", "Espec_sample", "Espec_pop", "EbetaH_pop", "Enaive_bias"),
    measured = c(mean(D$n_harm), mean(D$nPg_pop), mean(D$ppv), mean(D$PQg_pop), mean(D$sens), mean(D$sens_pop), mean(D$spec), mean(D$spec_pop), mean(D$beta_or), mean(D$That - D$beta_or)),
    se = c(se_mean(D$n_harm), se_mean(D$nPg_pop), se_mean(D$ppv), se_mean(D$PQg_pop), se_mean(D$sens), se_mean(D$sens_pop), se_mean(D$spec), se_mean(D$spec_pop), se_mean(D$beta_or), se_mean(D$That - D$beta_or)),
    n = nrow(D))
}
fc_row <- function(p) c(EnH_sample = p$EnH, EnH_pop = p$EnH, Eppv_sample = p$Eppv, Eppv_pop = p$Eppv, Esens_sample = p$Esens, Esens_pop = p$Esens,
                        Espec_sample = p$Espec, Espec_pop = p$Espec, EbetaH_pop = p$EbetaH, Enaive_bias = p$Enaive_bias)
t45 <- row_meas(D); t45$forecast_resample <- fc_row(FC$star$resample)[t45$quantity]; t45$forecast_split <- fc_row(FC$star$split)[t45$quantity]
t45$gap <- t45$measured - t45$forecast_resample; t45$z <- t45$gap / t45$se
say("[4.5] the declared rule given declaration at c1* (%s, n = %d):", measured_cell, nrow(D)); print(t45, digits = 5, row.names = FALSE)
# the between-rule gap at this harm (confs-compare §3 definition): measured-frequency-weighted population size minus the analytic sel_c-weighted expectation
pop_realized <- mean(D$nPg_pop)
between <- data.frame(gate = c("resample", "split"), analytic_EnH = c(FC$star$resample$EnH, FC$star$split$EnH), pop_of_realized = pop_realized,
                      between = pop_realized - c(FC$star$resample$EnH, FC$star$split$EnH), within = mean(D$n_harm) - pop_realized, measured_EnH = mean(D$n_harm),
                      md40_between = c(2.111, 2.116))
# the sel_c-reweighted check (analytic weights on the measured signatures' population sizes)
sig_an <- vapply(fam$lab, rule_sig, character(1), USE.NAMES = FALSE)
size_sig <- tapply(D$nPg_pop, D$sig, mean)
rw <- function(selc) { A <- tapply(selc, sig_an, sum); w <- A[names(size_sig)]; w[is.na(w)] <- 0; sum(w * size_sig) / sum(w) }
between$reweighted <- c(rw(FC$star$resample$sel_c), rw(FC$star$split$sel_c))
say("[4.5] between-rule size gap at this harm, beside MD40's +2.11 (resample) / +2.12 (split):"); print(between, digits = 5, row.names = FALSE)
# PPV / sensitivity population-versus-sample offsets, for the record
say("[4.5] population-vs-sample offsets (sample - population, paired): size %+.3f (SE %.3f) ; PPV %+.4f (%.4f) ; sens %+.4f (%.4f) ; spec %+.4f (%.4f)",
    mean(D$n_harm - D$nPg_pop), se_mean(D$n_harm - D$nPg_pop), mean(D$ppv - D$PQg_pop), se_mean(D$ppv - D$PQg_pop),
    mean(D$sens - D$sens_pop), se_mean(D$sens - D$sens_pop), mean(D$spec - D$spec_pop), se_mean(D$spec - D$spec_pop))

# ---- 4.6 composition ----------------------------------------------------------------------
nD <- nrow(D)
m95 <- mean(D$PQg_pop >= 0.95); m95_se <- bse(m95, nD)
say("[4.6] measured selection mass on rules with PQg >= 0.95: %d / %d = %.4f (SE %.4f) vs forecast %.4f (resample) / %.4f (split)",
    sum(D$PQg_pop >= 0.95), nD, m95, m95_se, FC$sel_star$resample$mass_PQg95, FC$sel_star$split$mass_PQg95)
say("[4.6] measured mass on rules with PQg == 1: %.4f ; Q itself (age > 34 & preanti <= 744.5) selected: %d times", mean(D$PQg_pop == 1), sum(D$PQg_pop == 1 & abs(D$sens_pop - 1) < 1e-12))
# closest analytic labels: by Jaccard on df_super membership
memb_fam <- fam$memb
closest <- function(rule) {
  g <- forestsearch:::.fs_resolve_membership(frame120, rule)$in_region
  j <- colSums(memb_fam & g) / colSums(memb_fam | g)
  k <- which.max(j); c(lab = fam$lab[k], jaccard = sprintf("%.3f", j[k]), same_sig = identical(sig_an[k], rule_sig(rule)))
}
tv <- sort(table(D$sg_def), decreasing = TRUE)
top_real <- data.frame(sg_def = names(tv)[1:15], n = as.integer(tv[1:15]), freq = as.numeric(tv[1:15]) / nD)
ev_top <- eval_rules(top_real$sg_def, frame120, FALSE)
top_real$Pg <- ev_top$Pg_pop[match(top_real$sg_def, ev_top$sg_def)]; top_real$PQg <- ev_top$PQg_pop[match(top_real$sg_def, ev_top$sg_def)]; top_real$beta_g <- ev_top$beta_pop[match(top_real$sg_def, ev_top$sg_def)]
cl <- t(sapply(top_real$sg_def, closest)); top_real$closest_analytic <- cl[, "lab"]; top_real$jaccard <- cl[, "jaccard"]
o <- order(FC$star$resample$sel_c, decreasing = TRUE)[1:15]
top_an <- data.frame(lab = fam$lab[o], sel_c = FC$star$resample$sel_c[o], Pg = fam$Pg[o], PQg = fam$PQg[o], beta_g = fam$beta_g[o])
# measured frequency of each analytic top-15 (by signature)
Mf_sig <- table(D$sig) / nD
top_an$measured_sig_freq <- as.numeric(Mf_sig[sig_an[o]]); top_an$measured_sig_freq[is.na(top_an$measured_sig_freq)] <- 0
say("[4.6] top 15 realized rules by frequency (verbatim), with population Pg/PQg/beta and the closest analytic label:"); print(top_real, digits = 4, row.names = FALSE)
say("[4.6] top 15 analytic candidates by sel_c (resample) with the measured frequency of their signature:"); print(top_an, digits = 4, row.names = FALSE)
# signature-level distribution, TVD, multinomial floor (confs-compare, verbatim)
A_r <- tapply(FC$star$resample$sel_c, sig_an, sum); A_s <- tapply(FC$star$split$sel_c, sig_an, sum)
keys <- union(names(A_r), names(Mf_sig)); A <- align(A_r, keys); As <- align(A_s, keys); Mf <- align(Mf_sig, keys)
tvd_r <- tvd(A, Mf); tvd_s <- tvd(As, Mf)
floor_sim <- function(A, B = 2000L) { p <- A[A > 0]; p <- p / sum(p); v <- replicate(B, { x <- stats::rmultinom(1, nD, p)[, 1] / nD; tvd(x, p) }); c(mean = mean(v), q025 = unname(stats::quantile(v, 0.025)), q975 = unname(stats::quantile(v, 0.975))) }
fl <- floor_sim(A)
say("[4.6] TVD over signatures: resample %.4f, split %.4f ; multinomial floor at %d declarations: mean %.4f (2.5-97.5%%: %.4f-%.4f) ; excess %.4f ; support analytic %d / measured %d ; analytic mass on never-selected signatures %.4f ; measured mass absent from the family %.4f",
    tvd_r, tvd_s, nD, fl[["mean"]], fl[["q025"]], fl[["q975"]], tvd_r - fl[["mean"]], sum(A > 0), sum(Mf > 0), sum(A[Mf == 0]), sum(Mf[A == 0]))
top_sig <- names(sort(A, decreasing = TRUE))[1:15]
sig_tab <- data.frame(signature = top_sig, analytic_resample = round(A[top_sig], 4), analytic_split = round(As[top_sig], 4), measured = round(Mf[top_sig], 4), measured_n = as.integer(round(Mf[top_sig] * nD)))
say("[4.6] top 15 analytic signatures beside measured frequency:"); print(sig_tab, row.names = FALSE)
top_msig <- names(sort(Mf, decreasing = TRUE))[1:15]
say("[4.6] top 15 measured signatures beside analytic:"); print(data.frame(signature = top_msig, measured = round(Mf[top_msig], 4), measured_n = as.integer(round(Mf[top_msig] * nD)), analytic_resample = round(A[top_msig], 4)), row.names = FALSE)

# ---- 4.7 at the driver's c1 = 30 ----------------------------------------------------------
measured_from_payload <- function(pl) {   # oc_wrapper_grid_corrected_2026-08-30.R, verbatim
  oc <- pl$oc; res <- pl$results
  idc <- oc$identification[oc$identification$convention == "conditional", ]
  est <- oc$estimation; orH <- est[est$block == "H" & est$estimator == "oracle", ]; nvH <- est[est$block == "H" & est$estimator == "naive",  ]
  det_rows <- res[res$detected %in% TRUE, ]
  list(det_rate = idc$detection, EnH = idc$mean_size_H, Esens = idc$sens, Espec = idc$spec, Eppv = idc$ppv, Enpv = idc$npv,
       EbetaH = orH$avg - orH$bias_beta, EbetaH_results = mean(oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE), Enaive_bias = nvH$bias_beta,
       det_rate_se = sqrt(idc$detection * (1 - idc$detection) / pl$meta$n_sims),
       EnH_se = se_mean(det_rows$n_harm), Esens_se = se_mean(det_rows$sens), Espec_se = se_mean(det_rows$spec), Eppv_se = se_mean(det_rows$ppv), Enpv_se = se_mean(det_rows$npv),
       EbetaH_se = se_mean(det_rows$betaHhat_H), Enaive_bias_se = se_mean(det_rows$nv_H_est - oc$targets$orient * det_rows$betaHhat_H))
}
m120 <- measured_from_payload(P1); m40 <- CORR$alt$measured[["500"]]; m40_re <- measured_from_payload(P40)
stopifnot(isTRUE(all.equal(m40$det_rate, m40_re$det_rate)), isTRUE(all.equal(m40$EnH, m40_re$EnH)), isTRUE(all.equal(m40$Enaive_bias, m40_re$Enaive_bias)))
q7 <- c("det_rate", "EnH", "Eppv", "Esens", "Espec", "EbetaH", "Enaive_bias")
t47 <- data.frame(quantity = q7,
                  md120_measured = sapply(q7, function(q) m120[[q]]), md120_se = sapply(q7, function(q) m120[[paste0(q, "_se")]]),
                  md120_forecast_resample = sapply(q7, function(q) FC$named$resample$driver[[q]]), md120_forecast_split = sapply(q7, function(q) FC$named$split$driver[[q]]),
                  md40_measured = sapply(q7, function(q) m40[[q]]), md40_se = sapply(q7, function(q) m40[[paste0(q, "_se")]]),
                  md40_analytic_resample = sapply(q7, function(q) CORR$alt$runs$n500_resample$pred[[q]]))
t47$gap120 <- t47$md120_measured - t47$md120_forecast_resample; t47$z120 <- t47$gap120 / t47$md120_se
say("[4.7] at c1 = 30 (saturated): measured MD120 beside the forecast rung's row and MD40's measured / analytic columns (corrected report):"); print(t47, digits = 5, row.names = FALSE)
# the same population-based rows at c1 = 30, for symmetry with 4.5
D30 <- R1[R1$declared_30, ]; t47b <- row_meas(D30); t47b$forecast_resample <- fc_row(FC$named$resample$driver)[t47b$quantity]; t47b$gap <- t47b$measured - t47b$forecast_resample; t47b$z <- t47b$gap / t47b$se
between30 <- data.frame(analytic_EnH = FC$named$resample$driver$EnH, pop_of_realized = mean(D30$nPg_pop), between = mean(D30$nPg_pop) - FC$named$resample$driver$EnH, within = mean(D30$n_harm) - mean(D30$nPg_pop))
say("[4.7] population rows at c1 = 30 (n = %d):", nrow(D30)); print(t47b, digits = 5, row.names = FALSE); print(between30, digits = 5, row.names = FALSE)

# ---- §5 verdicts ------------------------------------------------------------------------
verdict <- function(gap, se_comb, kind = "noise") if (abs(gap) <= 2 * se_comb) "within noise" else if (kind == "pvs") "population-versus-sample" else sprintf("a gap (%+.4g)", gap)
V <- data.frame(item = c("power at c1*", "type-I at c1*", "EbetaH at c1*", "specificity at c1* (pop)", "specificity at c1* (sample)", "PPV at c1* (pop)", "PPV at c1* (sample)",
                         "sensitivity at c1* (pop)", "sensitivity at c1* (sample)", "E|H| at c1* (pop)", "E|H| at c1* (sample)", "naive bias at c1*", "between-rule gap vs +2.11"),
                measured = c(pw, t1, t45$measured[t45$quantity == "EbetaH_pop"], t45$measured[t45$quantity == "Espec_pop"], t45$measured[t45$quantity == "Espec_sample"],
                             t45$measured[t45$quantity == "Eppv_pop"], t45$measured[t45$quantity == "Eppv_sample"], t45$measured[t45$quantity == "Esens_pop"], t45$measured[t45$quantity == "Esens_sample"],
                             t45$measured[t45$quantity == "EnH_pop"], t45$measured[t45$quantity == "EnH_sample"], t45$measured[t45$quantity == "Enaive_bias"], between$between[1]),
                predicted = c(fc_pw, fc_t1, FC$star$resample$EbetaH, FC$star$resample$Espec, FC$star$resample$Espec, FC$star$resample$Eppv, FC$star$resample$Eppv,
                              FC$star$resample$Esens, FC$star$resample$Esens, FC$star$resample$EnH, FC$star$resample$EnH, FC$star$resample$Enaive_bias, 2.111),
                se_comb = c(sqrt(pw_se^2 + fc_pw_se^2), sqrt(t1_se^2 + fc_t1_se^2), t45$se[t45$quantity == "EbetaH_pop"], t45$se[t45$quantity == "Espec_pop"], t45$se[t45$quantity == "Espec_sample"],
                            t45$se[t45$quantity == "Eppv_pop"], t45$se[t45$quantity == "Eppv_sample"], t45$se[t45$quantity == "Esens_pop"], t45$se[t45$quantity == "Esens_sample"],
                            t45$se[t45$quantity == "EnH_pop"], t45$se[t45$quantity == "EnH_sample"], t45$se[t45$quantity == "Enaive_bias"], t45$se[t45$quantity == "EnH_pop"]))
V$gap <- V$measured - V$predicted; V$z <- V$gap / V$se_comb
V$verdict <- mapply(function(g, s, it) verdict(g, s, if (grepl("sample", it) && grepl("PPV|sens", it)) "pvs" else "noise"), V$gap, V$se_comb, V$item)
say("[5] verdicts (within noise = |measured - predicted| <= 2 combined SE):"); print(V, digits = 5, row.names = FALSE)

saveRDS(list(c1_star = c1_star, k_inter = k_inter, forecast_commit = "c9cb0ca2", measured_cell = measured_cell,
             power = list(rate = pw, se = pw_se, n = n1, k = sum(R1$declared_star), forecast = fc_pw, forecast_se = fc_pw_se),
             curve = curve, run2 = list(rate = d2, se = d2_se, agree_decl = sum(agree_decl), both = sum(both), same_def = sum(same_def, na.rm = TRUE), memb_same = sum(memb_same, na.rm = TRUE), disagreements = dis),
             typeI = list(rate = t1, se = t1_se, n = n0, k = sum(R0$declared_star), forecast = fc_t1, forecast_se = fc_t1_se), null_curve = null_curve, null_pts = null_pts,
             rule = t45, between = between, composition = list(m95 = m95, m95_se = m95_se, top_real = top_real, top_an = top_an, sig_tab = sig_tab, tvd = c(resample = tvd_r, split = tvd_s), floor = fl, A = A, As = As, Mf = Mf),
             c1_30 = list(table = t47, pop = t47b, between = between30), verdicts = V,
             reps = list(run1 = R1[, c("sim_id", "status", "declared_30", "declared_star", "That", "sg_def", "sig", "n_harm", "ppv", "sens", "spec", "beta_or", "Pg_pop", "nPg_pop", "PQg_pop", "sens_pop", "spec_pop", "beta_pop")],
                         run2 = R2[, c("sim_id", "status", "declared_30", "That", "sg_def", "n_harm")], null = R0[, c("sim_id", "status", "declared_30", "declared_star", "That", "sg_def", "n_harm")]),
             pkg_version = as.character(packageVersion("forestsearch")), built_at = Sys.time()),
        file.path(DIR, "oc_breadth_stage2_2026-08-31_score.rds"))
say("[score] written")
