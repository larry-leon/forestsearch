# =============================================================================
# oc_applied_eval_summary_2026-08-31.R
# Applied OC evaluation §4 post-processing: from the eleven stored rung
# checkpoints (grid tables + named-point predicts), build every table the
# application document reads -- type-I, the ladder, what is declared, the
# calibration curve, coherence checks -- plus the §5 se_g bands, into one
# object: oc_applied_eval_summary_2026-08-31.rds.  No draws are run here.
#
# MC SEs: det_rate carries the package's binomial SE.  The selection-weighted
# functionals E = sum(sel_c * x) are conditional means over the multinomial
# selection counts, so their delta-method MC SE is
#   sqrt( (sum(sel_c * x^2) - E^2) / (draws * det_rate) )
# (exact for EnH, Esens, Espec, Eppv, EbetaH, mass_below, whose x are fixed
# population quantities).  Enaive_bias averages a continuous per-draw noise
# whose within-winner moments are not stored; its SE is approximated as
#   sqrt( sum(sel_c * se_g^2) / (draws * det_rate) )
# treating the winner's noise scale as se_g (ignores the variance reduction
# from max-selection truncation, so mildly conservative) -- labeled approx.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))

out_dir <- "dev/glm-continuous-sims"
build <- readRDS(file.path(out_dir, "oc_applied_eval_build_2026-08-31.rds"))
res0  <- readRDS(file.path(out_dir, "stage0_oc_applied_2026-08-31.rds"))
q_rungs <- build$q_rungs
T_obs   <- build$context$T_obs
n       <- build$context$N

rungs <- lapply(seq_along(q_rungs), function(i)
  readRDS(file.path(out_dir,
                    sprintf("oc_applied_eval_rung%02d_2026-08-31.rds", i))))
stopifnot(!any(vapply(rungs, `[[`, logical(1), "missing_named")))
np <- rungs[[1L]]$named_points
c1_05 <- np$c1_05; c1_10 <- np$c1_10

# ---- the curves (long): rung x c1 -> det_rate --------------------------------
curves <- do.call(rbind, lapply(seq_along(rungs), function(i) {
  tb <- rungs[[i]]$table
  data.frame(rung = i, q = q_rungs[i], c1 = tb$c1,
             det_rate = tb$det_rate, det_rate_se = tb$det_rate_se)
}))

rate_at <- function(i, c1v) {
  tb <- rungs[[i]]$table
  j <- which(abs(tb$c1 - c1v) < 1e-9)
  stopifnot(length(j) == 1L)
  tb[j, c("det_rate", "det_rate_se")]
}
# det_rate is nonincreasing in c1, so "the c1 with rate >= target" read off
# the grid is the LARGEST integer c1 still holding the target (the breadth
# ladder's c1* convention; the smallest such c1 is trivially 0).
largest_int_c1_ge <- function(i, target) {
  tb <- rungs[[i]]$table
  ok <- tb$c1 %% 1 == 0 & tb$det_rate >= target
  if (!any(ok)) NA_real_ else max(tb$c1[ok])
}

# ---- §4.1 type-I -------------------------------------------------------------
tb1 <- rungs[[1L]]$table
type1 <- list(
  curve   = tb1[, c("c1", "det_rate", "det_rate_se")],
  at_analyst = rate_at(1L, 10),
  c1_05 = c1_05, rate_at_c1_05 = rate_at(1L, c1_05),
  c1_10 = c1_10, rate_at_c1_10 = rate_at(1L, c1_10))

# ---- §4.2 the ladder ---------------------------------------------------------
ladder <- do.call(rbind, lapply(seq_along(rungs), function(i) {
  r10  <- rate_at(i, 10)
  r05  <- rate_at(i, c1_05)
  r10t <- rate_at(i, c1_10)
  data.frame(
    rung = i, q = q_rungs[i],
    det_at_analyst = r10$det_rate,  det_at_analyst_se = r10$det_rate_se,
    power_at_c1_05 = r05$det_rate,  power_at_c1_05_se = r05$det_rate_se,
    power_at_c1_10 = r10t$det_rate, power_at_c1_10_se = r10t$det_rate_se,
    c1_pow80 = largest_int_c1_ge(i, 0.80),
    c1_pow90 = largest_int_c1_ge(i, 0.90),
    c1_pow95 = largest_int_c1_ge(i, 0.95))
}))

# ---- §4.3 what is declared at the analyst's c1 -------------------------------
fx_se <- function(p, x) {
  E <- sum(p$sel_c * x)
  sqrt(max(sum(p$sel_c * x^2) - E^2, 0) /
         (p$settings$draws * p$det_rate))
}
declared <- do.call(rbind, lapply(seq_along(rungs), function(i) {
  p <- rungs[[i]]$named$analyst
  fam <- p$family
  data.frame(
    rung = i, q = q_rungs[i],
    det_rate = p$det_rate, det_rate_se = p$det_rate_se,
    EnH = p$EnH,               EnH_se = fx_se(p, n * fam$Pg),
    Eppv = p$Eppv,             Eppv_se = fx_se(p, fam$PQg),
    Esens = p$Esens,           Esens_se = fx_se(p, fam$sens_g),
    Espec = p$Espec,           Espec_se = fx_se(p, fam$spec_g),
    EbetaH = p$EbetaH,         EbetaH_se = fx_se(p, fam$beta_g),
    Enaive_bias = p$Enaive_bias,
    Enaive_bias_se_approx = sqrt(sum(p$sel_c * fam$se_g^2) /
                                   (p$settings$draws * p$det_rate)),
    mass_below = p$mass_below,
    mass_below_se = fx_se(p, as.numeric(fam$beta_g < 10)))
}))

# ---- §4.4 the calibration curve ----------------------------------------------
calibration <- do.call(rbind, lapply(seq_along(rungs), function(i) {
  r <- rate_at(i, T_obs)
  data.frame(rung = i, q = q_rungs[i],
             p_ge_Tobs = r$det_rate, p_ge_Tobs_se = r$det_rate_se)
}))
# crossings by linear interpolation between adjacent rungs (stated as such)
cross_q <- function(level) {
  y <- calibration$p_ge_Tobs; x <- calibration$q
  hit <- which(diff(y >= level) != 0)
  if (!length(hit)) return(NA_real_)
  j <- hit[1L]
  x[j] + (level - y[j]) * (x[j + 1L] - x[j]) / (y[j + 1L] - y[j])
}
crossings <- data.frame(level = c(0.05, 0.5, 0.95),
                        q_cross = vapply(c(0.05, 0.5, 0.95), cross_q,
                                         numeric(1)))

# ---- §4.5 coherence checks ---------------------------------------------------
coherence <- list(
  power_monotone_analyst = all(diff(ladder$det_at_analyst) >= 0),
  power_monotone_c1_05   = all(diff(ladder$power_at_c1_05) >= 0),
  power_monotone_c1_10   = all(diff(ladder$power_at_c1_10) >= 0),
  zero_plus_EbetaH       = declared$EbetaH[1L],
  zero_plus_det          = declared$det_rate[1L])

# ---- timings -----------------------------------------------------------------
timings <- data.frame(
  rung = seq_along(rungs), q = q_rungs,
  t_grid_secs = vapply(rungs, `[[`, numeric(1), "t_grid"))

# ---- §5 se_g bands (if present) ----------------------------------------------
band_file <- file.path(out_dir, "oc_applied_eval_seg_band_2026-08-31.rds")
seg_band <- if (file.exists(band_file)) readRDS(band_file) else NULL

# ---- §5 triggered sensitivity (se_g -> se_direct at the top rung) ------------
# Beside-the-primary comparison at the top rung; a sensitivity, not adopted.
sens_file <- file.path(out_dir, "oc_applied_eval_sens_seg_2026-08-31.rds")
sensitivity <- NULL
if (file.exists(sens_file)) {
  sn <- readRDS(sens_file)
  stopifnot(!sn$missing_named)
  sn_rate <- function(c1v) {
    j <- which(abs(sn$table$c1 - c1v) < 1e-9)
    sn$table[j, c("det_rate", "det_rate_se")]
  }
  cmp <- function(src, r10, r05, r10t, rT, p) data.frame(
    run = src,
    det_at_analyst = r10$det_rate, det_at_analyst_se = r10$det_rate_se,
    power_at_c1_05 = r05$det_rate, power_at_c1_05_se = r05$det_rate_se,
    power_at_c1_10 = r10t$det_rate, power_at_c1_10_se = r10t$det_rate_se,
    p_ge_Tobs = rT$det_rate, p_ge_Tobs_se = rT$det_rate_se,
    EnH = p$EnH, EbetaH = p$EbetaH, Eppv = p$Eppv)
  sensitivity <- list(
    label = sn$label,
    table = rbind(
      cmp("primary (se_g)", rate_at(11L, 10), rate_at(11L, c1_05),
          rate_at(11L, c1_10), rate_at(11L, T_obs),
          rungs[[11L]]$named$analyst),
      cmp("sensitivity (se_direct)", sn_rate(10), sn_rate(c1_05),
          sn_rate(c1_10), sn_rate(T_obs), sn$named$analyst)),
    t_grid = sn$t_grid)
}

summary_obj <- list(
  meta = list(
    n = n, M = 4508L, q_rungs = q_rungs, T_obs = T_obs,
    beta_treat = build$context$beta_treat,
    analyst_c1 = 10, analyst_c2 = 5, pconsistency = 0.90,
    consistency_method = "resample", draws = 2e5, block = 5e4,
    seed = 8316951L, c1_05 = c1_05, c1_10 = c1_10,
    anchor = list(def = res0$anchor$def, n_H = res0$anchor$n_H,
                  T_obs = res0$anchor$T_obs, p_cons = res0$anchor$p_cons),
    purity = res0$purity, itt = res0$itt,
    stage0_null_det = 1.000,   # stage 0 §5: the null branch's oriented-abs
    scale = res0$scale,        # sigma and the q = 20 region table
    gates = build$gates,
    forestsearch_version = as.character(packageVersion("forestsearch")),
    built_at = Sys.time()),
  curves = curves,
  type1 = type1,
  ladder = ladder,
  declared = declared,
  calibration = calibration,
  crossings = crossings,
  coherence = coherence,
  timings = timings,
  seg_band = seg_band,
  sensitivity = sensitivity)

saveRDS(summary_obj,
        file.path(out_dir, "oc_applied_eval_summary_2026-08-31.rds"))
cat("[done] summary written.\n")
cat(sprintf("type-I at (10, 5): %.5f (SE %.5f)\n",
            type1$at_analyst$det_rate, type1$at_analyst$det_rate_se))
cat(sprintf("c1_05 = %s (rate %.5f); c1_10 = %s (rate %.5f)\n",
            format(c1_05), type1$rate_at_c1_05$det_rate,
            format(c1_10), type1$rate_at_c1_10$det_rate))
print(ladder, digits = 4)
print(calibration, digits = 4)
print(crossings, digits = 5)
