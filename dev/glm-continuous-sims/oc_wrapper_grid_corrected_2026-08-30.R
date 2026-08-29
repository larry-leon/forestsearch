# =============================================================================
# oc_wrapper_grid_corrected_2026-08-30.R -- the analytic family re-enumerated on
# the DRIVER'S confounder list (13 variables, str2 included)
# -----------------------------------------------------------------------------
# Task: dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md, §3.
#
# The 08-29 scripts (oc_wrapper_grid_2026-08-29.R, oc_wrapper_grid_sweep_2026-
# 08-29.R, oc_wrapper_null_2026-08-29.R) enumerated the analytic family on a
# 12-variable list that omitted `str2`; the three drivers passed 13 (§2 of the
# task report).  This script is those three scripts with ONE input changed --
# `fs_args$confounders.name` is the driver's `confs` vector, verbatim -- and the
# same seed (20260825) and draw count (2e5), so the candidate space is the only
# difference.  The 08-29 scripts and their .rds are left untouched as the record.
#
# Per cell it produces, exactly as the 08-29 scripts did:
#   alt n in {500, 700} x gate in {resample, split}:
#     fs_oc_predict() with MC SEs; M and stage counts; floor; settings; seed;
#     wall-clock; the c1 sweep (fs_oc_grid, block = Inf) at the driver's c2;
#     the inversions (fs_oc_invert, vector targets -- the order-statistic path)
#   null n = 500 x both gates:
#     false declaration and MC SE, per-candidate P1 range, L_eff, EnH, EbetaH,
#     naive bias, mass_below, M; the c1 sweep; the inversions
#   measured columns recomputed from payload$oc / payload$results (never typed)
#
# Usage (from the repository root; parts are independent and were run as three
# concurrent processes, then merged):
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.R alt500
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.R alt700
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.R null500
#   OC_PART_DIR=<dir> Rscript dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.R merge
# Output of merge: dev/glm-continuous-sims/oc_wrapper_grid_corrected_2026-08-30.rds
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))

part <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(part %in% c("alt500", "alt700", "null500", "merge"))
PART_DIR <- Sys.getenv("OC_PART_DIR", unset = tempdir())
part_path <- function(p) file.path(PART_DIR, sprintf("oc_corrected_part_%s.rds", p))

PAY_DIR <- file.path("quarto", "simulations", "actg175", "continuous", "mr_md_harm")
pay_alt <- function(n) file.path(
  PAY_DIR, sprintf("fs_maxeffCons_mr_md40_knoise0_n%d_s1000_d5000", n),
  sprintf("fs_maxeffCons_mr_md40_knoise0_n%d_res_1_1000.rds", n))
PAY_NULL <- file.path(PAY_DIR, "fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000",
                      "fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds")
OUT  <- file.path("dev", "glm-continuous-sims", "oc_wrapper_grid_corrected_2026-08-30.rds")
G_OLD  <- readRDS(file.path("dev", "glm-continuous-sims", "oc_wrapper_grid_2026-08-29.rds"))
N0_OLD <- readRDS(file.path("dev", "glm-continuous-sims", "oc_wrapper_null_2026-08-29.rds"))

SEED  <- G_OLD$seed;  DRAWS <- G_OLD$draws          # 20260825, 2e5 -- unchanged
stopifnot(SEED == 20260825L, DRAWS == 2e5)
SWEEP_C1     <- G_OLD$sweep$c1                     # seq(20, 120, by = 5)
TARGETS_ALT  <- G_OLD$invert$targets               # 0.80 0.90 0.95
TARGETS_NULL <- N0_OLD$invert$targets              # 0.05 ... 0.95

# ---- the ONE corrected input: the driver's confounders.name ------------------
# sim_fs_maxeffCons_mr_md40_knoise0_n{500,700}_batch_1_1000.qmd:
#   analysis_continuous_vars <- c("age","preanti","wtkg","karnof","cd40","cd80")   (L242)
#   analysis_binary_vars     <- c("hemo","homo","drugs","race","gender","symptom") (L243)
#   include_str2 <- TRUE; analysis_binary_vars <- c(analysis_binary_vars, "str2") (L267-269)
#   confounders_analysis <- c(analysis_continuous_vars, analysis_binary_vars)     (L302)
#   confs <- c(confounders_analysis, <no noise: k_random_noise = 0>)              (L550-551)
#   forestsearch(confounders.name = confs, ...)                                   (L562)
fs_args <- G_OLD$fs_args
fs_args$confounders.name <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                              "hemo", "homo", "drugs", "race", "gender", "symptom",
                              "str2")
stopifnot(length(fs_args$confounders.name) == 13L,
          identical(fs_args$confounders.name[1:12], G_OLD$fs_args$confounders.name),
          identical(fs_args[-1], G_OLD$fs_args[-1]))   # every other argument unchanged

# ---- the ACTG175 frame, as every 08-29 script builds it ----------------------
actg_frame <- function() {
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- 1L - ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$z1  <- as.factor(ifelse(actg_df$age > 34, 1L, 0L))
  actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= 744.5, 1L, 0L))
  actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
  actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
  actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
  actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
  actg_df$z7  <- as.factor(actg_df$hemo);  actg_df$z8  <- as.factor(actg_df$homo)
  actg_df$z9  <- as.factor(actg_df$drugs); actg_df$z10 <- as.factor(actg_df$race)
  actg_df$z11 <- as.factor(actg_df$gender);actg_df$z12 <- as.factor(actg_df$symptom)
  actg_df
}
build_dgm <- function(model, k_inter) generate_glm_dgm(
  data = actg_frame(), factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = model, k_treat = 1, k_inter = k_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)

read_pay <- function(p) { z <- readRDS(p); nm <- names(z); nm[is.na(nm)] <- ""; setNames(z, nm) }

# ---- measured columns from a payload (oc_wrapper_grid_2026-08-29.R, verbatim) -
measured_from_payload <- function(pl) {
  oc <- pl$oc; res <- pl$results
  idc <- oc$identification[oc$identification$convention == "conditional", ]
  est <- oc$estimation
  orH <- est[est$block == "H" & est$estimator == "oracle", ]
  nvH <- est[est$block == "H" & est$estimator == "naive",  ]
  det_rows <- res[res$detected %in% TRUE, ]
  list(
    det_rate    = idc$detection,
    EnH         = idc$mean_size_H,
    Esens       = idc$sens,
    Espec       = idc$spec,
    Eppv        = idc$ppv,
    Enpv        = idc$npv,
    EbetaH      = orH$avg - orH$bias_beta,
    EbetaH_results = mean(oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE),
    Enaive_bias = nvH$bias_beta,
    n_used      = idc$n_used,
    n_sims      = pl$meta$n_sims,
    det_rate_se = sqrt(idc$detection * (1 - idc$detection) / pl$meta$n_sims),
    EnH_se      = stats::sd(det_rows$n_harm, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$n_harm))),
    Esens_se    = stats::sd(det_rows$sens, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$sens))),
    Espec_se    = stats::sd(det_rows$spec, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$spec))),
    Eppv_se     = stats::sd(det_rows$ppv,  na.rm = TRUE) / sqrt(sum(!is.na(det_rows$ppv))),
    Enpv_se     = stats::sd(det_rows$npv,  na.rm = TRUE) / sqrt(sum(!is.na(det_rows$npv))),
    EbetaH_se   = stats::sd(det_rows$betaHhat_H, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$betaHhat_H))),
    Enaive_bias_se = stats::sd(det_rows$nv_H_est - oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE) /
                     sqrt(sum(!is.na(det_rows$nv_H_est))),
    source = "payload$oc$identification[conditional], payload$oc$estimation[H: oracle, naive], payload$results",
    M = NA_integer_)
}
# null measured column (oc_wrapper_null_2026-08-29.R, verbatim)
measured_null_from_payload <- function(pl) {
  oc <- pl[["oc"]]; res <- pl[["results"]]; meta <- pl[["meta"]]
  idc <- oc$identification[oc$identification$convention == "conditional", ]
  est <- oc$estimation
  nvH <- est[est$block == "H" & est$estimator == "naive", ]
  det_rows <- res[res$detected %in% TRUE, ]
  list(
    det_rate = idc$detection, det_rate_se = sqrt(idc$detection * (1 - idc$detection) / meta$n_sims),
    EnH = idc$mean_size_H, EnH_se = stats::sd(det_rows$n_harm, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$n_harm))),
    Esens = idc$sens,
    Espec = idc$spec, Espec_se = stats::sd(det_rows$spec, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$spec))),
    Eppv = idc$ppv, Enpv = idc$npv,
    EbetaH = mean(oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE),
    EbetaH_se = stats::sd(det_rows$betaHhat_H, na.rm = TRUE) / sqrt(sum(!is.na(det_rows$betaHhat_H))),
    Enaive_bias = nvH$bias_beta,
    Enaive_bias_se = stats::sd(det_rows$nv_H_est - oc$targets$orient * det_rows$betaHhat_H, na.rm = TRUE) /
                     sqrt(sum(!is.na(det_rows$nv_H_est))),
    n_used = idc$n_used, n_sims = meta$n_sims,
    source = "null payload $oc$identification[conditional], $oc$estimation[H, naive], $results")
}

family_record <- function(fam, n, secs) list(
  M = fam$M, counts = fam$counts, n.min = fam$args_used$n.min,
  floor_Pg = fam$args_used$n.min / n, Pg_range = range(fam$Pg),
  PQ = fam$PQ, cuts = fam$cuts, lab = fam$lab, Pg = fam$Pg, PQg = fam$PQg,
  sens_g = fam$sens_g, spec_g = fam$spec_g,
  beta_g = fam$beta_g, se_g = fam$se_g, null = isTRUE(fam$null),
  args_used = fam$args_used, secs = secs)

# =============================================================================
# ALT cell (n = 500 or 700): predict x 2 gates, sweep, inversions
# =============================================================================
run_alt <- function(n) {
  pl <- readRDS(pay_alt(n)); truth <- pl$truth
  dgm <- build_dgm("alt", truth$beta_inter)
  gate <- c(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q),
            abs(dgm$subgroup_info$proportion - truth$prevalence_Q),
            abs(dgm$model_params$beta_inter - truth$beta_inter))
  stopifnot(all(gate < 1e-9))
  sc <- fs_dgm_scale(dgm)
  stopifnot(isTRUE(all.equal(sc$regions, pl$scale$regions)),
            isTRUE(all.equal(sc$regions, G_OLD$scale$regions)))
  pm <- pl$meta
  stopifnot(pm$effect_threshold == fs_args$effect.threshold,
            pm$consistency_threshold == fs_args$consistency.threshold,
            pm$sg_focus == fs_args$sg_focus, pm$n_sample == n,
            pm$consistency_method == "resample")
  stopifnot(all(fs_args$confounders.name %in% names(dgm$df_super)))
  cat(sprintf("[alt%d] fixture rebuilt and gated; scale identical to payload and to the 08-29 grid\n", n))
  measured <- measured_from_payload(pl)
  # the measured column must be what the 08-29 grid stored (same payload)
  stopifnot(isTRUE(all.equal(measured, G_OLD$measured[[as.character(n)]])))

  # -- family on 13 variables
  t0 <- proc.time()[3]
  fam <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = TRUE)
  t_fam <- proc.time()[3] - t0
  cat(sprintf("[alt%d] corrected family: M = %d (08-29: %d); cut columns %d (08-29: %d); %.1f s\n",
              n, fam$M, G_OLD$families[[as.character(n)]]$M,
              fam$counts[["cut_columns"]], G_OLD$families[[as.character(n)]]$counts[["cut_columns"]], t_fam))
  stopifnot("str2" %in% fam$cuts, any(grepl("str2", fam$lab)))
  family <- family_record(fam, n, t_fam)

  # -- the two single-point evaluations at the driver's (c1, c2)
  runs <- list()
  for (g in c("resample", "split")) {
    gc(); t0 <- proc.time()[3]
    pred <- fs_oc_predict(family = fam, n = n,
                          c1 = fs_args$effect.threshold, c2 = fs_args$consistency.threshold,
                          consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                          draws = DRAWS, seed = SEED)
    secs <- proc.time()[3] - t0
    peak <- sum(gc()[, 6])
    pred$family <- NULL
    runs[[sprintf("n%d_%s", n, g)]] <- list(
      n = n, gate = g, pred = pred, secs = secs, peak_MB = peak,
      settings = pred$settings, seed = SEED, draws = DRAWS)
    cat(sprintf("[alt%d] gate = %-8s M = %d det_rate = %.4f EnH = %.2f Esens = %.4f Eppv = %.4f naive = %.2f  %.0f s peak %.0f MB\n",
                n, g, fam$M, pred$det_rate, pred$EnH, pred$Esens, pred$Eppv, pred$Enaive_bias, secs, peak))
  }

  # -- the c1 sweep at the driver's c2, both gates, block = Inf
  t0 <- proc.time()[["elapsed"]]
  sw <- fs_oc_grid(family = fam, n = n, c1 = SWEEP_C1, c2 = fs_args$consistency.threshold,
                   consistency_method = c("resample", "split"),
                   pconsistency = fs_args$pconsistency.threshold,
                   draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE)
  cat(sprintf("[alt%d] sweep: %d rows in %.0f s\n", n, nrow(sw$table), proc.time()[["elapsed"]] - t0))
  for (g in c("resample", "split")) {
    i <- which(sw$table$n == n & sw$table$consistency_method == g & sw$table$c1 == fs_args$effect.threshold)
    new <- sw$results[[i]]; new$family <- NULL
    ok <- identical(new, runs[[sprintf("n%d_%s", n, g)]]$pred)
    cat(sprintf("[alt%d] %-8s sweep point at c1 = %g identical() to fs_oc_predict(): %s\n", n, g, fs_args$effect.threshold, ok))
    stopifnot(ok)
  }
  print(sw$timing)

  # -- the inversions, vector targets (order-statistic reduction: one draw set per gate)
  inv <- list(); inv_rows <- list()
  for (g in c("resample", "split")) {
    t0 <- proc.time()[["elapsed"]]
    iv <- fs_oc_invert(family = fam, n = n, target = TARGETS_ALT, solve_for = "c1",
                       c2 = fs_args$consistency.threshold, consistency_method = g,
                       pconsistency = fs_args$pconsistency.threshold, draws = DRAWS, seed = SEED)
    secs <- proc.time()[["elapsed"]] - t0
    tb <- attr(iv, "table"); tb$n <- n; tb$consistency_method <- g; tb$secs <- secs
    inv[[g]] <- iv; inv_rows[[g]] <- tb
    cat(sprintf("[alt%d] invert %-8s: %d targets in %.0f s\n", n, g, length(TARGETS_ALT), secs))
  }
  inv_tab <- do.call(rbind, inv_rows); rownames(inv_tab) <- NULL
  print(inv_tab, digits = 5)

  out <- list(n = n, family = family, runs = runs, measured = measured,
              sweep = list(table = sw$table, c1 = SWEEP_C1, c2 = fs_args$consistency.threshold,
                           block = Inf, timing = sw$timing),
              invert = list(table = inv_tab, targets = TARGETS_ALT, objects = inv),
              scale = sc, payload = pay_alt(n), truth = truth)
  saveRDS(out, part_path(sprintf("alt%d", n)))
  cat("[alt", n, "] part written: ", part_path(sprintf("alt%d", n)), "\n", sep = "")
}

# =============================================================================
# NULL cell (n = 500): predict x 2 gates, sweep, inversions
# =============================================================================
run_null <- function() {
  n <- 500
  pl <- read_pay(PAY_NULL)
  truth <- pl[["truth"]]; meta <- pl[["meta"]]
  stopifnot(is.na(truth$effect_Q), truth$prevalence_Q == 0, truth$beta_inter == 0,
            meta$consistency_method == "resample", meta$n_sample == n)
  dgm <- build_dgm("null", 0)
  gate <- c(Q_empty      = sum(dgm$df_super$flag_harm) == 0L,
            effect_Q_NA  = is.na(dgm$hazard_ratios$harm_subgroup),
            effect_Qc    = abs(dgm$hazard_ratios$no_harm_subgroup - truth$effect_Qc) < 1e-9,
            beta_inter   = isTRUE(all.equal(dgm$model_params$beta_inter, truth$beta_inter)),
            effect_ITT   = abs(dgm$hazard_ratios$overall - truth$effect_ITT) < 1e-9,
            model_null   = identical(dgm$model, "null"))
  cat("[null500] rebuild gate:\n"); print(gate); stopifnot(all(gate))
  stopifnot(meta$effect_threshold == fs_args$effect.threshold,
            meta$consistency_threshold == fs_args$consistency.threshold)
  c1 <- fs_args$effect.threshold; c2 <- fs_args$consistency.threshold
  measured <- measured_null_from_payload(pl)
  stopifnot(isTRUE(all.equal(measured, N0_OLD$measured)))

  t0 <- proc.time()[["elapsed"]]
  fam <- fs_oc_family_enumerate(dgm, fs_args, n = n, max_M = 5000L, verbose = TRUE)
  t_enum <- proc.time()[["elapsed"]] - t0
  stopifnot(isTRUE(fam$null), all(fam$PQg == 0), all(is.na(fam$sens_g)))
  stopifnot("str2" %in% fam$cuts, any(grepl("str2", fam$lab)))
  cat(sprintf("[null500] corrected null family: M = %d (08-29: %d); common effect %.6f; %.1f s\n",
              fam$M, N0_OLD$family$M, fam$beta_g[1], t_enum))
  family <- family_record(fam, n, t_enum)
  family$beta_common <- fam$beta_g[1]; family$scale <- fam$scale

  runs <- list()
  for (g in c("resample", "split")) {
    gc(); t0 <- proc.time()[["elapsed"]]
    p <- fs_oc_predict(family = fam, n = n, c1 = c1, c2 = c2,
                       consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                       draws = DRAWS, seed = SEED)
    secs <- proc.time()[["elapsed"]] - t0
    peak <- sum(gc()[, 6])
    p1 <- p$P1
    L_eff <- log(1 - p$det_rate) / log(1 - max(p1))
    p$family <- NULL
    runs[[g]] <- list(pred = p, secs = secs, peak_MB = peak, p1_range = range(p1), L_eff = L_eff)
    cat(sprintf("[null500] %-8s false declaration %.4f (SE %.4f); P1 range %.4f..%.4f; L_eff %.2f; EnH %.2f; EbetaH %.3f; naive %.2f; mass_below %.3f; %.0f s\n",
                g, p$det_rate, p$det_rate_se, min(p1), max(p1), L_eff, p$EnH, p$EbetaH, p$Enaive_bias, p$mass_below, secs))
  }

  t0 <- proc.time()[["elapsed"]]
  sw <- fs_oc_grid(family = fam, n = n, c1 = SWEEP_C1, c2 = c2,
                   consistency_method = c("resample", "split"),
                   pconsistency = fs_args$pconsistency.threshold,
                   draws = DRAWS, block = Inf, seed = SEED, verbose = TRUE)
  cat(sprintf("[null500] sweep: %d rows in %.0f s\n", nrow(sw$table), proc.time()[["elapsed"]] - t0))
  for (g in c("resample", "split")) {
    i <- which(sw$table$consistency_method == g & sw$table$c1 == c1)
    new <- sw$results[[i]]; new$family <- NULL
    ok <- identical(new, runs[[g]]$pred)
    cat(sprintf("[null500] %-8s sweep point at c1 = %g identical() to fs_oc_predict(): %s\n", g, c1, ok))
    stopifnot(ok)
  }
  print(sw$timing)

  inv <- list(); inv_rows <- list()
  for (g in c("resample", "split")) {
    t0 <- proc.time()[["elapsed"]]
    iv <- fs_oc_invert(family = fam, n = n, target = TARGETS_NULL, solve_for = "c1", c2 = c2,
                       consistency_method = g, pconsistency = fs_args$pconsistency.threshold,
                       draws = DRAWS, seed = SEED)
    secs <- proc.time()[["elapsed"]] - t0
    tb <- attr(iv, "table"); tb$consistency_method <- g; tb$secs <- secs
    inv[[g]] <- iv; inv_rows[[g]] <- tb
    cat(sprintf("[null500] invert %-8s: %d targets in %.0f s\n", g, length(TARGETS_NULL), secs))
  }
  inv_tab <- do.call(rbind, inv_rows); rownames(inv_tab) <- NULL
  print(inv_tab, digits = 5)

  cmp <- merge(sw$table[sw$table$consistency_method == "resample", c("c1", "det_rate", "det_rate_se")],
               sw$table[sw$table$consistency_method == "split", c("c1", "det_rate")],
               by = "c1", suffixes = c("_resample", "_split"))
  cmp$diff <- cmp$det_rate_split - cmp$det_rate_resample

  out <- list(n = n, family = family, runs = runs, measured = measured,
              sweep = list(table = sw$table, c1 = SWEEP_C1, c2 = c2, block = Inf, timing = sw$timing),
              invert = list(table = inv_tab, targets = TARGETS_NULL, objects = inv),
              gate_compare = cmp, truth = truth, meta = meta, payload = PAY_NULL)
  saveRDS(out, part_path("null500"))
  cat("[null500] part written: ", part_path("null500"), "\n", sep = "")
}

# =============================================================================
# MERGE: one .rds, shaped like the 08-29 grid + null .rds side by side
# =============================================================================
run_merge <- function() {
  a500 <- readRDS(part_path("alt500")); a700 <- readRDS(part_path("alt700")); nul <- readRDS(part_path("null500"))
  stopifnot(a500$n == 500, a700$n == 700, nul$n == 500)
  alt <- list(
    runs = c(a500$runs, a700$runs),
    families = list(`500` = a500$family, `700` = a700$family),
    measured = list(`500` = a500$measured, `700` = a700$measured),
    sweep = list(table = rbind(a500$sweep$table, a700$sweep$table), c1 = SWEEP_C1,
                 c2 = fs_args$consistency.threshold, block = Inf, seed = SEED, draws = DRAWS),
    sweep_timing = rbind(a500$sweep$timing, a700$sweep$timing),
    invert = list(table = rbind(a500$invert$table, a700$invert$table), targets = TARGETS_ALT,
                  objects = list(`500` = a500$invert$objects, `700` = a700$invert$objects)),
    scale = a500$scale, payloads = c(`500` = a500$payload, `700` = a700$payload))
  null <- list(
    runs = nul$runs, family = nul$family, measured = nul$measured,
    sweep = nul$sweep, invert = nul$invert, gate_compare = nul$gate_compare,
    truth = nul$truth, meta = nul$meta, n = nul$n, payload = nul$payload)
  out <- list(
    alt = alt, null = null,
    fs_args = fs_args, fs_args_superseded = G_OLD$fs_args,
    seed = SEED, draws = DRAWS,
    superseded = c(alt = "dev/glm-continuous-sims/oc_wrapper_grid_2026-08-29.rds",
                   null = "dev/glm-continuous-sims/oc_wrapper_null_2026-08-29.rds"),
    task = "dev/tasks/cc_task_oc_wrapper_confs_2026-08-30.md",
    built_at = Sys.time(),
    pkg_version = as.character(utils::packageVersion("forestsearch")))
  saveRDS(out, OUT)
  cat("written:", OUT, "\n")
  cat(sprintf("M corrected: alt500 %d (08-29 %d) alt700 %d (08-29 %d) null500 %d (08-29 %d)\n",
              alt$families[["500"]]$M, G_OLD$families[["500"]]$M,
              alt$families[["700"]]$M, G_OLD$families[["700"]]$M,
              null$family$M, N0_OLD$family$M))
}

switch(part,
       alt500  = run_alt(500),
       alt700  = run_alt(700),
       null500 = run_null(),
       merge   = run_merge())
