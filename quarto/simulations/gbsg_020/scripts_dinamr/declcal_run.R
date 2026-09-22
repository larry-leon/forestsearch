# declcal_run.R -- one cell of the declaration-calibration finite-sample
# evaluation (TASK_declcal_CAMPAIGN_2026-09-22_v2).  Called by declcal.sh.
#
# TRANSPLANT, NOT AUTHORED (task section 10).  Source: the committed template
# sim_fs_maxeffCons_fb_mr_field_m1_template.qmd, the per-replicate engine every
# nullid / nullc125 / nullmr cell rendered.  Carried over, with the template's
# line numbers at HEAD 66d1ad6c:
#   T:244-246    .env_chr / .env_int / .env_num
#   T:249, 294   n_sims, sim_id_start
#   T:318, 342   sg_focus, effect_neighborhood (FS_S7_FOCUS / FS_S7_NBHD)
#   T:357-397    target HR, n, z1_quantile, FS_S7_DGM and its guard
#   T:489-501    FS_S7_C1 / FS_S7_C2 and their guards
#   T:584-586    analysis_time, cens_adjust, n_super
#   T:603, 624-648  consistency_method, selection_rule, stop_threshold, c1/c2,
#                p*, fs_splits, maxk, n_min, d0/d1_min, use_* flags, forced
#                cuts, er grid
#   T:752-771    column names, confounders_base, seed_base, n_workers
#   T:798-983    build-dgm chunk, including the structural-null design-point
#                gate, verbatim except change (1) below
#   T:1224-1300  the replicate's data draw and the forestsearch() call
#   T:1671-1676  .safe_record's never-abort pattern
#   T:1689-1705  the "sims" foreach %dofuture% loop
# Named changes (and nothing else):
#   (1) DGM effect.  Under FS_S7_DGM=null with FS_S7_HR=1 the template's
#       uniroot on k_treat cannot run: gamma["treat"] <- k_treat * gamma["treat"]
#       (R/sim_aft_gbsg.R:417) reaches HR 1 only at k_treat = 0, and
#       .create_gbsg_dgm_() requires k_treat > 0 (:256).  The complete null
#       therefore sets k_treat = 1e-10 (log HR ~1e-11); the template's own
#       three-check design-point gate still runs and must pass.
#   (2) c1 / c2 = 1.0 / 1.0 via FS_S7_C1 / FS_S7_C2 (the template knobs).
#   (3) n via FS_S7_N; FS only (FS_S7_METHOD unset = "consistency").
#   (4) Capture formals: forestsearch() runs with mr_inference = FALSE; the
#       declaration field is captured on EVERY replicate by calling
#       fs_mr_inference(keep_declaration_field = TRUE, keep_field_matrix =
#       FALSE) on the replicate's enumerated pre-reduction family -- the same
#       enumeration forestsearch() hands MR (R/forestsearch_main.R:3713-3726),
#       over the cut matrix Z the search itself used (captured by trace() on
#       subgroup.search; no R/ edit), with the admission set the fit resolved
#       and the MR seed forestsearch() would use (seedit = seed_base + sim_id).
#       Verified identical to forestsearch()'s internal route (family, beta_hat,
#       sigma_D, Mstar) on a declaring replicate before the campaign.
#   (5) The per-replicate record schema of task section 5, plus a side table
#       (aux) holding the search's own declaration indicator (the fidelity
#       gate's other side), the condition message, and the wall split.
#   (6) Per-replicate time cap (setTimeLimit), chunked loop with the per-cell
#       1% abort/error gate, per-chunk checkpoint.
suppressPackageStartupMessages({
  library(forestsearch); library(survival); library(data.table)
  library(foreach); library(doFuture); library(future)
})

.env_chr <- function(k, d) { v <- Sys.getenv(k, unset = NA); if (is.na(v) || !nzchar(v)) d else v }
.env_int <- function(k, d) as.integer(.env_chr(k, d))
.env_num <- function(k, d) as.numeric(.env_chr(k, d))

# ---- knobs (template names) -------------------------------------------------
n_sims       <- .env_int("FS_S7_NSIMS", 2000L)
sim_id_start <- .env_int("FS_S7_START", 1L)
sg_focus     <- .env_chr("FS_S7_FOCUS", "effMaxSG")
stopifnot(identical(sg_focus, "effMaxSG"))
effect_neighborhood <- .env_num("FS_S7_NBHD", 0.20)
er_jcuts     <- .env_int("FS_S7_ER_JCUTS", 10L)
target_hr_harm <- .env_num("FS_S7_HR", 1.0)
n_sample     <- .env_int("FS_S7_N", 1000L)
harm_z1_quantile <- .env_num("FS_S7_Z1Q", 0.25)
dgm_model    <- .env_chr("FS_S7_DGM", "alt")
stopifnot(dgm_model %in% c("alt", "null"))
if (identical(dgm_model, "null") && nzchar(Sys.getenv("FS_S7_Z1Q")))
  stop("FS_S7_Z1Q is set but FS_S7_DGM = 'null' plants no region.", call. = FALSE)
thr_c1 <- .env_num("FS_S7_C1", 0.90)
thr_c2 <- .env_num("FS_S7_C2", 0.80)
stopifnot(is.finite(thr_c1), is.finite(thr_c2), thr_c1 > 0, thr_c2 > 0, thr_c2 <= thr_c1)
campaign_tag <- .env_chr("FS_S7_CAMPAIGN", "declcal")
stopifnot(grepl("^[A-Za-z0-9_]+$", campaign_tag))
subgroup_method <- .env_chr("FS_S7_METHOD", "consistency")
stopifnot(identical(subgroup_method, "consistency"))      # FS only (task section 4)

# declcal knobs
cell_id  <- .env_chr("DECLCAL_CELL", "X")
dgm_lab  <- .env_chr("DECLCAL_DGM_LABEL", dgm_model)
B_cal    <- .env_int("DECLCAL_BCAL", 2000L)
run_kind <- .env_chr("DECLCAL_MODE", "stage2")            # "pilot" | "stage2"
stopifnot(run_kind %in% c("pilot", "stage2"))
cap_s    <- .env_num("DECLCAL_CAP_S", 0)                  # per-replicate hard cap; 0 = none
chunk_n  <- .env_int("DECLCAL_CHUNK", 500L)
out_path <- .env_chr("DECLCAL_OUT", sprintf("results/%s_%s_res_%d_%d.rds",
                     campaign_tag, cell_id, sim_id_start, sim_id_start + n_sims - 1L))
pilot_B  <- c(500L, 1000L, 2000L)
if (identical(run_kind, "pilot")) stopifnot(B_cal == max(pilot_B))

analysis_time <- 84; cens_adjust <- log(1.5); n_super <- 100000L

consistency_method <- "resample"
selection_rule     <- "neighborhood"
stop_threshold     <- NULL
hr_threshold       <- thr_c1
hr_consistency     <- thr_c2
pconsistency       <- 0.90
fs_splits          <- 400L
maxk               <- 2L
n_min              <- NULL
d0_min             <- 10L
d1_min             <- 10L
use_lasso <- FALSE; use_dina <- FALSE; use_grf <- FALSE; use_twostage <- TRUE
fs_conf_force      <- c("meno == 0", "er <= 0", "pgr <= 0")
fs_conf.cont_jcuts <- list(er = er_jcuts)
mult_law <- "poisson"                                     # production law (task section 6)

outcome_name <- "y_sim"; event_name <- "event_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"
confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")
seed_base  <- 8316951L
n_workers  <- min(.env_int("FS_S7_WORKERS", 60L),
                  max(1L, parallel::detectCores(logical = FALSE) - 1L))
z975 <- qnorm(0.975)

# The floors in force, as the template sets them (n_min NULL -> adaptive
# max(60, ceiling(0.10 n)); d0 / d1 10).  Asserted against every fit.
nmin_template <- max(60L, as.integer(ceiling(0.10 * n_sample)))
floors_id <- sprintf("nmin%d_d0%d_d1%d", nmin_template, d0_min, d1_min)

cat(sprintf("declcal cell %s (%s): dgm=%s hr=%.3f n=%d c1=%.2f c2=%.2f p*=%.2f focus=%s nbhd=%.2f B_cal=%d mode=%s cap=%.1fs workers=%d reps %d-%d floors=%s\n",
            cell_id, dgm_lab, dgm_model, target_hr_harm, n_sample, thr_c1, thr_c2,
            pconsistency, sg_focus, effect_neighborhood, B_cal, run_kind, cap_s,
            n_workers, sim_id_start, sim_id_start + n_sims - 1L, floors_id))

# ---- build-dgm (T:798-983) --------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
k_inter <- 1
k_treat <- 1
if (identical(dgm_model, "alt")) {
  k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm,
                               model = dgm_model, use_ahr = FALSE,
                               z1_quantile = harm_z1_quantile)
} else if (abs(target_hr_harm - 1) < 1e-12) {
  # Named change (1): the complete null.  See the header.
  k_treat <- 1e-10
} else {
  .hr_at <- function(kt)
    setup_gbsg_dgm(model = "null", k_treat = kt,
                   z1_quantile = harm_z1_quantile,
                   n_super = n_super, seed = seed_base)$hr_causal
  k_treat <- stats::uniroot(function(kt) .hr_at(kt) - target_hr_harm,
                            interval = c(0.1, 5), extendInt = "yes",
                            tol = 1e-12)$root
}
dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, k_treat = k_treat,
                      z1_quantile = harm_z1_quantile,
                      n_super = n_super, seed = seed_base)
harm_prevalence_super <- mean(dgm$df_super$flag_harm)
cat(sprintf("DGM: %s ; super-population prevalence %.4f ; k_inter %.5f ; k_treat %.3e ; hr_causal %.10f\n",
            dgm_model, harm_prevalence_super, k_inter, k_treat, dgm$hr_causal))
dgm <- compute_dgm_cde(dgm)
truth <- list(hr_causal = dgm$hr_causal %||% NA_real_,
              marg_H = dgm$hr_H_true %||% NA_real_, marg_Hc = dgm$hr_Hc_true %||% NA_real_)

if (identical(dgm_model, "null")) {
  eval_df <- fs_build_eval_frame(dgm, outcome_type = "survival",
                                 eval_seed = 20260628L,
                                 analysis_time = analysis_time,
                                 cens_adjust = cens_adjust)
  .tol <- 1e-8
  cat("\nSTRUCTURAL-NULL DESIGN-POINT GATE (tolerance 1e-8)\n")
  .prev <- mean(dgm$df_super$flag_harm)
  .lab_absent <- is.null(dgm$subgroup_info$fs_harm_true) &&
                 is.null(dgm$subgroup_info$grf_harm_true)
  cat(sprintf("  [1] planted-region prevalence %.10f ; truth labels absent %s -> %s\n",
              .prev, .lab_absent, if (.prev == 0 && .lab_absent) "PASS" else "FAIL"))
  .fsd <- get_FSdata(df.analysis = eval_df, use_lasso = FALSE, use_grf = FALSE,
                     confounders.name = intersect(confounders_base, names(eval_df)),
                     conf_force = fs_conf_force, conf.cont_jcuts = fs_conf.cont_jcuts,
                     outcome.name = outcome_name, event.name = event_name,
                     details = FALSE)
  .Z  <- as.matrix(.fsd$df[, .fsd$confs_names, drop = FALSE])
  .lp <- eval_df$loghr_po
  .lhr_unif <- unname(dgm$model_params$b_hr["treat"])
  .L  <- ncol(.Z)
  .idx <- c(lapply(seq_len(.L), function(i) i),
            if (maxk >= 2L && .L >= 2L)
              utils::combn(.L, 2L, simplify = FALSE) else list())
  .dev <- vapply(.idx, function(ii) {
    m <- if (length(ii) == 1L) .Z[, ii] == 1L else
      Reduce(`&`, lapply(ii, function(j) .Z[, j] == 1L))
    if (!any(m)) NA_real_ else abs(mean(.lp[m]) - .lhr_unif)
  }, numeric(1))
  .maxdev <- max(.dev, na.rm = TRUE)
  cat(sprintf("  [2] candidate family: %d factors, %d enumerated conjunctions; max_g |beta(g) - log HR_uniform| = %.3e -> %s\n",
              .L, length(.idx), .maxdev, if (.maxdev <= .tol) "PASS" else "FAIL"))
  .tol3_fb <- 1e-6
  .dev3 <- abs(dgm$hr_causal - target_hr_harm)
  .ok3_primary <- .dev3 <= .tol
  .ok3 <- .ok3_primary || .dev3 <= .tol3_fb
  cat(sprintf("  [3] super-population marginal Cox HR %.12f vs target %.12f ; |diff| %.3e -> %s\n",
              dgm$hr_causal, target_hr_harm, .dev3,
              if (.ok3_primary) "PASS (1e-8)" else if (.ok3) "PASS UNDER AMENDMENT (1e-6)" else "FAIL"))
  cat(sprintf("      uniform patient-level HR exp(b0[treat]) = %.12f ; k_treat = %.3e\n",
              exp(.lhr_unif), k_treat))
  if (!(.prev == 0 && .lab_absent) || !(.maxdev <= .tol) || !.ok3)
    stop("STRUCTURAL-NULL DESIGN-POINT GATE FAILED", call. = FALSE)
  cat("  GATE: PASS\n\n")
  rm(eval_df, .fsd, .Z, .lp)
}

# ---- per-replicate record ---------------------------------------------------
.schema <- c("cell_id", "dgm", "n", "c1", "c2", "rep", "seed",
             "G_pre", "G_post", "max_T_pre", "max_T_post",
             "kappa_hat_05", "kappa_hat_10", "pstar_implied_05", "pstar_implied_10",
             "alpha_FW_hat_1645", "alpha_FW_hat_1621",
             "Mstar_q90", "Mstar_q95", "Mstar_q99",
             "declared_conv", "declared_conv_exact", "declared_cal05", "declared_cal10",
             "n_admitted_conv", "n_admitted_cal05", "n_band",
             "sg_size_declared", "sg_size_argmax",
             "B_cal", "mult_law", "pconsistency_digits", "floors_id", "wall_sec", "status")
z_exact <- qnorm((1 + 0.90) / 2)          # 1.644854
z_round <- qnorm((1 + 0.895) / 2)         # 1.621 -- effective cutoff of the rounded rule

.row_template <- function(sim_id) {
  r <- as.list(setNames(rep(NA, length(.schema)), .schema))
  r$cell_id <- cell_id; r$dgm <- dgm_lab; r$n <- n_sample; r$c1 <- thr_c1; r$c2 <- thr_c2
  r$rep <- sim_id; r$seed <- seed_base + sim_id; r$B_cal <- B_cal; r$mult_law <- mult_law
  r$floors_id <- floors_id; r$status <- NA_character_
  r
}

.enum_family <- function(Z, maxk, n.min) {
  ns <- asNamespace("forestsearch")
  L <- ncol(Z); combo <- ns$generate_combination_indices(L, maxk)
  tot <- ns$calculate_max_combinations(L, maxk)
  fam <- list()
  for (kk in seq_len(tot)) {
    covs.in <- ns$get_covs_in(kk, maxk, L, combo$counts_1, combo$indices_1,
                              combo$counts_2, combo$indices_2,
                              combo$counts_3, combo$indices_3)
    k_sel <- sum(covs.in)
    if (k_sel < 1L || k_sel > maxk) next
    mem <- which(ns$get_subgroup_membership(Z, covs.in))
    if (length(mem) >= n.min)
      fam[[paste(colnames(Z)[covs.in == 1], collapse = " & ")]] <- mem
  }
  fam
}

record_replicate <- function(sim_id) {
  ns <- asNamespace("forestsearch")
  if (!isTRUE(get0(".declcal_traced", envir = .GlobalEnv, ifnotfound = FALSE))) {
    suppressMessages(trace("subgroup.search", where = ns, print = FALSE,
      tracer = quote(assign(".declcal_cap", list(Z = Z, Y = Y, Event = Event, Treat = Treat),
                            envir = .GlobalEnv))))
    assign(".declcal_traced", TRUE, envir = .GlobalEnv)
  }
  if (exists(".declcal_cap", envir = .GlobalEnv)) rm(".declcal_cap", envir = .GlobalEnv)
  r <- .row_template(sim_id)
  aux <- list(rep = sim_id, search_declared = NA_integer_, msg = NA_character_,
              wall_search = NA_real_, wall_field = NA_real_, n_unmatched = NA_integer_,
              replay_check = NA, nmin_fit = NA_integer_, d0_fit = NA_integer_,
              d1_fit = NA_integer_, digits_src = NA_character_)
  t0 <- proc.time()[3]

  df <- simulate_from_dgm(dgm, n = n_sample, analysis_time = analysis_time,
                          cens_adjust = cens_adjust, seed = seed_base + sim_id)
  df[[id_name]] <- seq_len(nrow(df))
  confs <- intersect(confounders_base, names(df))
  base_args <- list(
    df.analysis = df, outcome.name = outcome_name, event.name = event_name,
    treat.name = treat_name, id.name = id_name, flag_harm.name = harm_col,
    confounders.name = confs, is.RCT = TRUE, seedit = seed_base + sim_id,
    quiet = TRUE, sg_focus = sg_focus, subgroup_method = subgroup_method,
    hr.threshold = hr_threshold, hr.consistency = hr_consistency,
    pconsistency.threshold = pconsistency, n.min = n_min,
    selection_rule = selection_rule, effect_neighborhood = effect_neighborhood,
    stop_threshold = stop_threshold, parallel_args = list(plan = "sequential"),
    mr_inference = FALSE)
  method_args <- list(
    consistency_method = consistency_method,
    use_lasso = use_lasso, use_grf = use_grf, use_twostage = use_twostage,
    use_dina = use_dina, conf_force = fs_conf_force,
    conf.cont_jcuts = fs_conf.cont_jcuts, fs.splits = fs_splits, maxk = maxk,
    d0.min = d0_min, d1.min = d1_min)
  fit <- do.call(forestsearch, c(base_args, method_args))
  aux$wall_search <- proc.time()[3] - t0
  if (!inherits(fit, "forestsearch")) stop("forestsearch() returned no fit object")
  cap <- get0(".declcal_cap", envir = .GlobalEnv)
  if (is.null(cap)) stop("subgroup.search was not reached (no Z captured)")

  aca <- fit$args_call_all
  aux$nmin_fit <- as.integer(aca$n.min); aux$d0_fit <- as.integer(aca$d0.min)
  aux$d1_fit <- as.integer(aca$d1.min)
  if (!identical(aux$nmin_fit, nmin_template) || !identical(aux$d0_fit, d0_min) ||
      !identical(aux$d1_fit, d1_min))
    stop(sprintf("floors in force (n.min %s, d0 %s, d1 %s) differ from the template's %s",
                 aux$nmin_fit, aux$d0_fit, aux$d1_fit, floors_id))
  digits <- aca$pconsistency.digits
  aux$digits_src <- "fit$args_call_all"
  if (is.null(digits)) {
    digits <- eval(formals(ns$subgroup.consistency)$pconsistency.digits)
    aux$digits_src <- "subgroup.consistency() formal default (forestsearch() does not pass it)"
  }
  r$pconsistency_digits <- as.integer(digits)
  search_declared <- as.integer(!is.null(fit$sg.harm))
  aux$search_declared <- search_declared
  r$sg_size_declared <- if (search_declared == 1L)
    as.integer(sum(fit$grp.consistency$sg.harm.id == 1L, na.rm = TRUE)) else NA_integer_

  # ---- field capture on the pre-reduction family ----------------------------
  t1 <- proc.time()[3]
  fam <- .enum_family(cap$Z, maxk, aux$nmin_fit)
  if (!length(fam)) stop("empty pre-reduction family")
  dff <- data.frame(Y = cap$Y, Event = cap$Event, Treat = cap$Treat)
  gspec <- list(outcome_type = "survival", effect_measure = "HR",
                treat.name = "Treat", outcome.name = "Y", event.name = "Event",
                offset.name = NULL, adjust_covariates = NULL, adverse_outcome = TRUE)
  adm <- fit$admission
  if (is.null(adm$consistency)) stop("fit's admission set has no consistency floor")
  mr <- NULL
  for (sel_i in order(lengths(fam), decreasing = TRUE)[1:3]) {
    # The selected member only feeds MR's own de-biasing, which is not used;
    # it is taken FROM the family so the family is never augmented.
    mr <- ns$fs_mr_inference(
      df = dff, candidates = fam, spec = gspec, selected_members = fam[[sel_i]],
      admission = adm,
      reselection = ns$.fs_mr_reselection_from_focus(sg_focus, engine = "consistency"),
      effect_neighborhood = effect_neighborhood, selection_rule = selection_rule,
      draws = B_cal, multiplier = mult_law, include_complement = FALSE,
      ci_method = "ij", seed = seed_base + sim_id, return_reselection = FALSE,
      field_complement = FALSE, keep_declaration_field = TRUE,
      keep_field_matrix = FALSE)
    if (!is.null(mr$declaration_field)) break
  }
  fld <- mr$declaration_field
  if (is.null(fld)) stop("declaration field not captured")
  if (isTRUE(fld$meta$selected_appended)) stop("family was augmented by the selected member")
  aux$wall_field <- proc.time()[3] - t1

  c_cons <- adm$consistency$c_cons; p_star <- adm$consistency$p_star
  c_screen <- adm$effect_floor
  if (abs(p_star - 0.90) > 1e-12) stop("p_star is not 0.90")
  bh <- fld$beta_hat; sdv <- fld$sigma_D; Ms <- fld$Mstar
  T_pre <- (bh - c_cons) / sdv
  r$G_pre <- length(bh)
  r$max_T_pre <- max(T_pre)
  am <- which.max(T_pre)
  r$sg_size_argmax <- as.integer(fld$meta$sizes[am])

  q <- function(p, M = Ms) stats::quantile(M, p, type = 1, names = FALSE)
  r$kappa_hat_05 <- q(0.95); r$kappa_hat_10 <- q(0.90)
  r$pstar_implied_05 <- 2 * pnorm(r$kappa_hat_05) - 1
  r$pstar_implied_10 <- 2 * pnorm(r$kappa_hat_10) - 1
  r$alpha_FW_hat_1645 <- mean(Ms > z_exact)
  r$alpha_FW_hat_1621 <- mean(Ms > z_round)
  r$Mstar_q90 <- q(0.90); r$Mstar_q95 <- q(0.95); r$Mstar_q99 <- q(0.99)
  adm_at <- function(k) { fl <- c_cons + k * sdv; if (!is.null(c_screen)) fl <- pmax(c_screen, fl); bh >= fl }
  r$declared_cal05 <- as.integer(any(adm_at(r$kappa_hat_05)))
  r$declared_cal10 <- as.integer(any(adm_at(r$kappa_hat_10)))
  r$n_admitted_cal05 <- as.integer(sum(adm_at(r$kappa_hat_05)))

  # ---- post-reduction family: the one the executed screen evaluated --------
  red <- ns$.fs_decl_reduction(fit, fld$family_id)
  post <- if (is.null(red)) character(0) else red$screened
  if (!is.null(red)) { aux$n_unmatched <- as.integer(red$n_unmatched); aux$replay_check <- red$replay_check }
  r$G_post <- length(post)
  if (length(post)) {
    Tp <- T_pre[post]
    rate <- pmax(0, 2 * pnorm(Tp) - 1)
    adm_round <- round(rate, digits) >= p_star
    r$max_T_post <- max(Tp)
    r$declared_conv <- as.integer(any(adm_round))
    r$declared_conv_exact <- as.integer(any(Tp >= z_exact))
    r$n_admitted_conv <- as.integer(sum(adm_round))
    r$n_band <- as.integer(sum(rate >= 0.895 & rate < 0.900))
  } else {
    r$max_T_post <- NA_real_; r$declared_conv <- 0L; r$declared_conv_exact <- 0L
    r$n_admitted_conv <- 0L; r$n_band <- 0L
  }

  if (identical(run_kind, "pilot")) {
    for (Bp in pilot_B) {
      Mb <- Ms[seq_len(Bp)]            # sub-sampled from the one assembled field
      k05 <- q(0.95, Mb); k10 <- q(0.90, Mb)
      aux[[sprintf("kappa05_B%d", Bp)]] <- k05
      aux[[sprintf("kappa10_B%d", Bp)]] <- k10
      aux[[sprintf("cal05_B%d", Bp)]] <- as.integer(any(adm_at(k05)))
      aux[[sprintf("cal10_B%d", Bp)]] <- as.integer(any(adm_at(k10)))
    }
  }
  r$wall_sec <- proc.time()[3] - t0
  r$status <- "ok"
  list(row = r, aux = aux)
}

.safe_record <- function(s) {
  t0 <- proc.time()[3]
  if (cap_s > 0) setTimeLimit(elapsed = cap_s, transient = TRUE)
  out <- tryCatch(record_replicate(s), error = function(e) {
    msg <- conditionMessage(e)
    r <- .row_template(s)
    r$status <- if (grepl("reached elapsed time limit", msg, fixed = TRUE)) "abort_time" else "error"
    r$wall_sec <- proc.time()[3] - t0
    list(row = r, aux = list(rep = s, msg = msg))
  })
  setTimeLimit(elapsed = Inf)
  out
}

.bind <- function(lst, key) data.table::rbindlist(lapply(lst, `[[`, key), fill = TRUE)

# ---- run loop (T:1689-1705), chunked for the per-cell gate -------------------
sim_ids <- sim_id_start:(sim_id_start + n_sims - 1L)
chunks <- split(sim_ids, ceiling(seq_along(sim_ids) / chunk_n))
gate_max <- floor(0.01 * n_sims)                          # > 1% trips the cell
plan("sequential"); gc()
plan("multisession", workers = n_workers)
t_all <- proc.time()[3]
res <- list(); cell_status <- "complete"
.save <- function(final) {
  rows <- .bind(res, "row"); aux <- .bind(res, "aux")
  if (nrow(rows)) data.table::setcolorder(rows, .schema)
  payload <- list(
    results = as.data.frame(rows), aux = as.data.frame(aux),
    meta = list(task = "TASK_declcal_CAMPAIGN_2026-09-22_v2", cell_id = cell_id,
                dgm_label = dgm_lab, dgm_model = dgm_model, target_hr = target_hr_harm,
                k_inter = k_inter, k_treat = k_treat, harm_z1_quantile = harm_z1_quantile,
                harm_prevalence_super = harm_prevalence_super, truth = truth,
                n = n_sample, c1 = thr_c1, c2 = thr_c2, p_star = pconsistency,
                sg_focus = sg_focus, effect_neighborhood = effect_neighborhood,
                selection_rule = selection_rule, er_jcuts = er_jcuts, maxk = maxk,
                floors = list(n_min = "NULL -> max(60, ceiling(0.10 n))",
                              n_min_resolved = nmin_template, d0_min = d0_min, d1_min = d1_min),
                floors_id = floors_id, family = "prereduction (MR enumeration: n.min applied; d0/d1 not)",
                B_cal = B_cal, mult_law = mult_law, run_kind = run_kind,
                pilot_B = if (identical(run_kind, "pilot")) pilot_B else NULL,
                cap_s = cap_s, chunk_n = chunk_n, gate_max = gate_max,
                seed_convention = "data seed = seed_base + rep (simulate_from_dgm); forestsearch seedit = seed_base + rep; MR multiplier seed = seed_base + rep (forestsearch()'s own default); one offset, seed_base = 8316951",
                seed_base = seed_base, sim_id_start = sim_id_start, n_sims = n_sims,
                z_exact = z_exact, z_round = z_round,
                n_workers = n_workers, cell_status = cell_status, final = final,
                elapsed_s = proc.time()[3] - t_all,
                forestsearch_version = as.character(utils::packageVersion("forestsearch")),
                forestsearch_built = utils::packageDescription("forestsearch")$Built,
                r_version = as.character(getRversion()),
                hostname = unname(Sys.info()[["nodename"]]), written_at = Sys.time()))
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(payload, out_path)
  payload
}
for (ci in seq_along(chunks)) {
  part <- foreach(s = chunks[[ci]], .errorhandling = "pass",
                  .options.future = list(packages = c("forestsearch", "survival", "data.table"),
                                         seed = TRUE)) %dofuture% {
    .safe_record(s)
  }
  bad <- vapply(part, function(p) inherits(p, "error") || !is.list(p) || is.null(p$row), logical(1))
  if (any(bad)) {
    for (i in which(bad)) {
      s <- chunks[[ci]][i]; r <- .row_template(s); r$status <- "error"
      part[[i]] <- list(row = r, aux = list(rep = s, msg = paste("worker failure:", conditionMessage(part[[i]]))))
    }
  }
  res <- c(res, part)
  st <- vapply(res, function(p) p$row$status, character(1))
  n_fail <- sum(st %in% c("abort_time", "error"))
  cat(sprintf("chunk %d/%d done: %d reps, elapsed %.0fs, fail %d (abort %d, error %d), gate max %d\n",
              ci, length(chunks), length(res), proc.time()[3] - t_all, n_fail,
              sum(st == "abort_time"), sum(st == "error"), gate_max))
  if (n_fail > gate_max) { cell_status <- "tripped"; .save(FALSE); break }
  .save(ci == length(chunks))
}
plan("sequential")
payload <- .save(TRUE)
rows <- payload$results; aux <- payload$aux
cat(sprintf("CELL %s %s: %d rows, elapsed %.0fs, status table: %s\n", cell_id, cell_status,
            nrow(rows), payload$meta$elapsed_s,
            paste(names(table(rows$status)), table(rows$status), collapse = ", ")))

# ---- fidelity gate (task section 2.2) -----------------------------------------
ok <- rows$status == "ok"
mis <- rows$rep[ok & rows$declared_conv != aux$search_declared[match(rows$rep, aux$rep)]]
cat(sprintf("FIDELITY: declared_conv (rounded, post-reduction) vs search indicator on %d ok reps: %d disagree%s\n",
            sum(ok), length(mis), if (length(mis)) paste0(" -> reps ", paste(mis, collapse = ",")) else ""))
cat(sprintf("RATES: conv %.4f exact %.4f cal05 %.4f cal10 %.4f ; median wall %.2fs ; median G_pre %d\n",
            mean(rows$declared_conv[ok]), mean(rows$declared_conv_exact[ok]),
            mean(rows$declared_cal05[ok]), mean(rows$declared_cal10[ok]),
            median(rows$wall_sec[ok]), as.integer(median(rows$G_pre[ok]))))
if (length(mis)) quit(status = 3L)
if (identical(cell_status, "tripped")) quit(status = 4L)
quit(status = 0L)
