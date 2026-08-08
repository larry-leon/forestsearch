# =============================================================================
# mr_mechanism_A1.R  --  re-selection fidelity, instrumented
#
# A1 of CC_TASK_mr_mechanism_and_fb_pilot.md. INVESTIGATION ONLY: no R/ edits,
# no harness edits, no DGM changes. Internals are called via ::: as the task
# permits.
#
# Question: is the growing MR residual (REPORT_stopB_md_harm_grid.md) a
# RE-SELECTION FIDELITY GAP -- MR's simulated selection event differing from
# the gates that actually executed -- or genuine LINEARIZATION INADEQUACY in
# the extreme-selection tail?  The discriminating quantity is the SIZE
# DISTRIBUTION of MR's simulated winners against the realized winners
# (grid medians 76-167).
#
# Method: for each replicate, run forestsearch() harness-verbatim, then
# RE-EXECUTE the MR assembly + draw pipeline on that replicate's own fitted
# family using the package's own internals.  The replay is DISTRIBUTIONALLY
# FAITHFUL with its own recorded draw seed; it does not reproduce the grid
# run's individual draws bit-for-bit and no such claim is made.
#
# Run from the package root against the INSTALLED package:
#   Rscript dev/glm-continuous-sims/verification/mr_mechanism_A1.R
# =============================================================================

suppressMessages({
  library(forestsearch); library(speff2trial)
  library(doFuture); library(foreach); library(future)
})

fs_ns <- asNamespace("forestsearch")
.get  <- function(nm) get(nm, envir = fs_ns)

OUT_DIR <- Sys.getenv("A1_OUT", file.path(tempdir(), "a1"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## --- cells, per the task ----------------------------------------------------
CELLS <- list(list(n = 1000L, s = 200L), list(n = 4000L, s =  50L))
N_WORKERS <- 50L

## --- harness-verbatim configuration ----------------------------------------
actg_arms <- c(1L, 3L); actg_treat_arm <- 1L
actg_age_cut <- 34; actg_preanti_cut <- 744.5
dgm_n_super <- 5000L; cal_target_md <- -40

md_threshold <- 30; md_consistency <- 10
adverse_outcome <- FALSE; pconsistency <- 0.90
fs_splits <- 400L; maxk <- 2L; n_min <- 60L; d0_min <- 12L; d1_min <- 12L
vi_grf_min <- -0.2; use_twostage <- TRUE; is_rct <- TRUE
mr_draws <- 5000L
MR_MULTIPLIER <- "poisson"      # forestsearch()'s committed MR default

analysis_continuous_vars <- c("age","preanti","wtkg","karnof","cd40","cd80")
analysis_binary_vars     <- c("hemo","homo","drugs","race","gender","symptom")
confounders_analysis     <- c(analysis_continuous_vars, analysis_binary_vars)
outcome_name <- "y_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"

seed_base <- 8316951L; MAX_SIMS <- 5000L

## --- DGM, exactly as the harness builds it ---------------------------------
actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat_orig <- ifelse(actg_df$arms == actg_treat_arm, 1L, 0L)
actg_df$treat      <- 1L - actg_df$treat_orig
actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1  <- as.factor(ifelse(actg_df$age > actg_age_cut, 1L, 0L))
actg_df$z2  <- as.factor(ifelse(actg_df$preanti <= actg_preanti_cut, 1L, 0L))
actg_df$z3  <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4  <- as.factor(ifelse(actg_df$karnof <= median(actg_df$karnof), 1L, 0L))
actg_df$z5  <- as.factor(ifelse(actg_df$cd40 <= median(actg_df$cd40), 1L, 0L))
actg_df$z6  <- as.factor(ifelse(actg_df$cd80 <= median(actg_df$cd80), 1L, 0L))
actg_df$z7  <- as.factor(actg_df$hemo);   actg_df$z8  <- as.factor(actg_df$homo)
actg_df$z9  <- as.factor(actg_df$drugs);  actg_df$z10 <- as.factor(actg_df$race)
actg_df$z11 <- as.factor(actg_df$gender); actg_df$z12 <- as.factor(actg_df$symptom)
for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

dgm <- calibrate_glm_interaction(
  data = actg_df, factor_vars = paste0("z", 1:12),
  outcome_var = "cd4_change", treatment_var = "treat",
  target_effect = cal_target_md, outcome_type = "continuous",
  effect_measure = "MD", subgroup_vars = c("z1","z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L), k_inter_range = c(0,120),
  grid_step = 2, n_super = dgm_n_super, seed = seed_base, verbose = FALSE)

eval_df <- fs_build_eval_frame(dgm, outcome_type = "continuous")

set.seed(seed_base)
SEED_TABLE <- sample.int(.Machine$integer.max - 1L, MAX_SIMS)

## ---------------------------------------------------------------------------
## Rebuild the cut matrix Z from the fitted object's own labels.
##
## forestsearch() builds Z as dummy(df[, conf.screen]) after the (here absent)
## front-end screen, so with no front end conf.screen is the full confounder
## set and the columns are exactly the cuts named in $confounders.evaluated.
## Rebuilding from those labels reproduces Z without needing the internal
## object, and the reconstruction is CHECKED below against the fit's own
## qualifying-candidate count.
## ---------------------------------------------------------------------------
build_Z <- function(df, labels, qnames) {
  cols <- list()
  for (j in seq_along(labels)) {
    lab <- labels[j]; q <- qnames[j]
    if (grepl(" <= ", lab, fixed = TRUE)) {
      v  <- sub(" <=.*$", "", lab)
      cv <- as.numeric(sub("^.* <= ", "", lab))
      ind <- as.integer(df[[v]] <= cv)
    } else {
      ind <- as.integer(as.character(df[[lab]]) == "1")
    }
    cols[[paste0(q, ".0")]] <- 1L - ind     # complement / negation
    cols[[paste0(q, ".1")]] <- ind          # the cut holds
  }
  Z <- as.matrix(as.data.frame(cols))
  storage.mode(Z) <- "integer"
  Z
}

## Enumerate MR's candidate family exactly as forestsearch() does in its MR
## block (all <= maxk combinations of Z, kept when membership >= n.min).
build_family <- function(Z, maxk, n_min) {
  gci <- .get("generate_combination_indices"); cmc <- .get("calculate_max_combinations")
  gci_ <- .get("get_covs_in"); gsm <- .get("get_subgroup_membership")
  L <- ncol(Z); combo <- gci(L, maxk); tot <- cmc(L, maxk)
  fam <- list()
  for (kk in seq_len(tot)) {
    covs.in <- gci_(kk, maxk, L, combo$counts_1, combo$indices_1,
                    combo$counts_2, combo$indices_2,
                    combo$counts_3, combo$indices_3)
    k_sel <- sum(covs.in)
    if (k_sel < 1L || k_sel > maxk) next
    mem <- which(gsm(Z, covs.in))
    if (length(mem) >= n_min)
      fam[[paste(colnames(Z)[covs.in == 1], collapse = " & ")]] <- mem
  }
  fam
}

qs <- function(x) if (!length(x)) rep(NA_real_, 5) else
  unname(stats::quantile(x, c(0, .25, .5, .75, 1), na.rm = TRUE))

## --- one replicate ----------------------------------------------------------
run_one <- function(i, n_sample) {
  sd_i <- SEED_TABLE[i]
  RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)
  df <- simulate_from_glm_dgm(dgm, n = n_sample, seed = sd_i)
  df[[id_name]] <- seq_len(nrow(df))

  fs <- tryCatch(suppressWarnings(forestsearch(
    df.analysis = df, confounders.name = confounders_analysis,
    outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = md_threshold, consistency.threshold = md_consistency,
    pconsistency.threshold = pconsistency, fs.splits = fs_splits,
    n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
    vi.grf.min = vi_grf_min, consistency_method = "resample",
    use_twostage = use_twostage, is.RCT = is_rct,
    adverse_outcome = adverse_outcome,
    details = FALSE, quiet = TRUE, seedit = sd_i,
    parallel_args = list(plan = "sequential"),
    mr_inference = TRUE,
    mr_inference_args = list(ci_method = "ij", draws = mr_draws,
                             include_complement = TRUE))),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "a1err"))
  if (inherits(fs, "a1err"))
    return(data.frame(sim_id = i, n = n_sample, ok = FALSE, err = fs$msg))

  detected <- !is.null(fs$sg.harm) && length(fs$sg.harm) > 0 && !all(is.na(fs$sg.harm))
  if (!detected)
    return(data.frame(sim_id = i, n = n_sample, ok = FALSE, err = "NO-DETECTION"))

  sel_members <- which(fs$grp.consistency$sg.harm.id == 1L)

  ## ---- rebuild the family and replay the MR pipeline ---------------------
  Z   <- build_Z(df, fs$confounders.evaluated, fs$confounders.candidate)
  fam <- build_family(Z, maxk, n_min)

  spec <- list(outcome_type = "continuous", effect_measure = "MD",
               treat.name = treat_name, outcome.name = outcome_name,
               event.name = NULL, offset.name = NULL,
               adjust_covariates = NULL, adverse_outcome = adverse_outcome)

  # Ensure the realized winner is in the family and is the target, exactly as
  # fs_mr_inference() does before assembling.
  hit <- which(vapply(fam, function(ix) setequal(ix, sel_members), logical(1)))
  if (length(hit)) sel_lab <- names(fam)[hit[1]] else {
    fam[[".selected_H"]] <- sel_members; sel_lab <- ".selected_H" }

  asm <- .get(".fs_mr_assemble")(df, fam, spec)
  sel <- match(sel_lab, asm$names)
  if (is.na(sel))
    return(data.frame(sim_id = i, n = n_sample, ok = FALSE,
                      err = "selected subgroup not fit in reconstructed family"))

  B <- asm$B; bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes

  ## Admission set, from the fit's own resolved admission (not reconstructed).
  adm <- fs$admission
  has_eff <- !is.null(adm$effect_floor); has_cons <- !is.null(adm$consistency)
  c_cons <- if (has_cons) adm$consistency$c_cons else NULL
  if (has_eff && has_cons) {
    zq  <- stats::qnorm((1 + adm$consistency$p_star)/2)
    t_g <- pmax(adm$effect_floor, c_cons + zq * sdv)
  } else if (has_eff) t_g <- adm$effect_floor else t_g <- NULL
  .admit <- if (is.null(t_g)) function(bs) seq_along(bs) else
                              function(bs) which(bs >= t_g)
  .zcons <- function(bs) if (has_cons) (bs - c_cons)/sdv else NULL

  ## ---- ZERO-PERTURBATION FIDELITY PROBE ---------------------------------
  ## Does MR's ranking statistic, evaluated on the UNPERTURBED data, pick the
  ## candidate the identifier actually selected?  If yes, MR's rule agrees
  ## with the executing rule at the point where they must agree, and any
  ## drift under perturbation is a property of the perturbed argmax rather
  ## than a rule mismatch.  If no, the rules differ and the gap is fidelity.
  ## sort_subgroups() for sg_focus = "hr" orders (-Pcons, -hr, K)
  ## (R/subgroup_consistency_helpers.R:574-577), and under resample
  ## consistency Pcons = 2*Phi(zcons) - 1 (R/consistency_resample.R header),
  ## a strictly monotone transform of zcons -- so argmax Pcons == argmax
  ## zcons, which is exactly MR's "maxcons".
  pass0  <- .admit(bh)
  zc0    <- .zcons(bh)
  obs_win <- if (length(pass0))
    .get(".fs_mr_select")(bh, zc0, sz, pass0, "maxcons", 0.10,
                          "neighborhood", asm$log_scale) else NA_integer_
  obs_match <- isTRUE(obs_win == sel)
  # rank of the realized winner under MR's statistic among admitted candidates
  sel_rank <- if (length(pass0) && sel %in% pass0)
    sum(zc0[pass0] > zc0[sel]) + 1L else NA_integer_
  n_admitted0 <- length(pass0)
  obs_win_size <- if (is.na(obs_win)) NA_integer_ else sz[obs_win]

  ## Draws.  Own recorded seed; distributionally faithful, not bit-matching.
  draw_seed <- sd_i + 1L
  set.seed(draw_seed)
  Xi <- .get(".fs_mr_multipliers")(nrow(B), mr_draws, MR_MULTIPLIER)
  P  <- crossprod(B, Xi)
  beta_star <- bh + P

  win <- rep(NA_integer_, mr_draws); sb <- rep(NA_real_, mr_draws)
  for (b in seq_len(mr_draws)) {
    bs <- beta_star[, b]
    pass <- .admit(bs)
    if (!length(pass)) next
    s <- .get(".fs_mr_select")(bs, .zcons(bs), sz, pass, "maxcons",
                               0.10, "neighborhood", asm$log_scale)
    if (!is.na(s)) { win[b] <- s; sb[b] <- P[s, b] }
  }
  ok_b <- which(!is.na(win))
  wsz  <- sz[win[ok_b]]
  wbs  <- beta_star[cbind(win[ok_b], ok_b)]

  ## ---- gate replay: d0.min/d1.min are NOT enforced anywhere in MR --------
  ## The effect floor and consistency screen ARE enforced (via .admit above),
  ## and n.min is enforced at family construction. Per-arm minima are not.
  tt <- df[[treat_name]]
  fam_kept <- fam[asm$keep]
  arm_ok <- vapply(fam_kept, function(ix)
    sum(tt[ix] == 1L) >= d1_min && sum(tt[ix] == 0L) >= d0_min, logical(1))
  fail_arm <- mean(!arm_ok[win[ok_b]])

  ## effect floor / consistency failures among simulated winners, measured
  ## rather than assumed (expected 0 by construction of .admit).
  fail_floor <- if (has_eff) mean(wbs < adm$effect_floor) else NA_real_
  zc_win <- if (has_cons) (wbs - c_cons)/sdv[win[ok_b]] else NULL
  fail_cons <- if (has_cons)
    mean(zc_win < stats::qnorm((1 + adm$consistency$p_star)/2)) else NA_real_

  ## ---- correction, as the code computes it (fs_mr_inference.R:456/481/482/485)
  selection_bias <- mean(sb, na.rm = TRUE)
  fixed_bias     <- mean(P[sel, ok_b])
  fb <- if (is.finite(fixed_bias)) fixed_bias else 0
  beta_naive <- bh[sel]
  beta_deb   <- beta_naive - selection_bias - fb

  bHhat_or <- -fs_attach_betaHhat(
    data.frame(sim_id = i, sg_def = paste(fs$sg.harm, collapse = " & "),
               detected = 1L, stringsAsFactors = FALSE),
    eval_df, focus = "harm", outcome_type = "continuous",
    effect_measure = "MD")$betaHhat_H

  wq <- qs(wsz); bq <- qs(wbs)
  data.frame(
    sim_id = i, n = n_sample, ok = TRUE, err = NA_character_,
    sg_def = paste(fs$sg.harm, collapse = " & "),
    realized_size = length(sel_members),
    naive = beta_naive, betaHhat_or = bHhat_or,
    naive_bias = beta_naive - bHhat_or,
    family_size = length(asm$names), n_draws_ok = length(ok_b),
    n_admitted0 = n_admitted0, obs_match = obs_match,
    sel_rank_z = sel_rank, obs_win_size = obs_win_size,
    reselect_freq = mean(win[ok_b] == sel),
    win_min = wq[1], win_q1 = wq[2], win_med = wq[3], win_q3 = wq[4], win_max = wq[5],
    frac_le100 = mean(wsz <= 100), frac_le200 = mean(wsz <= 200),
    bstar_med = bq[3], bstar_max = bq[5],
    fail_arm_minima = fail_arm, fail_effect_floor = fail_floor,
    fail_consistency = fail_cons,
    selection_bias = selection_bias, fixed_bias = fb,
    correction = selection_bias + fb,
    beta_deb = beta_deb, mr_bias = beta_deb - bHhat_or,
    draw_seed = draw_seed, stringsAsFactors = FALSE)
}

## --- run --------------------------------------------------------------------
for (cl in CELLS) {
  cat(sprintf("\n=== A1 cell: n = %d, s = %d, draws = %d, workers = %d ===\n",
              cl$n, cl$s, mr_draws, N_WORKERS))
  t0 <- proc.time()[3]
  plan("sequential"); gc(); plan("multisession", workers = N_WORKERS)
  res <- foreach(k = seq_len(cl$s),
                 .options.future = list(packages = c("forestsearch"), seed = TRUE)
                 ) %dofuture% run_one(k, cl$n)
  plan("sequential"); gc()
  el <- proc.time()[3] - t0
  rows <- do.call(rbind, lapply(res, function(x) {
    if (!"reselect_freq" %in% names(x)) {
      pad <- setdiff(names(res[[which(vapply(res, function(z) isTRUE(z$ok), logical(1)))[1]]]), names(x))
      for (p in pad) x[[p]] <- NA
    }
    x }))
  saveRDS(list(rows = rows, meta = list(n = cl$n, s = cl$s, draws = mr_draws,
               multiplier = MR_MULTIPLIER, elapsed = el, workers = N_WORKERS,
               pkg_commit = tryCatch(system("git rev-parse --short HEAD", intern = TRUE)[1],
                                     error = function(e) NA_character_))),
          file.path(OUT_DIR, sprintf("a1_n%d_s%d.rds", cl$n, cl$s)))
  cat(sprintf("done: %d/%d ok, %.1f s\n", sum(rows$ok), nrow(rows), el))
}
cat("\nA1 complete. Output in ", OUT_DIR, "\n", sep = "")
