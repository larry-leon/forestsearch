# =============================================================================
# overnight_funnel50_O1.R
#
# O1 of CC_TASK_md_overnight_funnel50.md, AS AMENDED by Larry's "O1-O3
# definitions and go" message (that message supersedes the document where they
# conflict).
#
# The document's O1 says "the SAME configuration the pilot ran, verbatim". The
# amendment replaces that: under the pilot's verbatim configuration the MR
# guard fires in ~2 ms on every replicate and there is no funnel to
# characterize at all (REPORT_blockers_1-3_ADDENDUM.md sections D.1, D.3).
# O1 therefore runs the MR-COMPATIBLE configuration:
#
#     use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE
#     mr_inference = FALSE        (MR stays OFF in O1)
#     everything else pilot-verbatim, including the seed table
#
# and each replicate executes under RNGkind("L'Ecuyer-CMRG") with an explicit
# set.seed(seed_for(i)), which section D.2 established is what reproduces the
# pilot's draws (CMRG gives n_true = 341 on replicate 1; Mersenne-Twister gives
# 334, and 341 is the bundle's value).
#
# CHARACTERIZATION ONLY. No harness edit, no fix applied, no R/ change, no DGM
# change, no re-pilot. The flags-off configuration here is a measurement
# instrument, not an applied fix: the harness on disk is untouched and still
# carries use_lasso <- TRUE at :125.
#
# Run from the package root against the INSTALLED package:
#   Rscript dev/glm-continuous-sims/verification/overnight_funnel50_O1.R
# =============================================================================

suppressMessages({
  library(forestsearch); library(speff2trial)
  library(doFuture); library(foreach); library(future)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

OUT_DIR <- Sys.getenv("O1_OUT", file.path(tempdir(), "o1"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUT_RDS <- file.path(OUT_DIR, "overnight_funnel50_O1.rds")

N_SIMS    <- 50L
N_SAMPLE  <- 1000L
N_WORKERS <- 25L

## --- pilot constants, transcribed from the harness qmd ---------------------
actg_arms        <- c(1L, 3L)
actg_treat_arm   <- 1L
actg_age_cut     <- 34
actg_preanti_cut <- 744.5
dgm_n_super      <- 5000L
cal_target_md    <- -40

md_threshold    <- 30
md_consistency  <- 10
adverse_outcome <- FALSE
pconsistency    <- 0.90
fs_splits <- 400L; maxk <- 2L; n_min <- 60L; d0_min <- 12L; d1_min <- 12L
vi_grf_min <- -0.2
use_twostage <- TRUE; is_rct <- TRUE

analysis_continuous_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
analysis_binary_vars     <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders_analysis     <- c(analysis_continuous_vars, analysis_binary_vars)

outcome_name <- "y_sim"; treat_name <- "treat_sim"
id_name <- "id"; harm_col <- "flag_harm"

seed_base <- 8316951L
MAX_SIMS  <- 5000L

## --- DGM, exactly as the qmd builds it -------------------------------------
build_dgm <- function() {
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
  calibrate_glm_interaction(
    data = actg_df, factor_vars = paste0("z", 1:12),
    outcome_var = "cd4_change", treatment_var = "treat",
    target_effect = cal_target_md, outcome_type = "continuous",
    effect_measure = "MD", subgroup_vars = c("z1", "z2"),
    subgroup_cuts = list(z1 = 1L, z2 = 1L), k_inter_range = c(0, 120),
    grid_step = 2, n_super = dgm_n_super, seed = seed_base, verbose = FALSE)
}

dgm <- build_dgm()

set.seed(seed_base)
SEED_TABLE <- sample.int(.Machine$integer.max - 1L, MAX_SIMS)
seed_for <- function(i) SEED_TABLE[i]

## --- helpers ---------------------------------------------------------------

# Oriented MD: adverse_outcome = FALSE negates Y internally, so the quantity
# every gate compares against the floor is the MD on the negated scale, where
# positive means harm.
oriented_md <- function(df, idx) {
  if (!length(idx)) return(c(n = 0, md = NA_real_, se = NA_real_))
  d <- df[idx, , drop = FALSE]
  if (sum(d[[treat_name]] == 1L) < 2L || sum(d[[treat_name]] == 0L) < 2L)
    return(c(n = length(idx), md = NA_real_, se = NA_real_))
  d$.y <- -d[[outcome_name]]
  m <- stats::lm(.y ~ get(treat_name), data = d)
  c(n = length(idx), md = unname(stats::coef(m)[2]),
    se = unname(summary(m)$coefficients[2, 2]))
}

numeric_cuts <- function(labels, var) {
  hit <- grep(paste0("^", var, " <= "), labels, value = TRUE)
  sort(unique(as.numeric(sub(paste0("^", var, " <= "), "", hit))))
}

# Nearest cut on each side of a boundary, and the nearest overall.
bracket_info <- function(cuts, boundary) {
  below <- cuts[cuts <= boundary]; above <- cuts[cuts > boundary]
  nb <- if (length(below)) max(below) else NA_real_
  na <- if (length(above)) min(above) else NA_real_
  nearest <- if (!length(cuts)) NA_real_ else cuts[which.min(abs(cuts - boundary))]
  list(below = nb, above = na, brackets = length(below) > 0 && length(above) > 0,
       nearest = nearest, gap = if (is.na(nearest)) NA_real_ else abs(nearest - boundary))
}

## --- one replicate ---------------------------------------------------------
run_one <- function(i) {
  sd_i <- seed_for(i)

  # D.2: the pilot's draws live under the worker RNG kind future installs.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(sd_i)

  df <- simulate_from_glm_dgm(dgm, n = N_SAMPLE, seed = sd_i)
  df[[id_name]] <- seq_len(nrow(df))
  n_true <- sum(df[[harm_col]] == 1L, na.rm = TRUE)

  t0 <- proc.time()[3]
  fs <- tryCatch(suppressWarnings(forestsearch(
    df.analysis = df, confounders.name = confounders_analysis,
    outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = md_threshold, consistency.threshold = md_consistency,
    pconsistency.threshold = pconsistency, fs.splits = fs_splits,
    n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
    vi.grf.min = vi_grf_min,
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,   # amended O1
    use_twostage = use_twostage, is.RCT = is_rct,
    adverse_outcome = adverse_outcome,
    details = FALSE, quiet = TRUE, seedit = sd_i,
    parallel_args = list(plan = "sequential"),
    mr_inference = FALSE)),                                  # MR OFF in O1
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "o1_error"))
  secs <- proc.time()[3] - t0

  base <- list(sim_id = i, seed = sd_i, n_true = n_true, secs = secs,
               errored = inherits(fs, "o1_error"),
               err_msg = if (inherits(fs, "o1_error")) fs$msg else NA_character_)

  if (inherits(fs, "o1_error")) return(list(row = base, grids = NULL))

  labs <- fs$confounders.evaluated
  qs   <- fs$confounders.candidate
  fc   <- fs$find.grps$filter_counts
  tc   <- fs$threshold_config

  age_cuts <- numeric_cuts(labs, "age")
  pre_cuts <- numeric_cuts(labs, "preanti")
  bi_age <- bracket_info(age_cuts, actg_age_cut)
  bi_pre <- bracket_info(pre_cuts, actg_preanti_cut)

  ## --- the TRUE region Q ---------------------------------------------------
  q_row <- oriented_md(df, which(df[[harm_col]] == 1L))

  ## --- the nearest grid-analogue conjunction -------------------------------
  # Q is {age > 34} AND {preanti <= 744.5}.  On the grid: the negation of the
  # age cut nearest 34, conjoined with the preanti cut nearest 744.5.
  ga_idx <- if (is.na(bi_age$nearest) || is.na(bi_pre$nearest)) integer(0)
            else which(df$age > bi_age$nearest & df$preanti <= bi_pre$nearest)
  ga_row <- oriented_md(df, ga_idx)

  ## --- single-boundary analogues (for the per-boundary clearance counts) ---
  age_only <- if (is.na(bi_age$nearest)) integer(0) else which(df$age > bi_age$nearest)
  pre_only <- if (is.na(bi_pre$nearest)) integer(0) else which(df$preanti <= bi_pre$nearest)
  age_row <- oriented_md(df, age_only)
  pre_row <- oriented_md(df, pre_only)

  ## --- the selected rule ---------------------------------------------------
  detected <- !is.null(fs$sg.harm) && length(fs$sg.harm) > 0 &&
              !all(is.na(fs$sg.harm))
  sel_idx <- if (!is.null(fs$grp.consistency$sg.harm.id))
               which(fs$grp.consistency$sg.harm.id == 1L) else integer(0)
  sel_row <- oriented_md(df, sel_idx)

  row <- c(base, list(
    n_cut_labels   = length(labs),
    n_qs           = length(qs),
    n_evaluated    = fc$n_evaluated       %||% NA_integer_,
    n_variance     = fc$n_passed_variance %||% NA_integer_,
    n_prevalence   = fc$n_passed_prevalence %||% NA_integer_,
    n_redundancy   = fc$n_passed_redundancy %||% NA_integer_,
    n_events       = fc$n_passed_events   %||% NA_integer_,
    n_sample_size  = fc$n_passed_sample_size %||% NA_integer_,
    n_model        = fc$n_passed_cox      %||% NA_integer_,
    n_effect_floor = fc$n_passed_hr       %||% NA_integer_,
    thr_screening  = tc$screening         %||% NA_real_,
    thr_consistency= tc$consistency       %||% NA_real_,
    thr_scale      = tc$scale             %||% NA_character_,
    adverse_exec   = isTRUE(fs$args_call_all$adverse_outcome),
    admission_floor= fs$admission$effect_floor %||% NA_real_,
    family_status  = fs$family_status     %||% NA_character_,
    sg_focus       = fs$sg_focus          %||% NA_character_,

    age_below = bi_age$below, age_above = bi_age$above,
    age_brackets = bi_age$brackets, age_nearest = bi_age$nearest,
    age_gap = bi_age$gap,
    pre_below = bi_pre$below, pre_above = bi_pre$above,
    pre_brackets = bi_pre$brackets, pre_nearest = bi_pre$nearest,
    pre_gap = bi_pre$gap,

    Q_n = q_row[["n"]], Q_md = q_row[["md"]], Q_se = q_row[["se"]],
    Q_clears = isTRUE(q_row[["md"]] > md_threshold),

    GA_n = ga_row[["n"]], GA_md = ga_row[["md"]], GA_se = ga_row[["se"]],
    GA_clears = isTRUE(ga_row[["md"]] > md_threshold),

    AGE_n = age_row[["n"]], AGE_md = age_row[["md"]],
    AGE_clears = isTRUE(age_row[["md"]] > md_threshold),
    PRE_n = pre_row[["n"]], PRE_md = pre_row[["md"]],
    PRE_clears = isTRUE(pre_row[["md"]] > md_threshold),

    detected = detected,
    sel_rule = if (detected) paste(fs$sg.harm, collapse = " & ") else NA_character_,
    sel_n    = sel_row[["n"]], sel_md = sel_row[["md"]], sel_se = sel_row[["se"]]))

  list(row = row,
       grids = list(sim_id = i, labels = labs,
                    age_cuts = age_cuts, pre_cuts = pre_cuts))
}

## --- run -------------------------------------------------------------------
cat(sprintf("O1: %d replicates, n = %d, %d workers\n", N_SIMS, N_SAMPLE, N_WORKERS))
cat("configuration: use_lasso=FALSE use_grf=FALSE use_dina=FALSE mr_inference=FALSE\n")
cat("RNG: L'Ecuyer-CMRG with explicit set.seed(seed_for(i)) per replicate\n\n")

t_all <- proc.time()[3]
plan("sequential"); gc(); plan("multisession", workers = N_WORKERS)
res <- foreach(s = seq_len(N_SIMS), .options.future =
                 list(packages = c("forestsearch"), seed = TRUE)) %dofuture%
  run_one(s)
plan("sequential"); gc()
elapsed <- proc.time()[3] - t_all

rows  <- do.call(rbind, lapply(res, function(x) as.data.frame(x$row,
                                                stringsAsFactors = FALSE)))
grids <- lapply(res, function(x) x$grids)

saveRDS(list(rows = rows, grids = grids,
             meta = list(n_sims = N_SIMS, n_sample = N_SAMPLE,
                         workers = N_WORKERS, elapsed_sec = elapsed,
                         rng_kind = "L'Ecuyer-CMRG",
                         config = "use_lasso=FALSE use_grf=FALSE use_dina=FALSE mr_inference=FALSE",
                         seed_base = seed_base,
                         md_threshold = md_threshold,
                         pkg_version = as.character(utils::packageVersion("forestsearch")),
                         pkg_commit = tryCatch(
                           system("git rev-parse --short HEAD", intern = TRUE)[1],
                           error = function(e) NA_character_),
                         r_version = as.character(getRversion()))),
        OUT_RDS)

cat(sprintf("\nDONE. %d rows in %.1f s (%.1f s/replicate wall, %d workers)\n",
            nrow(rows), elapsed, elapsed / N_SIMS, N_WORKERS))
cat("saved: ", OUT_RDS, "\n", sep = "")
cat("errored replicates: ", sum(rows$errored), "\n", sep = "")
