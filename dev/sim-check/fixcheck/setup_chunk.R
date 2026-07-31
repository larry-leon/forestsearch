.t0_doc <- Sys.time()    # start of document-level wall-clock (total reported at the very end)

# NOTE: forestsearch must be INSTALLED (devtools::install()), not load_all():
# the doFuture multisession workers spawn separate R processes that only see the
# installed package.
library(forestsearch)
library(survival)
library(data.table)
library(gt)
library(ggplot2)
library(foreach)
library(doFuture)
library(future)

# Single-source bundle enrichment: beta(Hhat)/beta(Hhat^c) truth targets at the
# realized region, plus covariate-name extraction for the identification figures.
# Must sit in the render working directory.  rule_covs() ships to the %dofuture%
# workers as a detected global of record_replicate (as .na_record/.classify do).
source("betaHhat_truth.R")

# Parallelization topology (works with run_mode / batching above):
#   "sims"  -> replicates fan out across workers; each replicate's forestsearch +
#              bootstrap run sequentially.  Best when n_sims >= n_workers (fills cores).
#   "boots" -> replicates run sequentially; the forestsearch search AND the nb_boots
#              bootstrap fan out across workers WITHIN each replicate.  Best when
#              n_sims is small (single-/few-trial deep dives) and "sims" would idle cores.
# Bootstrap resample indices are seeded on the main process, so results are
# plan-invariant: the two modes differ in WALL-CLOCK only, not in detection/estimates.
parallel_mode <- "boots"

# ── Smoke-test counts (the only knobs to change for production) ────────────
n_sims       <- 100L      # simulated trials                (prod: 1000 -> 20000)
diag_mode    <- FALSE     # DIAGNOSTIC: TRUE caps the run to 5 replicates and relays each
                          # worker's forestsearch() error (otherwise swallowed by the inner
                          # tryCatch).  Set FALSE for a real run.
nb_boots     <- 300L      # FB nonparametric bootstraps (prod: 200 -> 500).
                          # Set 0L (or NULL) to SKIP FB entirely and run only
k_random_noise <- 0       # Add k_random_noise confounders

# ── Batching across renders (split one long run into seed-disjoint pieces) ──
# Each render is ONE batch of `n_sims` replicates starting at `sim_id_start`.
# Every per-replicate seed is `seed_base + sim_id` (DGM, search, bootstrap, AND
# the noise covariates), so disjoint sim_id ranges across renders are independent
# and pool to exactly what one big run would have produced:
#
#   Batch 1 : run_mode "batch", sim_id_start = 1     -> ..._res_1_100.rds
#   Batch 2 : run_mode "batch", sim_id_start = 101   -> ..._res_101_200.rds   (re-render)
#   ...      (re-render per batch, advancing sim_id_start by n_sims each time)
#   Combine : run_mode "combine"                     -> reads ..._res_*.rds, rbinds,
#                                                       renders every table on the union
#
# DEFAULT is a single batch starting at 1: render once and you have the whole run.
run_mode     <- "batch"   # "batch" (run + save one piece) | "combine" (merge + report)
sim_id_start <- 1L        # first sim_id of THIS batch (101L, 201L, ... for later pieces)

# ── The identifier: engine + focus ─────────────────────────────────────────
# These two together ARE the estimator, and both feed the output stem below,
# so they are defined here rather than with the search-tuning knobs.  Defining
# sg_focus later would leave focus_tag undefined at stem-construction time.
subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
sg_focus           <- "maxcons"    # consistency argmax: highest consistency rate,
                                   # ties broken by effect.  Alias of "eff" -- both
                                   # normalize to one canonical focus, so this rename
                                   # is inert.  "maxcons" names what the rule does.

target_hr_harm  <- 1.0    # calibrate k_inter to this Cox HR in the harm subgroup
n_sample        <- 500L   # M1 sample size

# Output stem derives from the knobs above, so changing n_sample / target_hr_harm /
# k_random_noise / subgroup_method / sg_focus renames every saved file automatically
# -- no manual edit, and runs under a different engine or focus never overwrite each
# other. The leading token IS the subgroup method: "fs" for consistency, else
# "dina"/"grf"; focus_tag follows it. Batch range is appended below.
#
# NOT back-compatible: the stem changed from "<engine>_t1_t2_m1_..." to
# "<engine>_<focus>_fb_mr_m1_...".  Batch .rds files written under the old stem
# are invisible to combine_glob -- re-run them or rename them to match.
method_tag   <- if (identical(subgroup_method, "consistency")) "fs" else subgroup_method

# The IDENTIFIER is the estimator here; FB and MR are resampling procedures applied
# to it.  The identifier has two parts -- the engine (method_tag) and the selection
# focus -- and only the engine was in the stem, so re-running a cell under a
# different sg_focus silently overwrote the same .rds, and combine_glob pooled the
# two together.  focus_tag closes that.
#
# Aliases normalize FORWARD to the descriptive name rather than back through
# .normalize_sg_focus(), which maps eff -> "hr": the Cox-flavoured internal
# vocabulary.  "eff*" generalizes across GLM outcomes, "hr" does not.
#
# EVERY spelling of one rule must land on one tag, or the stem splits a single
# estimator across two pools -- the same failure focus_tag exists to prevent,
# just one level down.  forestsearch() collapses hr/eff/maxcons to the canonical
# "hr" internally (.normalize_sg_focus), so all three ARE one rule and must
# share a stem; likewise hrMaxSG/effMaxSG and hrMinSG/effMinSG.  The map below
# mirrors .normalize_sg_focus()'s alias table in the forward direction.
# maxSG / minSG / maxeff / maxeffCons are NOT aliases (see that helper's
# "NOT aliases" note) and pass through unchanged.
focus_tag    <- switch(sg_focus,
                       hr      = , eff      = , maxcons = "maxcons",
                       hrMaxSG = , effMaxSG = "effMaxSG",
                       hrMinSG = , effMinSG = "effMinSG",
                       sg_focus)

rds_stem     <- sprintf("%s_%s_fb_mr_m1_h%02d_knoise%d_n%d",
                        method_tag, focus_tag, round(10 * target_hr_harm),
                        k_random_noise, n_sample)

sim_id_end   <- sim_id_start + n_sims - 1L
rds_path     <- sprintf("%s_res_%d_%d.rds", rds_stem, sim_id_start, sim_id_end)  # THIS batch's file
combine_glob <- sprintf("%s_res_*.rds", rds_stem)   # combine mode: which batch files to merge
combine_files <- NULL   # optional explicit character vector; NULL -> use combine_glob
save_combined <- TRUE   # combine mode: ALSO write the pooled bundle to disk
combined_path <- NULL   # NULL -> auto "<rds_stem>_combined_<min>_<max>.rds" (deliberately
                        # outside combine_glob, so a re-combine never re-ingests the pool)

# ── Persist results? The FB bootstrap (B x n_sims, each re-running the
# engine) is the wall-clock bottleneck; saving the results bundle lets you
# regenerate the summary tables without re-running the simulation. 

# Scale-up scenarios for the timing projection (n_sims, nb_boots). Edit freely.
proj_scenarios <- list(
  c(n_sims = 100L,  nb_boots = 300L),
  c(n_sims = 500L,  nb_boots = 300L),
  c(n_sims = 10000L, nb_boots = 500L)
)
                          # the MR -- see run_fb below.
mr_draws   <- 5000L     # MR multiplier draws         (keep ~2000)

save_rds <- identical(run_mode, "batch") && !is.null(rds_path)  # write only when running a batch

# FB is the wall-clock bottleneck.  forestsearch_bootstrap_dofuture() itself
# requires nb_boots >= 1 (it errors otherwise), so a skip is enforced HERE -- the
# bootstrap is simply never called -- rather than by passing 0L/NULL into it.
# nb_boots = 0L (or NULL) => run_fb = FALSE => MR-only run.
run_fb <- is.numeric(nb_boots) && length(nb_boots) == 1L && nb_boots >= 1L
run_smoke    <- TRUE      # master flag; FALSE short-circuits the run loop


# ── DGM: León M1 cell (GBSG AFT super-population, alt hypothesis) ───────────
dgm_model       <- "alt"
analysis_time   <- 84     # administrative censoring (months)
cens_adjust     <- log(1.5)
n_super         <- 100000L  # fixed super-population = the eligible-patient universe; trials AND beta(Hhat) draw from this one pool

# ── Subgroup-identification method ─────────────────────────────────────────
# The FB bootstrap (H2, León Eq.7) re-runs whichever engine fs.est used,
# regenerating the candidate family on each resample, so the FB-vs-MR
# comparison is produced for every engine.  Interpretation differs by engine:
# for "consistency" the family is fixed, so FB and MR target the same
# estimand and should agree (the dfbeta-linearisation check); for "dina"/"grf"
# the family regenerates per resample, so FB is the *unconditional*
# comparator against MR's *conditional-on-family* estimate, and a
# FB<->MR gap is the family-regeneration signal, not an error.

grf_selection   <- "frontier"      # only if subgroup_method == "grf": "tree" | "frontier"
grf_select_statistic <- "effect"   # only if grf frontier: "dr" | "effect" (no-op for tree)
grf_depth       <- 2L              # GRF policy-tree depth / max frontier conjunctions
dmin.grf        <- 0.0             # GRF DR-score screening floor
dina_args       <- list()          # extra args forwarded to the DINA selector (dina_args =)
dina_select_statistic <- "effect"  # only if subgroup_method == "dina": "dina" | "effect"

# ── Search configuration (León tuning; resample identification) ────────────
# sg_focus and subgroup_method are set earlier, with the output stem they feed.
consistency_method <- "resample"   # decided: resample identification (the new default)
# MR's de-biasing rule is no longer set explicitly: forestsearch() derives it from
# sg_focus ("maxcons" for the consistency engine), reproducing the legacy
# value while staying faithful per engine.  Override via `reselection =` in
# mr_inference_args (advanced use only).
selection_rule     <- "neighborhood"  # effMaxSG/effMinSG band rule.  INERT under
                                      # sg_focus = "maxcons", whose key is
                                      # consistency-primary with no band term.
                                      # Retained unchanged so the effMaxSG/effMinSG
                                      # arms can use it.
effect_neighborhood <- 0.10           # package default; band half-width for
                                      # effMaxSG/effMinSG.  Pinned rather than
                                      # inherited so the value is visible; inert
                                      # under sg_focus = "maxcons".
hr_threshold       <- 0.90
hr_consistency     <- 0.80
pconsistency       <- 0.90
fs_splits          <- 400L         # moot under "resample" (closed-form, not literal split)
maxk               <- 2L
n_min              <- NULL
d0_min             <- 10L
d1_min             <- 10L
use_lasso          <- FALSE
use_dina           <- FALSE
use_grf            <- FALSE         # FALSE = LASSO-only (fast, "FSl"); set TRUE for "FSlg"
use_twostage       <- TRUE

fs_conf_force      <- c("meno == 0", "er <= 0", "pgr <= 0")   # premenopausal (raw analog of z3 == 1)
fs_conf.cont_jcuts <- list(er = 10)     # J-quantile grid on raw ER (analog of z1; now meaningful, ER is continuous)

# ── Column names emitted by simulate_from_dgm() ────────────────────────────
outcome_name   <- "y_sim"
event_name     <- "event_sim"
treat_name     <- "treat_sim"
id_name        <- "id"
harm_col       <- "flag_harm"
# Raw GBSG covariate candidates (as in a real analysis): forestsearch discovers
# its own cuts. Continuous -> er, age, pgr, nodes, size; categorical -> meno (0/1),
# grade (1/2/3). The pre-dichotomized z1..z5 are NOT used. The true harm subgroup
# (er <= 25th pct & meno == 0) must be re-discovered from these raw variables.

confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")

# ── Reproducibility / parallelism ──────────────────────────────────────────
seed_base  <- 8316951L
n_workers  <- ceiling(0.90*max(1L, min(parallel::detectCores(logical = FALSE) - 1L)))

stopifnot(parallel_mode %in% c("sims", "boots"))
inner_parallel <- if (identical(parallel_mode, "boots"))
  list(plan = "multisession", workers = n_workers) else list(plan = "sequential")
cat(sprintf("Parallel mode: %s (%s)\n", parallel_mode,
            if (identical(parallel_mode, "boots")) "replicates sequential, bootstrap parallel"
            else "replicates parallel, bootstrap sequential"))
z975       <- qnorm(0.975)

cat(sprintf("Smoke test: method=%s%s, n_sims=%d, nb_boots=%s, mr_draws=%d, workers=%d\n",
            subgroup_method,
            if (identical(subgroup_method, "grf")) paste0("/", grf_selection) else "",
            n_sims, if (is.null(nb_boots)) "0 (skip)" else as.character(nb_boots),
            mr_draws, n_workers))
