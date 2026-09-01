# =============================================================================
# closeout_parallel_surv_2026-09-02.R -- section 4.3 of the bootstrap
# close-out task (dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md):
# ONE survival (F2 gbsg) bootstrap at nb_boots = 100 under
# plan = multisession, workers = 48 -- WALL ONLY.  Its reproducibility
# follows from the same mechanism (main-process index matrix + per-replicate
# streams) plus Settlement A; the ~7-minute sequential comparator is not
# paid.  Fixture block verbatim from the committed battery.
# =============================================================================
suppressPackageStartupMessages(library(forestsearch))
suppressPackageStartupMessages(library(survival))
D <- "dev/glm-continuous-sims"
cat(sprintf("installed forestsearch %s | %s\n",
            as.character(utils::packageVersion("forestsearch")), R.version.string))
cat(sprintf("git HEAD: %s\n", system("git rev-parse --short HEAD", intern = TRUE)))

# ---- F2 gbsg fixture (verbatim from medians_baseline_2026-09-02.R) ----------
df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
args_surv <- list(
  df.analysis = df_surv, outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
  conf_force = c("er <= 0", "pgr <= 0"),
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = FALSE,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  minp = 0.025, use_twostage = TRUE, stop_threshold = NULL,
  consistency_method = "resample", fs.splits = 1000, max.minutes = 3,
  vi.grf.min = NULL,
  show_candidate_summary = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L))

cat("---- survival bootstrap 100 reps: multisession48 (wall only)\n")
set.seed(999L)
fs_est <- do.call(forestsearch, args_surv)
t0 <- proc.time()[[3]]
fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs_est, nb_boots = 100L, seed = 8316952L,
                                                       details = FALSE,
                                                       parallel_args = list(plan = "multisession", workers = 48L,
                                                                            show_message = FALSE)))
wall <- proc.time()[[3]] - t0
cat(sprintf("  wall %.1f s; per replicate %.2f s; implied B = 1000 at 48 workers: %.1f min\n",
            wall, wall / 100, wall / 100 * 1000 / 60))
saveRDS(list(wall = wall, per_rep = wall / 100,
             implied_B1000_min = wall / 100 * 1000 / 60,
             pkg_version = as.character(utils::packageVersion("forestsearch")),
             built_at = Sys.time()),
        file.path(D, "closeout_parallel_surv_2026-09-02.rds"))
cat("written:", file.path(D, "closeout_parallel_surv_2026-09-02.rds"), "\n")
