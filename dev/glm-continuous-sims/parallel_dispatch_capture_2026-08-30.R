# =============================================================================
# parallel_dispatch_capture_2026-08-30.R -- pre/post reference capture for the
# consistency screen's dispatch change (task: dev/tasks/cc_task_parallel_dispatch_2026-08-30.md)
# -----------------------------------------------------------------------------
# Usage (repository root):
#   Rscript dev/glm-continuous-sims/parallel_dispatch_capture_2026-08-30.R pre      # before the R/ edit
#   Rscript dev/glm-continuous-sims/parallel_dispatch_capture_2026-08-30.R post     # after
#   Rscript dev/glm-continuous-sims/parallel_dispatch_capture_2026-08-30.R compare  # identity gates (S4)
#   Rscript dev/glm-continuous-sims/parallel_dispatch_capture_2026-08-30.R boot     # S5 bootstrap check (post)
# Writes dev/glm-continuous-sims/parallel_dispatch_<tag>_2026-08-30.rds
#
# Fixtures (verbatim from bootstrap_profile_2026-08-30.R):
#   cont -- MD40 DGM, one trial n = 500 (seed 8316951 + 1), the drivers' forestsearch() args
#   surv -- survival::gbsg, the effMaxSG application's args
# Five configurations (task S2):
#   C1 cont  plan sequential, workers 1        resample   no early stop
#   C2 cont  plan sequential, workers 1        resample   maxeffCons + stop_threshold = 0.95
#   C3 cont  plan sequential, workers 1        split, fs.splits = 50, set.seed(20260830) right before the call
#   C4 surv  plan sequential, workers 1        resample   no early stop
#   C5 cont  plan multisession, workers 2      resample   no early stop
# Captured per configuration: sg.harm; sort(which(sg.harm.id == 1)); hr.subgroups (rows);
#   grp.consistency counters; the full out_sg$result table and its top row (effect, CI, Pcons);
#   df.est$treat.recommend; wall-clock.  Comparisons are on subject membership, never on strings.
# =============================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
stopifnot(requireNamespace("speff2trial", quietly = TRUE))
D <- "dev/glm-continuous-sims"
tag <- commandArgs(trailingOnly = TRUE)[1]
# narrow task (cc_task_parallel_dispatch_narrow_2026-08-30.md): "pre6" captures C6 alone pre-change,
# "post_narrow" captures C1-C6 post-change, "compare_narrow" runs S5's gates, "boot_narrow" S6's check
stopifnot(tag %in% c("pre", "post", "compare", "boot", "pre6", "post_narrow", "compare_narrow", "boot_narrow"))
out_path <- function(t) file.path(D, sprintf("parallel_dispatch_%s_2026-08-30.rds", t))
cat(sprintf("forestsearch %s; HEAD %s; vi.grf.min default %s\n", as.character(utils::packageVersion("forestsearch")),
            system("git rev-parse --short HEAD", intern = TRUE), deparse(formals(forestsearch)$vi.grf.min)))

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
  for (v in c("hemo", "homo", "drugs", "race", "gender", "symptom", "str2")) actg_df[[v]] <- as.factor(actg_df[[v]])
  actg_df
}
PAY <- "quarto/simulations/actg175/continuous/mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds"
truth <- readRDS(PAY)$truth
dgm <- generate_glm_dgm(
  data = actg_frame(), factor_vars = paste0("z", 1:12), outcome_var = "cd4_change",
  treatment_var = "treat", outcome_type = "continuous", effect_measure = "MD",
  subgroup_vars = c("z1", "z2"), subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model = "alt", k_treat = 1, k_inter = truth$beta_inter, n_super = 5000L,
  seed = 8316951L, verbose = FALSE)
stopifnot(abs(dgm$hazard_ratios$harm_subgroup - truth$effect_Q) < 1e-9,
          abs(dgm$subgroup_info$proportion - truth$prevalence_Q) < 1e-9)
SEED_BASE <- 8316951L; SIM_ID <- 1L
df_cont <- simulate_from_glm_dgm(dgm, n = 500L, seed = SEED_BASE + SIM_ID)
confs_cont <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                "hemo", "homo", "drugs", "race", "gender", "symptom", "str2")
args_cont <- function(pa, consistency_method = "resample", fs.splits = 400L, stop_threshold = NULL) list(
  df.analysis = df_cont, confounders.name = confs_cont,
  outcome.name = "y_sim", treat.name = "treat_sim", id.name = "id",
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = 30, consistency.threshold = 10, pconsistency.threshold = 0.90,
  fs.splits = fs.splits, n.min = 60L, d0.min = 12L, d1.min = 12L, maxk = 2L,
  vi.grf.min = NULL, sg_focus = "maxeffCons", selection_rule = "neighborhood",
  effect_neighborhood = 0.10, stop_threshold = stop_threshold, consistency_method = consistency_method,
  conf.cont_jcuts = list(age = 10, preanti = 10),
  use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
  is.RCT = TRUE, adverse_outcome = FALSE, details = FALSE, quiet = TRUE,
  seedit = SEED_BASE + SIM_ID, parallel_args = pa, mr_inference = FALSE)
df_surv <- survival::gbsg
df_surv <- within(df_surv, { id <- seq_len(nrow(df_surv)); time_months <- rfstime / 30.4375; grade3 <- ifelse(grade == "3", 1, 0) })
args_surv <- function(pa) list(
  df.analysis = df_surv, outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
  conf_force = c("er <= 0", "pgr <= 0"),
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.10, mr_inference = FALSE,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  minp = 0.025, use_twostage = TRUE, stop_threshold = NULL,
  consistency_method = "resample", fs.splits = 1000, max.minutes = 3, vi.grf.min = NULL,
  show_candidate_summary = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE, parallel_args = pa)

SEQ1 <- list(plan = "sequential", workers = 1L)
CONFIGS <- list(
  C1 = list(args = args_cont(SEQ1), seed = NULL),
  C2 = list(args = args_cont(SEQ1, stop_threshold = 0.95), seed = NULL),
  C3 = list(args = args_cont(SEQ1, consistency_method = "split", fs.splits = 50L), seed = 20260830L),
  C4 = list(args = args_surv(SEQ1), seed = NULL),
  C5 = list(args = args_cont(list(plan = "multisession", workers = 2L)), seed = NULL),
  C6 = list(args = args_cont(list(plan = "multisession", workers = 1L)), seed = NULL))   # narrow task S3: must stay on the batched path
CONFIG_SET <- switch(tag, pre = , post = , compare = , boot = c("C1", "C2", "C3", "C4", "C5"),
                     pre6 = "C6", post_narrow = , compare_narrow = , boot_narrow = c("C1", "C2", "C3", "C4", "C5", "C6"))

capture <- function(fs, secs, warns) {
  gc <- fs$grp.consistency
  res <- gc$out_sg$result
  list(
    sg.harm = fs$sg.harm,
    members = sort(which(gc$sg.harm.id == 1L)),
    treat.recommend = fs$df.est$treat.recommend,
    hr.subgroups = as.data.frame(fs$find.grps$out.found$hr.subgroups),
    counters = gc[intersect(c("n_candidates_total", "n_candidates_evaluated", "n_passed",
                              "early_stop_triggered", "early_stop_candidate", "algorithm", "stop_threshold"), names(gc))],
    result = if (!is.null(res)) as.data.frame(res) else NULL,
    top = if (!is.null(res)) as.data.frame(res)[1, , drop = FALSE] else NULL,
    secs = secs, warnings = warns)
}
run_all <- function() {
  out <- list()
  for (nm in CONFIG_SET) {
    cf <- CONFIGS[[nm]]
    warns <- character(0)
    if (!is.null(cf$seed)) set.seed(cf$seed)
    gc(); t0 <- proc.time()[["elapsed"]]
    fs <- withCallingHandlers(do.call(forestsearch, cf$args),
                              warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    secs <- proc.time()[["elapsed"]] - t0
    out[[nm]] <- capture(fs, secs, warns)
    cat(sprintf("[%s] %s | %6.2f s | selected %s | members %d | hr.subgroups %d | evaluated %s / %s, passed %s, early_stop %s at %s | warnings: %s\n",
                tag, nm, secs, paste(fs$sg.harm, collapse = " & "), length(out[[nm]]$members), nrow(out[[nm]]$hr.subgroups),
                out[[nm]]$counters$n_candidates_evaluated, out[[nm]]$counters$n_candidates_total, out[[nm]]$counters$n_passed,
                out[[nm]]$counters$early_stop_triggered, out[[nm]]$counters$early_stop_candidate,
                if (length(warns)) paste(unique(warns), collapse = " || ") else "none"))
    if (!is.null(out[[nm]]$top)) { cat("     top row: "); print(out[[nm]]$top[, intersect(c("hr", "L(HR)", "U(HR)", "Pcons", "N", "E", "K"), names(out[[nm]]$top))], row.names = FALSE) }
  }
  out$meta <- list(tag = tag, when = Sys.time(), head = system("git rev-parse --short HEAD", intern = TRUE),
                   pkg_version = as.character(utils::packageVersion("forestsearch")))
  saveRDS(out, out_path(tag)); cat("written:", out_path(tag), "\n")
}

compare <- function(pre_tag = "pre", post_tag = "post") {
  A <- readRDS(out_path(pre_tag)); B <- readRDS(out_path(post_tag))
  if (post_tag == "post_narrow") { A6 <- readRDS(out_path("pre6")); A$C6 <- A6$C6 }
  cat(sprintf("pre: HEAD %s (%s) | post: HEAD %s (%s)\n", A$meta$head, A$meta$pkg_version, B$meta$head, B$meta$pkg_version))
  comps <- c("sg.harm", "members", "treat.recommend", "hr.subgroups", "counters", "result", "top")
  gate <- TRUE
  for (nm in CONFIG_SET) {
    a <- A[[nm]]; b <- B[[nm]]
    id <- vapply(comps, function(k) identical(a[[k]], b[[k]]), logical(1))
    cat(sprintf("\n[%s] identical: %s | wall %.2f s -> %.2f s (x%.2f) | warnings pre: %s | post: %s\n", nm,
                paste(sprintf("%s=%s", comps, ifelse(id, "T", "F")), collapse = " "), a$secs, b$secs, a$secs / b$secs,
                if (length(a$warnings)) paste(unique(a$warnings), collapse = " || ") else "none",
                if (length(b$warnings)) paste(unique(b$warnings), collapse = " || ") else "none"))
    if (nm != "C3") { if (!all(id)) gate <- FALSE }
    if (nm == "C3" && !all(id)) {
      # characterise: per-candidate Pcons differences, keyed by candidate membership columns (grp id), and the selection
      ra <- a$result; rb <- b$result
      cat(sprintf("  C3 characterisation: selected same membership: %s; result rows %d vs %d\n",
                  identical(a$members, b$members), nrow(ra), nrow(rb)))
      key <- intersect(c("g", "grp"), names(ra))[1]
      if (!is.na(key)) {
        m <- merge(ra[, c(key, "Pcons")], rb[, c(key, "Pcons")], by = key, suffixes = c("_pre", "_post"))
        d <- m$Pcons_post - m$Pcons_pre
        cat(sprintf("  candidates in both: %d (pre-only %d, post-only %d); Pcons diff: mean %+.4f, sd %.4f, max|d| %.4f; MC sd of a rate at fs.splits = 50 and p ~ 0.9: %.4f\n",
                    nrow(m), nrow(ra) - nrow(m), nrow(rb) - nrow(m), mean(d), sd(d), max(abs(d)), sqrt(0.9 * 0.1 / 50)))
        print(summary(d))
      }
      if (!identical(a$members, b$members)) gate <- FALSE
    }
  }
  if (post_tag == "post_narrow") {
    b6 <- B$C6; if (b6$secs < 15) { cat(sprintf("C6 post wall %.2f s: plain-loop order -- the short-circuit caught a parallel plan\n", b6$secs)); gate <- FALSE }
    else cat(sprintf("C6 post wall %.2f s: batched-path order, as required\n", b6$secs))
    # the first task's post-change capture: C1-C5 must be identical to it, C3 included
    if (file.exists(out_path("post"))) {
      P1 <- readRDS(out_path("post"))
      for (nm in c("C1", "C2", "C3", "C4", "C5")) {
        id <- vapply(comps, function(k) identical(P1[[nm]][[k]], B[[nm]][[k]]), logical(1))
        cat(sprintf("[%s] narrow post vs first task's post: %s\n", nm, if (all(id)) "identical on all seven" else paste("DIFFERS:", paste(comps[!id], collapse = " "))))
        if (!all(id)) gate <- FALSE
      }
    } else cat("first task's post-change capture not in the tree; comparison skipped\n")
  }
  cat(sprintf("\nGATE (%s): %s\n", if (post_tag == "post_narrow") "C1, C2, C4, C5, C6 identical; C3 identical or same selection; C6 batched-order; C1-C5 identical to the first post" else "C1, C2, C4, C5 identical; C3 identical or same selection", if (gate) "PASS" else "FAIL"))
  invisible(gate)
}

boot_check <- function(out_tag = "boot") {
  NB <- 20L
  fs <- do.call(forestsearch, CONFIGS$C1$args)
  gc(); t0 <- proc.time()[["elapsed"]]
  fb <- suppressWarnings(forestsearch_bootstrap_dofuture(fs, nb_boots = NB, seed = SEED_BASE + SIM_ID, details = FALSE,
                                                         parallel_args = list(plan = "sequential", workers = 1L)))
  wall <- proc.time()[["elapsed"]] - t0
  bt <- fb$results
  cat(sprintf("bootstrap nb_boots = %d, plan sequential: wall %.1f s; per replicate %.2f s; bootstrap's own tmins_search mean %.2f s, median %.2f s; identified %d / %d\n",
              NB, wall, wall / NB, 60 * mean(bt$tmins_search, na.rm = TRUE), 60 * median(bt$tmins_search, na.rm = TRUE), sum(!is.na(bt$tmins_search)), nrow(bt)))
  cat(sprintf("projection B = 1000, sequential: %.0f s = %.1f min (profile task: 32.2 s per replicate, 490-537 min)\n", 1000 * wall / NB, 1000 * wall / NB / 60))
  saveRDS(list(nb = NB, wall = wall, per_rep = wall / NB, results = bt, when = Sys.time(),
               head = system("git rev-parse --short HEAD", intern = TRUE)), out_path(out_tag))
  cat("written:", out_path(out_tag), "\n")
}

switch(tag, pre = run_all(), post = run_all(), compare = compare(), boot = boot_check(),
       pre6 = run_all(), post_narrow = run_all(), compare_narrow = compare("pre", "post_narrow"), boot_narrow = boot_check("boot_narrow"))
