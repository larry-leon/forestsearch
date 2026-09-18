# ============================================================================
# accept_directive_b_probe.R  --  Directive B resolution probe (no fits)
#
# Records, per cell, the estimand and thresholds that forestsearch() resolves
# BEFORE any search runs: resolved effect_measure, screening and consistency
# thresholds, and the comparison scale.
#
# The probe does not re-implement the resolution logic.  It extracts the real
# top-level expressions from body(forestsearch) -- the argument-resolution
# prefix plus the whole `if (outcome_type != "survival")` block that builds
# threshold_config -- and evaluates them in a function whose formals ARE
# forestsearch()'s formals.  missing(hr.threshold) therefore behaves exactly
# as in a real call, and the probe cannot drift from the source it measures.
#
# Usage:  Rscript dev/accept_directive_b_probe.R <out_prefix>
# ============================================================================

suppressMessages(devtools::load_all(".", quiet = TRUE))

out_prefix <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_prefix)) out_prefix <- "dev/directive_b_probe"

# -- Build the probe function from real source expressions -------------------
# 2      outcome_type <- match.arg(outcome_type)
# 12     resolution site 1 (the live one)
# 13     adverse_outcome default
# 14-17  user_set_* flags and effect.threshold/consistency.threshold aliasing
# 19-20  args_names / args_call_all <- mget(formals)   [capture point]
# 33     sg_focus = "maxeff" override block (touches pconsistency.threshold)
# 39-46  match.arg re-run, outcome/threshold validation, NULL inits
# 47     the GLM/survival block: resolution site 2, thresholds, threshold_config
# 48     admission_resolved (method-dependent; carries the comparison scale)
.probe_expr_idx <- c(2L, 12L:17L, 19L:20L, 33L, 39L:48L)

.build_probe_fn <- function() {
  b <- body(forestsearch)
  stopifnot(identical(as.character(b[[1L]]), "{"))
  exprs <- lapply(.probe_expr_idx, function(i) b[[i]])
  record <- quote(list(
    effect_measure        = if (is.null(effect_measure)) NA_character_
                            else effect_measure,
    threshold_config      = threshold_config,
    effect_threshold      = effect_threshold,
    consistency_threshold = consistency_threshold,
    admission_resolved    = admission_resolved,
    aca_effect_measure    = if (is.null(args_call_all$effect_measure))
                              NA_character_ else args_call_all$effect_measure,
    aca_hr_threshold      = args_call_all$hr.threshold,
    aca_hr_consistency    = args_call_all$hr.consistency
  ))
  fn <- as.function(c(formals(forestsearch),
                      as.call(c(list(as.name("{")), exprs, list(record)))))
  environment(fn) <- environment(forestsearch)
  fn
}

# -- Fixed-seed data (the test-file factories) -------------------------------
sys.source(file.path("tests", "testthat", "helper-synthetic-dgm.R"),
           envir = environment())
dat <- list(
  survival   = .make_survival_data(seed = 42L),
  binary     = .make_binary_data(seed = 42L),
  continuous = .make_continuous_data(seed = 42L),
  count      = .make_count_data(seed = 42L)
)
CONF <- list(survival   = c("age", "stage", "sex"),
             binary     = c("age", "biomarker_hi", "sex"),
             continuous = c("age", "biomarker_hi", "sex"),
             count      = c("age", "biomarker_hi", "sex"))

# -- Cells -------------------------------------------------------------------
# hr.threshold / hr.consistency are deliberately NOT passed: the probe measures
# the inherited defaults, so missing() must stay TRUE.
.cell_args <- function(outcome_type, effect_measure, subgroup_method) {
  a <- list(
    df.analysis      = dat[[outcome_type]],
    confounders.name = CONF[[outcome_type]],
    outcome_type     = outcome_type,
    subgroup_method  = subgroup_method,
    treat.name       = "treat",
    id.name          = "id",
    quiet            = TRUE,
    details          = FALSE
  )
  if (outcome_type == "survival") {
    a$outcome.name <- "time"; a$event.name <- "event"
  } else {
    a$outcome.name <- "y"
    a$event.name   <- if (outcome_type == "binary") "y" else NULL
  }
  if (outcome_type == "count") a$offset.name <- "ftime"
  if (!is.na(effect_measure)) a$effect_measure <- effect_measure
  a
}

cells <- list()
for (em in c(NA, "OR", "RD", "RR"))
  for (m in c("consistency", "dina", "grf"))
    cells[[length(cells) + 1L]] <- list(ot = "binary", em = em, m = m)
for (m in c("consistency", "dina", "grf"))
  cells[[length(cells) + 1L]] <- list(ot = "survival", em = NA, m = m)
for (em in c(NA, "MD"))
  cells[[length(cells) + 1L]] <- list(ot = "continuous", em = em,
                                      m = "consistency")
for (em in c(NA, "IRR"))
  cells[[length(cells) + 1L]] <- list(ot = "count", em = em,
                                      m = "consistency")

# -- Run ---------------------------------------------------------------------
probe_fn <- .build_probe_fn()

.fmt <- function(x) if (is.null(x)) NA_character_ else
  paste(format(x, digits = 15), collapse = ",")

rows <- lapply(cells, function(cl) {
  args <- .cell_args(cl$ot, cl$em, cl$m)
  res <- tryCatch(suppressWarnings(do.call(probe_fn, args)),
                  error = function(e) conditionMessage(e))
  if (is.character(res))
    return(data.frame(outcome_type = cl$ot,
                      effect_measure_in = if (is.na(cl$em)) "<unset>" else cl$em,
                      subgroup_method = cl$m, resolved = "<error>",
                      screening = NA_character_, consistency = NA_character_,
                      screening_natural = NA_character_,
                      consistency_natural = NA_character_, scale = NA_character_,
                      adm_effect_floor = NA_character_,
                      adm_c_cons = NA_character_, adm_p_star = NA_character_,
                      aca_measure = NA_character_, note = res,
                      stringsAsFactors = FALSE))
  tc <- res$threshold_config
  data.frame(
    outcome_type        = cl$ot,
    effect_measure_in   = if (is.na(cl$em)) "<unset>" else cl$em,
    subgroup_method     = cl$m,
    resolved            = tc$effect_measure,
    screening           = .fmt(tc$screening),
    consistency         = .fmt(tc$consistency),
    screening_natural   = .fmt(tc$screening_natural),
    consistency_natural = .fmt(tc$consistency_natural),
    scale               = tc$scale,
    adm_effect_floor    = .fmt(res$admission_resolved$effect_floor),
    adm_c_cons          = .fmt(res$admission_resolved$consistency$c_cons),
    adm_p_star          = .fmt(res$admission_resolved$consistency$p_star),
    aca_measure         = res$aca_effect_measure,
    note                = "",
    stringsAsFactors    = FALSE
  )
})
probe <- do.call(rbind, rows)

sha <- system("git rev-parse HEAD", intern = TRUE)
dirty <- length(system("git status --porcelain --untracked-files=no",
                       intern = TRUE)) > 0L
attr(probe, "sha") <- sha
attr(probe, "tracked_dirty") <- dirty

saveRDS(probe, paste0(out_prefix, ".rds"))
con <- file(paste0(out_prefix, ".txt"), open = "wt")
writeLines(c(sprintf("# Directive B resolution probe"),
             sprintf("# SHA: %s%s", sha, if (dirty) " (tracked tree DIRTY)" else ""),
             sprintf("# R: %s", R.version.string), ""), con)
utils::capture.output(print(probe, right = FALSE), file = con)
close(con)
cat(readLines(paste0(out_prefix, ".txt")), sep = "\n")
