#!/usr/bin/env Rscript
# =============================================================================
# §5A evidence probe -- stage 2
#   Q4  three engines (consistency / dina / grf)
#   Q5  no subgroup identified
#   Q6  sg_focus = "maxeff"
#   Q7  GLM path with consistency_method != "resample"  (silent NULL)
#   Q8  reselection rule per sg_focus  (.fs_dg_reselection_from_focus)
#   Q9  gate family vs identifier admissible set (n.min/d0.min/d1.min/maxk,
#       max_subgroups_search)
#   Q10 does the gate run inside bootstrap replicates / CV folds?
# READ-ONLY with respect to the package.  Writes only under dev/terminology-work/.
# =============================================================================

suppressMessages(library(forestsearch))
source("tests/testthat/helper-synthetic-dgm.R")
dir.create("dev/terminology-work/out", showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a
hr  <- function(s) cat("\n===== ", s, " =====\n", sep = "")
ns  <- asNamespace("forestsearch")
get_ <- function(nm) get(nm, envir = ns)

CONF <- c("age", "stage", "sex", "noise")
d    <- .make_survival_data(N = 400L, HR_harm = 3.0)
mk <- function(...) .fs_args_for("survival", confounders = CONF,
  extra = utils::modifyList(
    list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE), list(...)))

DGA <- list(draws = 200L, seed = 11L)

# ---------------------------------------------------------------------------
hr("Q8: .fs_dg_reselection_from_focus() for every accepted sg_focus")
rf   <- get_(".fs_dg_reselection_from_focus")
foci <- c("hr", "eff", "maxcons", "maxeff", "maxeffCons",
          "maxSG", "minSG", "hrMaxSG", "hrMinSG", "effMaxSG", "effMinSG")
map <- data.frame(
  sg_focus    = foci,
  consistency = vapply(foci, function(f)
    tryCatch(rf(f, engine = "consistency"), error = function(e) "<error>"), ""),
  effect      = vapply(foci, function(f)
    tryCatch(rf(f, engine = "effect"),      error = function(e) "<error>"), ""),
  row.names = NULL, stringsAsFactors = FALSE)
print(map)

# ---------------------------------------------------------------------------
hr("Q4/Q6: engines x sg_focus -- did the gate run, and with which reselection?")
cases <- list(
  list(lab = "consistency / maxSG",  args = mk(subgroup_method = "consistency", sg_focus = "maxSG")),
  list(lab = "consistency / maxeff", args = mk(subgroup_method = "consistency", sg_focus = "maxeff")),
  list(lab = "consistency / maxeffCons", args = mk(subgroup_method = "consistency", sg_focus = "maxeffCons")),
  list(lab = "consistency / hr",     args = mk(subgroup_method = "consistency", sg_focus = "hr")),
  list(lab = "dina / maxSG",         args = mk(subgroup_method = "dina",  sg_focus = "maxSG")),
  list(lab = "grf / maxSG",          args = mk(subgroup_method = "grf",   sg_focus = "maxSG", use_grf = TRUE))
)
rows <- list()
for (cs in cases) {
  a  <- utils::modifyList(cs$args, list(debias_gate = TRUE, debias_gate_args = DGA))
  rr <- tryCatch(.run_fs_capture(d, a),
                 error = function(e) list(value = NULL, warnings = paste("ERROR:", conditionMessage(e))))
  fs <- rr$value; g <- fs$debias_gate
  rows[[length(rows) + 1L]] <- data.frame(
    case        = cs$lab,
    found       = !is.null(fs$sg.harm) && length(fs$sg.harm) > 0,
    gate_ran    = !is.null(g),
    reselection = g$gate$reselection %||% NA_character_,
    n_family    = g$n_family    %||% NA_integer_,
    n_selected  = g$n_selected  %||% NA_integer_,
    harm_flag   = format(g$harm_flag %||% NA),
    hfd         = format(fs$harm_flag_debiased),
    warns       = if (!length(rr$warnings)) "" else
                    substr(paste(unique(rr$warnings), collapse = " | "), 1, 60),
    stringsAsFactors = FALSE)
}
eng <- do.call(rbind, rows); print(eng, right = FALSE)
saveRDS(eng, "dev/terminology-work/out/engines_x_focus.rds")

# ---------------------------------------------------------------------------
hr("Q5: no subgroup identified -- does the gate run, skip, or error?")
d_null <- .make_survival_data(N = 300L, HR_harm = 1.0)   # no harm subgroup
a_null <- utils::modifyList(mk(hr.threshold = 3.0, pconsistency.threshold = 0.99),
                            list(debias_gate = TRUE, debias_gate_args = DGA))
rn <- tryCatch(.run_fs_capture(d_null, a_null),
               error = function(e) list(value = NULL, warnings = paste("ERROR:", conditionMessage(e))))
fsn <- rn$value
cat("sg.harm NULL?      :", is.null(fsn$sg.harm), "\n")
cat("debias_gate NULL?  :", is.null(fsn$debias_gate), "\n")
cat("harm_flag_debiased :", format(fsn$harm_flag_debiased), "\n")
cat("errored?           :", is.null(fsn), "\n")
cat("warnings           :", if (!length(rn$warnings)) "(none)" else
      substr(paste(unique(rn$warnings), collapse = " | "), 1, 200), "\n")

# ---------------------------------------------------------------------------
hr("Q7: the .dg_glm_ok / .dg_cox_ok eligibility test on the consistency path")
cat("Source (forestsearch_main.R:2817-2821):\n")
cat("  .dg_glm_ok <- consistency_method == 'resample' && !is.null(estimator_fn)\n")
cat("  .dg_cox_ok <- outcome_type == 'survival' && is.null(estimator_fn)\n\n")
db <- .make_binary_data(N = 400L, OR_harm = 3.0)
for (cm in c("resample", "sequential")) {
  ab <- .fs_args_for("binary", confounders = c("age", "biomarker_hi"),
    extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
                 consistency_method = cm,
                 debias_gate = TRUE, debias_gate_args = DGA))
  rb <- tryCatch(.run_fs_capture(db, ab),
                 error = function(e) list(value = NULL,
                   warnings = paste("ERROR:", conditionMessage(e))))
  fb <- rb$value
  cat(sprintf("consistency_method=%-11s found=%-5s gate_ran=%-5s hfd=%-5s warn=%s\n",
      cm, !is.null(fb$sg.harm) && length(fb$sg.harm) > 0,
      !is.null(fb$debias_gate), format(fb$harm_flag_debiased),
      if (!length(rb$warnings)) "(none)" else
        substr(paste(unique(rb$warnings), collapse = " | "), 1, 70)))
}

# ---------------------------------------------------------------------------
hr("Q10: does the gate run inside bootstrap replicates?")
# Count fs_debias_gate() invocations while a 2-replicate bootstrap runs.
.calls <- new.env(parent = emptyenv()); .calls$n <- 0L
a_boot <- utils::modifyList(mk(sg_focus = "maxSG"),
                            list(debias_gate = TRUE, debias_gate_args = DGA))
fs_boot <- .run_fs_capture(d, a_boot)$value
cat("original analysis: gate ran =", !is.null(fs_boot$debias_gate), "\n")
cat("args_call_all$debias_gate  =", format(fs_boot$args_call_all$debias_gate), "\n")
if (!is.null(fs_boot$sg.harm)) {
  suppressMessages(trace(fs_debias_gate, tracer = function() {
    .calls$n <- .calls$n + 1L }, print = FALSE, where = ns))
  .calls$n <- 0L
  bb <- tryCatch(suppressWarnings(suppressMessages(
    forestsearch_bootstrap_dofuture(fs_boot, nb_boots = 2L,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)))),
    error = function(e) { cat("bootstrap ERROR:", conditionMessage(e), "\n"); NULL })
  suppressMessages(try(untrace(fs_debias_gate, where = ns), silent = TRUE))
  cat("fs_debias_gate() calls during a 2-replicate bootstrap:", .calls$n, "\n")
  cat("  (0 = gate confined to the original analysis;",
      ">0 = it also runs per replicate)\n")
}
cat("\n[stage 2 done]\n")
