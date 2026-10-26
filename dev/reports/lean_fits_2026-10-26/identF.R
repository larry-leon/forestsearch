# identF.R -- Part-F identity pass (binary OR search fit through glm.fit()).
# Task: dev/tasks/cc_task_lean_fits_2026-10-26.md, section 4.
# Run from this directory against the installed A+F build: Rscript identF.R
source("../vi_deprecation_2026-10-26/vi_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Built")]
print(pd)

SUMMARY <- file.path(OUT, "identF_summary.txt")
cat(sprintf("Part-F identity pass -- installed forestsearch %s (Built %s)\n", pd$Version, pd$Built),
    file = SUMMARY)
sline <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = SUMMARY, append = TRUE, sep = "") }
failed <- character(0)
gate <- function(ok, label) {
  sline(if (isTRUE(ok)) "PASS" else "FAIL", ": ", label)
  if (!isTRUE(ok)) failed <<- c(failed, label)
}
with_warnings <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}
kind <- function(v) if (is.null(v)) "NULL" else if (isTRUE(v$nonestimable)) "nonestimable" else "fit"

fx <- FIXTURES$binary()
args <- variant_args(fx$args, "ii")
treat.name <- args$treat.name; outcome.name <- args$outcome.name
adverse_outcome <- args$adverse_outcome

# ---- per candidate: capture via trace() --------------------------------------
cap <- new.env(); cap$k <- 0L; cap$idx <- list(); cap$slow <- 0L
fast_tracer <- bquote({
  .k <- get("k", envir = .(cap)) + 1L
  assign("k", .k, envir = .(cap))
  .l <- get("idx", envir = .(cap)); .l[[.k]] <- id.x; assign("idx", .l, envir = .(cap))
  if (.k == 1L) {
    assign("yy", yy, envir = .(cap)); assign("tt", tt, envir = .(cap))
    assign("estimator_fn", estimator_fn, envir = .(cap))
  }
})
slow_tracer <- bquote(assign("slow", get("slow", envir = .(cap)) + 1L, envir = .(cap)))
ns <- asNamespace("forestsearch")
suppressMessages({
  trace("fit_glm_for_subgroup_fast", tracer = fast_tracer, where = ns, print = FALSE)
  trace("fit_glm_for_subgroup", tracer = slow_tracer, where = ns, print = FALSE)
})
rc_trace <- run_collect(args)
suppressMessages({
  untrace("fit_glm_for_subgroup_fast", where = ns)
  untrace("fit_glm_for_subgroup", where = ns)
})
sline(sprintf("INFO: captured %d fit_glm_for_subgroup_fast() calls (Stage-0b count 2577); slow fit_glm_for_subgroup() calls = %d",
              cap$k, cap$slow))
# Run 1 of this script gated on the task's expected 2577 (Stage-0b) and failed
# with 2584 (identF_log_run1_count_premise.txt): 2577 is n_passed_cox of the
# VI-screened (vi.grf.min = -0.2) binary run on CRAN 0.4.0 (committed
# stageA_binary_i.rds); the "ii" family has n_passed_cox = 2584 on both 0.4.0
# and the deprecation build.  The count is anchored to that committed record.
n_expect <- readRDS("../vi_deprecation_2026-10-26/stageC_binary_ii.rds")$filter_counts$n_passed_cox
sline(sprintf("INFO: Stage-0b count 2577 is the VI-screened (i) family (stageA_binary_i.rds n_passed_cox = %d); the ii family's committed n_passed_cox = %d",
              readRDS("../vi_deprecation_2026-10-26/stageA_binary_i.rds")$filter_counts$n_passed_cox, n_expect))
gate(cap$k == n_expect, sprintf("captured candidate count equals the committed ii record's n_passed_cox (%d)", n_expect))
gate(cap$slow == 0L, "routing: zero calls to the slow fit_glm_for_subgroup() during the run")

estimator_fn <- cap$estimator_fn
gate(isTRUE(attr(estimator_fn, "fs_fast_unadjusted")), "routing: attr(estimator_fn, 'fs_fast_unadjusted') is TRUE")
gate(identical(environment(estimator_fn)$fs_fast_core, "logit_or"), "routing: environment(estimator_fn)$fs_fast_core is 'logit_or'")

yy <- cap$yy; tt <- cap$tt
df_ref <- data.frame(yy, tt); names(df_ref) <- c(outcome.name, treat.name)
slow_fn <- forestsearch:::.make_glm_binary_estimator(treat.name, outcome.name, effect_measure = "OR",
                                                     adverse_outcome = adverse_outcome)
fast_entry <- forestsearch:::fit_glm_for_subgroup_fast
slow_entry <- forestsearch:::fit_glm_for_subgroup

n <- cap$k
same_val <- same_warn <- logical(n); kinds_fast <- kinds_slow <- character(n)
n_warn_cand <- 0L
clo_checked <- 0L; clo_same <- 0L
for (i in seq_len(n)) {
  id.x <- cap$idx[[i]]
  a <- with_warnings(fast_entry(yy, tt, id.x, estimator_fn))
  b <- with_warnings(slow_entry(df_ref, id.x, slow_fn))
  same_val[i]  <- identical(a$value, b$value)
  same_warn[i] <- identical(a$warnings, b$warnings)
  kinds_fast[i] <- kind(a$value); kinds_slow[i] <- kind(b$value)
  if (length(b$warnings)) n_warn_cand <- n_warn_cand + 1L
  if (i %% 50L == 0L || kinds_slow[i] != "fit") {
    sl <- df_ref[id.x == 1, , drop = FALSE]
    ca <- with_warnings(estimator_fn(sl)); cb <- with_warnings(slow_fn(sl))
    clo_checked <- clo_checked + 1L
    if (identical(ca$value, cb$value) && identical(ca$warnings, cb$warnings)) clo_same <- clo_same + 1L
  }
}
kt <- table(fast = kinds_fast, slow = kinds_slow)
sline("INFO: kinds table (fast entry x slow entry): ",
      paste(sprintf("%s/%s = %d", rep(rownames(kt), ncol(kt)), rep(colnames(kt), each = nrow(kt)), as.vector(kt)),
            collapse = "; "))
sline(sprintf("INFO: per candidate: %d/%d identical in value, %d/%d identical warning sets; candidates with warnings (slow) = %d",
              sum(same_val), n, sum(same_warn), n, n_warn_cand))
gate(all(same_val), "per candidate (entry level): all values identical to the slow chain")
gate(all(same_warn), "per candidate (entry level): all warning sets identical")
gate(identical(kinds_fast, kinds_slow), "per candidate: kinds (NULL / nonestimable / fit) equal")
sline(sprintf("INFO: closure level: %d/%d checked identical (every 50th + every non-fit candidate)", clo_same, clo_checked))
gate(clo_same == clo_checked && clo_checked > 0L, "closure level: full returned list and warnings identical")

# ---- end to end vs the committed deprecation record (all three fixtures) ----
for (fx_name in c("binary", "continuous", "survival")) {
  f <- FIXTURES[[fx_name]]()
  rc <- run_collect(variant_args(f$args, "ii"))   # untraced run
  x <- extract_run(rc, f$args$id.name)
  p <- file.path(OUT, sprintf("identF_%s_ii.rds", fx_name))
  saveRDS(x, p)
  xs <- readRDS(p)
  xr <- readRDS(sprintf("../vi_deprecation_2026-10-26/stageC_%s_ii.rds", fx_name))
  for (k in c(IDENTITY_COMPONENTS, "n_fits", "n_warnings", "warnings"))
    gate(identical(xs[[k]], xr[[k]]), sprintf("end to end %s ii: %s identical to committed stageC_%s_ii.rds", fx_name, k, fx_name))
  sline("INFO: [report only] ", fx_name, " ii args_call_all identical = ", identical(xs$args_call_all, xr$args_call_all))
  sline("INFO: ", fx_name, " ii sg.harm = ", paste(xs$sg.harm, collapse = " & "), "; warnings = ", xs$n_warnings)
}

sline(if (length(failed)) sprintf("PART-F IDENTITY GATES: %d FAILED", length(failed))
      else "PART-F IDENTITY GATES: ALL PASS")
if (length(failed)) quit(status = 1)
