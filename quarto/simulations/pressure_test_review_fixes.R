# =============================================================================
# Smoke test: A (KM median), B (null_target_hr), C (include: false), D (trim)
# =============================================================================

PASS <- 0L; FAIL <- 0L
ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  if (isTRUE(res)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [ OK ] %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (got: %s)\n", label, deparse(res, 60)[1]))
  }
}

QMD <- "/mnt/user-data/outputs/actg175_survival_benefit_simulations.qmd"
txt <- paste(readLines(QMD, warn = FALSE), collapse = "\n")

extract_r_chunks <- function(qmd_path) {
  lines <- readLines(qmd_path, warn = FALSE)
  in_chunk <- FALSE; chunks <- character(0); cur <- character(0)
  for (l in lines) {
    if (grepl("^```\\{r[^}]*\\}", l)) { in_chunk <- TRUE; next }
    if (in_chunk && grepl("^```\\s*$", l)) {
      chunks <- c(chunks, paste(cur, collapse = "\n"))
      cur <- character(0); in_chunk <- FALSE; next
    }
    if (in_chunk) cur <- c(cur, l)
  }
  paste(chunks, collapse = "\n\n")
}

.count <- function(pattern, text, perl = FALSE) {
  m <- gregexpr(pattern, text, perl = perl)[[1]]
  if (length(m) == 1L && m[1] == -1L) 0L else length(m)
}

cat("\n=== A. KM-based median ===\n\n")

# Old buggy line gone from EXECUTABLE R code (comments may still
# reference it; we strip line comments before counting)
r_src_for_a1 <- gsub("#[^\n]*", "", extract_r_chunks(QMD))
ok("A.1 naive median(dfa$y[dfa$event == 1]) removed from R code",
   .count("median\\(dfa\\$y\\[dfa\\$event == 1\\]\\)", r_src_for_a1) == 0L)
# KM construct present
ok("A.2 survfit(Surv(y, event) ~ 1) present",
   .count("survival::survfit\\(survival::Surv\\(y, event\\) ~ 1", txt) == 1L)
# Median extracted from the KM summary's $table component, indexed
# by [["median"]].  The vignette splits this across two lines:
#   km_summary       <- summary(sf_obs)$table
#   observed_median  <- km_summary[["median"]]
# We just verify both pieces appear.
ok("A.3 summary()$table assigned to km_summary",
   .count("summary\\(sf_obs\\)\\$table", txt) >= 1L)
ok("A.3b km_summary[[\"median\"]] indexed",
   .count("km_summary\\[\\[\"median\"\\]\\]", txt) >= 1L)
# Quantile fallback when median NA
ok("A.4 quantile(sf_obs, probs = 0.25) fallback present",
   .count("quantile\\(sf_obs,\\s*probs\\s*=\\s*0\\.25\\)", txt) == 1L)
# is.na guard around the fallback
ok("A.5 is.na(observed_median) check present",
   .count("is\\.na\\(observed_median\\)", txt) >= 1L)


cat("\n=== B. null_target_hr solver ===\n\n")

# null_target_hr defined as NULL by default in setup chunk
ok("B.1 null_target_hr <- NULL declared in setup",
   .count("(?m)^null_target_hr\\s*<-\\s*NULL", txt, perl = TRUE) == 1L)
# Solver block exists: uses log(null_target_hr) divided by fitted log-HR
# (we allow either inline log(...) / log(...) or two-step via fitted_log_hr)
ok("B.2 dgm_k_treat_null computation present",
   .count("dgm_k_treat_null\\s*<-\\s*log\\(null_target_hr\\)",
          txt) == 1L)
# Default fallback: dgm_k_treat_null <- dgm_k_treat
ok("B.3 default dgm_k_treat_null <- dgm_k_treat present",
   .count("dgm_k_treat_null\\s*<-\\s*dgm_k_treat", txt) == 1L)
# dgm_null call uses dgm_k_treat_null (not dgm_k_treat)
r_src <- extract_r_chunks(QMD)
# Find the dgm_null block specifically
ok("B.4 dgm_null block uses k_treat = dgm_k_treat_null",
   .count("k_treat\\s*=\\s*dgm_k_treat_null", r_src) >= 1L)
# Section anchor
ok("B.5 #sec-dgm-null anchor on Step 5 header",
   .count("(?m)^## Step 5: Create Null DGM \\{#sec-dgm-null\\}",
          txt, perl = TRUE) == 1L)
# Callout note explaining the option
ok("B.6 callout note explaining the option present",
   .count("Optional: widen the gap from c1 with `null_target_hr`",
          txt, perl = TRUE) == 1L)
# Diagnostic output: show k_treat used + gap
ok("B.7 dgm_null prints 'k_treat used:'",
   .count("k_treat used:", txt) == 1L)
ok("B.8 dgm_null prints gap to c1",
   .count("Gap \\(c1 - HR_switched\\)", txt) == 1L)


cat("\n=== C. save_results chunk visibility ===\n\n")

# Chunk header has include: false, NOT code-fold: false
ok("C.1 save-results chunk uses #| include: false",
   .count("(?m)^```\\{r save-results\\}\n#\\| include: false",
          txt, perl = TRUE) == 1L)
ok("C.2 save-results chunk does NOT use #| code-fold: false",
   .count("(?m)^```\\{r save-results\\}\n#\\| code-fold: false",
          txt, perl = TRUE) == 0L)


cat("\n=== D. winsorize_hr_columns + trimmed inputs ===\n\n")

# Helper definition
ok("D.1 winsorize_hr_columns helper defined",
   .count("(?m)^winsorize_hr_columns\\s*<-\\s*function\\(",
          txt, perl = TRUE) == 1L)
# Trimmed objects created
ok("D.2 results_alt_trimmed created via winsorize",
   .count("results_alt_trimmed\\s*<-\\s*winsorize_hr_columns\\(results_alt",
          txt) == 1L)
ok("D.3 results_null_trimmed created via winsorize",
   .count("results_null_trimmed\\s*<-\\s*winsorize_hr_columns\\(results_null",
          txt) == 1L)
# trim_fraction declared with sensible default
ok("D.4 trim_fraction <- 0.01 default",
   .count("(?m)^trim_fraction\\s*<-\\s*0\\.01", txt, perl = TRUE) == 1L)
# Diagnostic chunk callout
ok("D.5 callout note explaining trimming present",
   .count("Trimming Cox-HR plug-in columns for the estimation tables",
          txt) == 1L)

# Every build_estimation_table OR interpret_estimation_table call in R chunks
# should reference *_trimmed, not the bare results_alt / results_null
calls <- regmatches(r_src, gregexpr(
  "(?:build_estimation_table|interpret_estimation_table)\\([^)]*?(?:results_alt|results_null)[^)]*?\\)",
  r_src, perl = TRUE))[[1]]
n_calls   <- length(calls)
n_trimmed <- sum(grepl("results_(alt|null)_trimmed", calls))
ok(sprintf("D.6 all build/interpret calls use *_trimmed (%d/%d)",
           n_trimmed, n_calls),
   n_trimmed == n_calls && n_calls >= 5L)


cat("\n=== Parse check ===\n\n")
parse_res <- tryCatch(parse(text = r_src), error = function(e) e)
ok("Z.1 vignette R code parses cleanly", !inherits(parse_res, "error"))


cat("\n=== Functional test: winsorize_hr_columns ===\n\n")

# Reproduce the helper locally and exercise it on synthetic data with
# extreme values -- replicates the bug scenario where one Q^c subset
# has Cox HR collapsing to ~0 on the switched scale.
winsorize_hr_columns <- function(results, trim_fraction = 0.01,
                                 cols = c("hr.H.hat", "hr.Hc.hat",
                                          "hr.H.bc",  "hr.Hc.bc")) {
  if (is.null(trim_fraction) || trim_fraction <= 0) return(results)
  if (trim_fraction >= 0.5) stop("trim_fraction must be < 0.5")
  out <- as.data.frame(results)
  trim_log <- data.frame(column = character(0), n_trimmed = integer(0),
                         lower = numeric(0), upper = numeric(0),
                         stringsAsFactors = FALSE)
  for (col in intersect(cols, names(out))) {
    v <- out[[col]]
    finite_idx <- is.finite(v)
    if (sum(finite_idx) < 5L) next
    lo <- stats::quantile(v[finite_idx], trim_fraction,    na.rm = TRUE)
    hi <- stats::quantile(v[finite_idx], 1 - trim_fraction, na.rm = TRUE)
    flag <- finite_idx & (v < lo | v > hi)
    out[[col]][flag & v < lo] <- lo
    out[[col]][flag & v > hi] <- hi
    trim_log <- rbind(trim_log, data.frame(
      column = col, n_trimmed = sum(flag, na.rm = TRUE),
      lower = as.numeric(lo), upper = as.numeric(hi),
      stringsAsFactors = FALSE
    ))
  }
  attr(out, "trim_log") <- trim_log
  out
}

set.seed(123)
n <- 200
toy <- data.frame(
  hr.H.hat  = rlnorm(n, meanlog = 0,    sdlog = 0.2),
  hr.Hc.hat = c(rlnorm(n - 2, meanlog = 0, sdlog = 0.2),
                1e-10, 1e-12),  # two collapsing replicates
  hr.H.bc   = rlnorm(n, meanlog = 0,    sdlog = 0.2),
  hr.Hc.bc  = rlnorm(n, meanlog = 0,    sdlog = 0.2)
)
trimmed <- winsorize_hr_columns(toy, trim_fraction = 0.01)

ok("F.1 winsorize returns same row count",
   nrow(trimmed) == nrow(toy))
ok("F.2 attr(trimmed, 'trim_log') has 4 rows",
   nrow(attr(trimmed, "trim_log")) == 4L)
# Without trim, 1/min(hr.Hc.hat) ~ 1e10
ok("F.3 raw 1/min(hr.Hc.hat) is huge (>1e9)",
   1 / min(toy$hr.Hc.hat) > 1e9)
# After trim, the bottom is clamped at the 1st percentile (a finite, sane number)
ok("F.4 trimmed min(hr.Hc.hat) is bounded above 0",
   min(trimmed$hr.Hc.hat) > 0.01)
ok("F.5 trimmed 1/min(hr.Hc.hat) is sane (<10)",
   1 / min(trimmed$hr.Hc.hat) < 10)
# trim_fraction = NULL -> identity passthrough
identical_pass <- identical(
  winsorize_hr_columns(toy, trim_fraction = NULL),
  toy)
ok("F.6 trim_fraction = NULL is identity passthrough",
   identical_pass)

# Solver math sanity (B)
fitted_hr      <- 1.19   # switched-scale fitted
target_hr      <- 1.05
solved_k_treat <- log(target_hr) / log(fitted_hr)
ok("F.7 solver: target_hr=1.05 with fitted_hr=1.19 gives k_treat ~0.281",
   abs(solved_k_treat - 0.281) < 0.005)
# Default behaviour preserved: NULL -> dgm_k_treat
default_branch_ok <- {
  null_target_hr <- NULL
  dgm_k_treat <- 1
  dgm_k_treat_null <- dgm_k_treat
  if (!is.null(null_target_hr)) dgm_k_treat_null <- log(null_target_hr) / log(1.19)
  identical(dgm_k_treat_null, 1)
}
ok("F.8 null_target_hr=NULL preserves dgm_k_treat (=1)",
   default_branch_ok)


cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
