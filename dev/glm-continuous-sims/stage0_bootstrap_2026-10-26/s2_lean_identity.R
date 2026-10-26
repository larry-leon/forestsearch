# s2_lean_identity.R -- candidate A sized and its bit-identity tested, with
# NO change to the package:
#
#  1. base::trace() on forestsearch:::fit_cox_for_subgroup (in this session
#     only) records the subset vectors (Y, E, Treat) of every candidate that
#     reaches the search-stage Cox fit during ONE forestsearch() call on
#     fixture V.  untrace() follows immediately.
#  2. base::trace() on survival::coxph.fit records the arguments coxph()
#     passes to it during ONE formula fit, so the lean call below uses the
#     installed survival's conventions rather than assumptions.  untrace().
#  3. For every captured candidate: the package's formula block
#     (R/subgroup_search.R lines 834-845, transplanted verbatim) and a lean
#     block that calls survival::coxph.fit() directly.  identical() on the
#     consumed CI row; agreement of the failure sets; timings (5 reps,
#     alternating order).
#
# Output: the machinery seconds the lean path would remove per call, the
# identity count, and the captured coxph.fit() conventions (the recipe).

source("stage0_common.R")
fx <- stage0_fixture()

# ---- 1. capture the search-stage Cox-fit inputs ------------------------------
.cap <- new.env()
.cap$n <- 0L
.cap$items <- vector("list", 20000L)
.cap$push <- function(item) {
  .cap$n <- .cap$n + 1L
  if (.cap$n > length(.cap$items)) length(.cap$items) <- 2L * length(.cap$items)
  .cap$items[[.cap$n]] <- item
}
assign(".cap", .cap, envir = globalenv())

stopifnot(exists("fit_cox_for_subgroup", envir = asNamespace("forestsearch"), inherits = FALSE))
invisible(capture.output(
  trace("fit_cox_for_subgroup", where = asNamespace("forestsearch"), print = FALSE,
        tracer = quote({
          .GlobalEnv$.cap$push(list(
            Y = yy[id.x == 1], E = dd[id.x == 1], Treat = tt[id.x == 1],
            adj = forestsearch:::.fs_adjust_terms(adjust_covariates)))
        }))
))
on.exit(try(suppressMessages(untrace("fit_cox_for_subgroup", where = asNamespace("forestsearch"))), silent = TRUE), add = TRUE)
t_traced <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
invisible(capture.output(untrace("fit_cox_for_subgroup", where = asNamespace("forestsearch"))))

K <- .cap$n
cands <- .cap$items[seq_len(K)]
cc <- candidate_counts(fs)
kv(capture_n_fits = K,
   capture_matches_filter_counts = identical(K, as.integer(cc$n_reached_cox_fit)),
   capture_n_adjusted = sum(vapply(cands, function(z) length(z$adj) > 0L, logical(1))),
   capture_median_n = median(vapply(cands, function(z) length(z$Y), numeric(1))),
   capture_call_s_under_trace = t_traced)
if (K == 0L) {
  cat("GATE FAILED: no fit_cox_for_subgroup() call was captured -- the search did not run",
      "in this process (future fallback changed?).  Nothing to measure.\n")
  quit(status = 4)
}
if (any(vapply(cands, function(z) length(z$adj) > 0L, logical(1)))) {
  cat("NOTE: adjusted candidates present; the lean path applies to the unadjusted branch only.",
      "They are excluded from the comparison below and counted above.\n")
  cands <- Filter(function(z) length(z$adj) == 0L, cands)
}
cands <- Filter(function(z) length(z$Y) > 0L, cands)
K <- length(cands)

# ---- 2. capture what coxph() passes to coxph.fit() -----------------------------
.capfit <- new.env(); assign(".capfit", .capfit, envir = globalenv())
# Every formal of the installed coxph.fit() is captured by name, so a formal
# this script does not know about (a newer survival) is still recorded and
# reused below rather than silently left at its default.
invisible(capture.output(
  trace("coxph.fit", where = asNamespace("survival"), print = FALSE,
        tracer = quote({
          .fr <- environment()
          .nm <- setdiff(names(formals(survival::coxph.fit)), "...")
          .GlobalEnv$.capfit$args <- lapply(setNames(.nm, .nm), function(.a)
            tryCatch(get(.a, envir = .fr, inherits = FALSE),
                     error = function(e) structure(list(), class = "stage0_capture_failed")))
        }))
))
probe <- cands[[which.max(vapply(cands, function(z) length(z$Y), numeric(1)))]]
Y <- probe$Y; E <- probe$E; Treat <- probe$Treat
fit_probe <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
invisible(capture.output(untrace("coxph.fit", where = asNamespace("survival"))))
a <- .capfit$args
stopifnot(!is.null(a), all(c("x", "y") %in% names(a)))
failed_caps <- names(a)[vapply(a, inherits, logical(1), "stage0_capture_failed")]
a <- a[setdiff(names(a), failed_caps)]
kv(conv_formals_not_capturable = if (length(failed_caps)) paste(failed_caps, collapse = ",") else "<none>")
ctl <- survival::coxph.control()
fit_formals <- names(formals(survival::coxph.fit))
y_raw <- survival::Surv(Y, E)
y_fix <- if (isTRUE(ctl$timefix)) survival::aeqSurv(y_raw) else y_raw
same_vals <- function(u, v) isTRUE(all.equal(unclass(unname(u)), unclass(unname(v)), check.attributes = FALSE))
conv <- list(
  coxph_fit_formals   = paste(fit_formals, collapse = ","),
  x_storage           = storage.mode(a$x),
  x_colnames          = paste(colnames(a$x), collapse = ","),
  x_values_equal_treat = same_vals(a$x[, 1L], as.double(Treat)),
  y_type              = attr(a$y, "type"),
  y_values_equal_raw  = same_vals(a$y, y_raw),
  y_values_equal_aeqSurv = same_vals(a$y, y_fix),
  timefix             = isTRUE(ctl$timefix),
  strata_is_null      = is.null(a$strata),
  offset_desc         = if (is.null(a$offset)) "NULL" else sprintf("numeric length %d, all zero: %s", length(a$offset), all(a$offset == 0)),
  init_desc           = if (is.null(a$init)) "NULL" else paste(format(a$init), collapse = ","),
  control_is_default  = identical(a$control, ctl),
  weights_is_null     = is.null(a$weights),
  method              = a$method,
  rownames_desc       = if (is.null(a$rownames)) "NULL" else sprintf("%s length %d", class(a$rownames)[1], length(a$rownames)),
  resid               = if (!"resid" %in% names(a)) "<no formal>" else as.character(a$resid),
  nocenter            = if (!"nocenter" %in% names(a)) "<no formal>" else if (is.null(a$nocenter)) "NULL" else paste(a$nocenter, collapse = ","),
  other_formals       = {
    oth <- setdiff(names(a), c("x", "y", "strata", "offset", "init", "control", "weights",
                               "method", "rownames", "resid", "nocenter"))
    if (length(oth)) paste(vapply(oth, function(k) paste0(k, "=", paste(deparse(a[[k]]), collapse = "")), character(1)), collapse = "; ") else "<none>"
  }
)
do.call(kv, setNames(conv, paste0("conv_", names(conv))))

# ---- 3. the two blocks ---------------------------------------------------------
# (a) the package's own block, R/subgroup_search.R:834-845 verbatim (CRAN 0.4.0)
formula_block <- function(Y, E, Treat) {
  try({
    fit.cox <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
    beta <- fit.cox$coefficients
    se_beta <- sqrt(diag(fit.cox$var))
    z_ci <- stats::qnorm((1 + 0.95) / 2, 0, 1)
    ci <- cbind(exp(beta), exp(-beta),
                exp(beta - z_ci * se_beta), exp(beta + z_ci * se_beta))
    dimnames(ci) <- list(names(beta),
                         c("exp(coef)", "exp(-coef)",
                           "lower .95", "upper .95"))
    ci
  }, silent = TRUE)
}
# (b) the lean block: coxph.fit() with the conventions captured in step 2.
#     Data-dependent arguments (x, y, offset, rownames) are rebuilt per
#     candidate; every other captured formal (strata, init, control, weights,
#     method, resid, nocenter, and anything a newer survival adds) is reused
#     exactly as captured.  names(beta) <- "Treat" restores the coefficient
#     name coxph() takes from the model-matrix column, so the CI row's
#     dimnames agree.
lean_args_fixed <- a[setdiff(names(a), c("x", "y", "offset", "rownames"))]
use_offset_zero <- !is.null(a$offset)
use_rownames    <- !is.null(a$rownames)
use_aeqSurv     <- isTRUE(ctl$timefix)
lean_block <- function(Y, E, Treat) {
  try({
    n <- length(Y)
    ys <- survival::Surv(Y, E)
    if (use_aeqSurv) ys <- survival::aeqSurv(ys)
    args <- c(list(x = matrix(as.double(Treat), ncol = 1L),
                   y = ys,
                   offset = if (use_offset_zero) rep(0, n) else NULL,
                   rownames = if (use_rownames) as.character(seq_len(n)) else NULL),
              lean_args_fixed)
    fit.cox <- do.call(survival::coxph.fit, args)
    beta <- fit.cox$coefficients
    names(beta) <- "Treat"
    se_beta <- sqrt(diag(fit.cox$var))
    z_ci <- stats::qnorm((1 + 0.95) / 2, 0, 1)
    ci <- cbind(exp(beta), exp(-beta),
                exp(beta - z_ci * se_beta), exp(beta + z_ci * se_beta))
    dimnames(ci) <- list(names(beta),
                         c("exp(coef)", "exp(-coef)",
                           "lower .95", "upper .95"))
    ci
  }, silent = TRUE)
}

run_all <- function(block) {
  res <- vector("list", K); nwarn <- integer(K)
  for (k in seq_len(K)) {
    z <- cands[[k]]
    res[[k]] <- withCallingHandlers(block(z$Y, z$E, z$Treat),
                 warning = function(w) { nwarn[k] <<- nwarn[k] + 1L; invokeRestart("muffleWarning") })
  }
  list(res = res, nwarn = nwarn)
}

# identity pass (unordered timing -- separate from the timed reps below)
ra <- run_all(formula_block); rb <- run_all(lean_block)
err_a <- vapply(ra$res, inherits, logical(1), "try-error")
err_b <- vapply(rb$res, inherits, logical(1), "try-error")
ident <- mapply(function(u, v) identical(u, v), ra$res, rb$res)
ident_vals <- mapply(function(u, v) {
  if (inherits(u, "try-error") || inherits(v, "try-error")) return(NA)
  identical(unname(u), unname(v)) }, ra$res, rb$res)
kv(ident_n_candidates = K,
   ident_ci_row_identical = sum(ident),
   ident_ci_values_identical = sum(ident_vals, na.rm = TRUE),
   ident_errors_formula = sum(err_a), ident_errors_lean = sum(err_b),
   ident_error_sets_agree = identical(err_a, err_b),
   ident_warn_formula_total = sum(ra$nwarn), ident_warn_lean_total = sum(rb$nwarn),
   ident_warn_sets_agree = identical(ra$nwarn > 0L, rb$nwarn > 0L),
   ident_identical_among_warned = if (any(ra$nwarn > 0L)) sum(ident[ra$nwarn > 0L]) else NA,
   ident_n_warned = sum(ra$nwarn > 0L))
if (any(!ident)) {
  k <- which(!ident)[1L]
  cat("\nFirst non-identical candidate (k =", k, "): n =", length(cands[[k]]$Y), "\n")
  cat("formula:\n"); print(ra$res[[k]]); cat("lean:\n"); print(rb$res[[k]])
  if (!inherits(ra$res[[k]], "try-error") && !inherits(rb$res[[k]], "try-error"))
    cat("max abs diff:", max(abs(ra$res[[k]] - rb$res[[k]])), "\n")
}

# ---- timings: 5 reps, alternating order --------------------------------------
time_block <- function(block) {
  system.time(for (k in seq_len(K)) { z <- cands[[k]]; suppressWarnings(block(z$Y, z$E, z$Treat)) })[["elapsed"]]
}
ta <- tb <- numeric(5)
for (r in 1:5) {
  if (r %% 2L == 1L) { ta[r] <- time_block(formula_block); tb[r] <- time_block(lean_block) }
  else               { tb[r] <- time_block(lean_block);    ta[r] <- time_block(formula_block) }
}
kv(time_formula_s_reps = ta, time_lean_s_reps = tb,
   time_formula_s_median = median(ta), time_lean_s_median = median(tb),
   time_ms_per_fit_formula = 1000 * median(ta) / K,
   time_ms_per_fit_lean = 1000 * median(tb) / K,
   machinery_s_per_call = median(ta) - median(tb),
   machinery_share_of_formula = (median(ta) - median(tb)) / median(ta))

save_result(list(K = K, conv = conv, capfit = a, ident = ident, ident_vals = ident_vals,
                 err_a = err_a, err_b = err_b, nwarn_a = ra$nwarn, nwarn_b = rb$nwarn,
                 ta = ta, tb = tb, cand_sizes = vapply(cands, function(z) length(z$Y), numeric(1))),
            "s2_lean_identity")
cat("\ns2 complete\n")
