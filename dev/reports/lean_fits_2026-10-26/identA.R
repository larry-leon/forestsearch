# identA.R -- Part-A identity pass (survival search fit through coxph.fit()).
# Task: dev/tasks/cc_task_lean_fits_2026-10-26.md, section 3.
# Run from this directory against the installed A build: Rscript identA.R
source("../vi_deprecation_2026-10-26/vi_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Built")]
print(pd)

SUMMARY <- file.path(OUT, "identA_summary.txt")
cat(sprintf("Part-A identity pass -- installed forestsearch %s (Built %s)\n", pd$Version, pd$Built),
    file = SUMMARY)
sline <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = SUMMARY, append = TRUE, sep = "") }
failed <- character(0)
gate <- function(ok, label) {
  sline(if (isTRUE(ok)) "PASS" else "FAIL", ": ", label)
  if (!isTRUE(ok)) failed <<- c(failed, label)
}

# CRAN 0.4.0 R/subgroup_search.R L834-L845, L873, L875-L881, L888-L894, verbatim
ref_fit_cox <- function(Y, E, Treat) {
  hr.cox <- try({
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
  if (inherits(hr.cox, "try-error")) return(NULL)
  if (is.null(dim(hr.cox))) {
    hr.cox <- matrix(hr.cox, nrow = 1L, dimnames = list("Treat", names(hr.cox)))
  }
  trow <- if ("Treat" %in% rownames(hr.cox)) hr.cox["Treat", ] else hr.cox[1, ]
  list(hr = trow[1], lower = trow[3], upper = trow[4], med0 = NA_real_, med1 = NA_real_)
}

with_warnings <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(value = v, warnings = w)
}

fx <- FIXTURES$survival()
args <- variant_args(fx$args, "ii")

# ---- per candidate: capture inputs via trace() ------------------------------
cap <- new.env(); cap$k <- 0L; cap$x <- list()
tracer <- bquote({
  .k <- get("k", envir = .(cap)) + 1L
  assign("k", .k, envir = .(cap))
  .x <- get("x", envir = .(cap))
  .x[[.k]] <- list(Y = yy[id.x == 1], E = dd[id.x == 1], Treat = tt[id.x == 1],
                   adj = forestsearch:::.fs_adjust_terms(adjust_covariates))
  assign("x", .x, envir = .(cap))
})
ns <- asNamespace("forestsearch")
suppressMessages(trace("fit_cox_for_subgroup", tracer = tracer, where = ns, print = FALSE))
rc_trace <- run_collect(args)
suppressMessages(untrace("fit_cox_for_subgroup", where = ns))
cands <- cap$x
n_adj <- sum(vapply(cands, function(z) length(z$adj) > 0L, logical(1)))
sline(sprintf("INFO: captured %d fit_cox_for_subgroup() calls (Stage-0 count 1410); adjusted calls = %d",
              length(cands), n_adj))
gate(length(cands) == 1410L, "captured candidate count equals the Stage-0 count 1410")
gate(n_adj == 0L, "all captured calls unadjusted")

fit_cox <- forestsearch:::fit_cox_for_subgroup
same_val <- same_warn <- logical(length(cands))
n_warn_cand <- 0L; n_null <- 0L
for (i in seq_along(cands)) {
  z <- cands[[i]]
  a <- with_warnings(fit_cox(z$Y, z$E, z$Treat, rep(1L, length(z$Y))))
  b <- with_warnings(ref_fit_cox(z$Y, z$E, z$Treat))
  same_val[i]  <- identical(a$value, b$value)
  same_warn[i] <- identical(a$warnings, b$warnings)
  if (length(b$warnings)) n_warn_cand <- n_warn_cand + 1L
  if (is.null(b$value)) n_null <- n_null + 1L
}
sline(sprintf("INFO: per candidate: %d/%d identical in value, %d/%d identical warning sets; candidates with warnings (reference) = %d; NULL returns = %d",
              sum(same_val), length(cands), sum(same_warn), length(cands), n_warn_cand, n_null))
gate(all(same_val), "per candidate: all values identical to the 0.4.0 formula block")
gate(all(same_warn), "per candidate: all warning sets identical")

# ---- end to end vs the committed deprecation record -------------------------
rc <- run_collect(args)   # untraced run
x <- extract_run(rc, fx$args$id.name)
saveRDS(x, file.path(OUT, "identA_survival_ii.rds"))
xs <- readRDS(file.path(OUT, "identA_survival_ii.rds"))
xr <- readRDS("../vi_deprecation_2026-10-26/stageC_survival_ii.rds")
for (k in c(IDENTITY_COMPONENTS, "n_fits", "n_warnings", "warnings"))
  gate(identical(xs[[k]], xr[[k]]), sprintf("end to end survival ii: %s identical to committed stageC_survival_ii.rds", k))
sline("INFO: [report only] args_call_all identical = ", identical(xs$args_call_all, xr$args_call_all))
sline("INFO: survival ii sg.harm = ", paste(xs$sg.harm, collapse = " & "), "; warnings = ", xs$n_warnings)

sline(if (length(failed)) sprintf("PART-A IDENTITY GATES: %d FAILED", length(failed))
      else "PART-A IDENTITY GATES: ALL PASS")
if (length(failed)) quit(status = 1)
