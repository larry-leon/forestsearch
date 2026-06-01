# =============================================================================
# Smoke test: GRF subgroup LABEL correctness (regression guard)
# =============================================================================
# Guards the fix for incorrect GRF subgroup definitions (sg.harm.id) at
# grf_depth >= 2.  The adversarial DGM has a true harm region that is a UNION
# of two conjunctions reached by RIGHT-turns -- the case both former label
# builders mishandled.  The test asserts, for survival and GLM:
#
#   (1) membership reconstructed from the PRINTED definition (via get_dfpred)
#       equals GRF's own treat.recommend  -> 0 subjects misclassified;
#   (2) the treatment effect (HR / OR) computed in the printed-label subgroup
#       equals the effect in GRF's actual subgroup.
#
# Run against the INSTALLED package (after devtools::install/load_all):
#   source("grf_label_smoke_test.R")
# Internal arg-builders are reached with forestsearch::: ; exported functions
# (grf.subg.harm.*, get_dfpred) are called directly.
# =============================================================================

suppressMessages({
  library(forestsearch); library(survival); library(data.table)
})

# ---- assertion accumulator -------------------------------------------------
.tc_results <- list()
.tc <- function(desc, pass, detail = "") {
  .tc_results[[length(.tc_results) + 1L]] <<-
    data.frame(check = desc, status = if (isTRUE(pass)) "PASS" else "FAIL",
               detail = detail, stringsAsFactors = FALSE)
  invisible(NULL)
}

# True harm region: (x1 > 0 & x2 > 0) OR (x1 <= 0 & x3 > 0) -- union, right-turns
.in_true_harm <- function(X) (X[, 1] > 0 & X[, 2] > 0) | (X[, 1] <= 0 & X[, 3] > 0)
covs <- paste0("x", 1:5)

# get_dfpred() argument from the structured definition: a brace-label vector
# for a single conjunction, or the disjunctive definition string for a union.
.def_arg <- function(grf_res) {
  d <- grf_res$sg_def
  if (!is.null(d$labels)) d$labels else d$definition
}

eff_in <- function(df, mask, kind) {
  dd <- df[mask, , drop = FALSE]
  if (length(unique(dd$trt)) < 2) return(NA_real_)
  if (kind == "cox") {
    if (sum(dd$trt == 1) < 5 || sum(dd$trt == 0) < 5) return(NA_real_)
    f <- tryCatch(coxph(Surv(time, status) ~ trt, dd), error = function(e) NULL)
  } else {
    f <- tryCatch(glm(Y ~ trt, binomial, dd), error = function(e) NULL)
  }
  if (is.null(f)) NA_real_ else exp(unname(coef(f)[length(coef(f))]))
}

# ---- SURVIVAL --------------------------------------------------------------
set.seed(101); n <- 3000
X <- matrix(runif(n * 5, -1, 1), n, 5); colnames(X) <- covs
W <- rbinom(n, 1, .5); inS <- .in_true_harm(X); tau <- ifelse(inS, 1.3, -0.3)
te <- rexp(n, exp(0.15 * X[, 1] + W * tau)); tc <- runif(n, 0, 1.6)
dS <- data.frame(time = pmin(te, tc), status = as.integer(te <= tc),
                 trt = W, as.data.frame(X), id = seq_len(n))

argsS <- forestsearch:::.build_grf_survival_args(
  data = dS, confounders.name = covs, outcome.name = "time",
  event.name = "status", id.name = "id", treat.name = "trt",
  frac.tau = 0.80, n.min = 60, dmin.grf = 0.0, is.RCT = TRUE,
  grf_depth = 2, seedit = 8316951, return_selected_cuts_only = TRUE)
resS <- do.call(grf.subg.harm.survival, argsS)

grf_mS <- resS$data$treat.recommend == 0
lab_mS <- get_dfpred(dS, .def_arg(resS), version = 2)$treat.recommend == 0
.tc("survival: label membership == GRF membership (0 misclassified)",
    identical(grf_mS, lab_mS),
    sprintf("GRF n=%d, label n=%d, misclassified=%d",
            sum(grf_mS), sum(lab_mS), sum(grf_mS != lab_mS)))
.tc("survival: HR(label) == HR(GRF)",
    isTRUE(all.equal(eff_in(dS, grf_mS, "cox"), eff_in(dS, lab_mS, "cox"))),
    sprintf("HR(GRF)=%.3f HR(label)=%.3f",
            eff_in(dS, grf_mS, "cox"), eff_in(dS, lab_mS, "cox")))

# ---- GLM (binary / OR) -----------------------------------------------------
set.seed(101); n <- 3000
X <- matrix(runif(n * 5, -1, 1), n, 5); colnames(X) <- covs
W <- rbinom(n, 1, .5); inS <- .in_true_harm(X)
Yb <- rbinom(n, 1, plogis(-0.2 + 0.1 * X[, 1] + W * ifelse(inS, 1.8, -0.6)))
dG <- data.frame(Y = Yb, trt = W, as.data.frame(X), id = seq_len(n))

argsG <- forestsearch:::.build_grf_glm_args(
  data = dG, confounders.name = covs, outcome.name = "Y", treat.name = "trt",
  id.name = "id", outcome_type = "binary", n.min = 60, dmin.grf = 0.0,
  is.RCT = TRUE, grf_depth = 2, seedit = 8316951,
  return_selected_cuts_only = TRUE, adverse_outcome = FALSE)
resG <- do.call(grf.subg.harm.glm, argsG)

grf_mG <- resG$data$treat.recommend == 0
lab_mG <- get_dfpred(dG, .def_arg(resG), version = 2)$treat.recommend == 0
.tc("glm: label membership == GRF membership (0 misclassified)",
    identical(grf_mG, lab_mG),
    sprintf("GRF n=%d, label n=%d, misclassified=%d",
            sum(grf_mG), sum(lab_mG), sum(grf_mG != lab_mG)))
.tc("glm: OR(label) == OR(GRF)",
    isTRUE(all.equal(eff_in(dG, grf_mG, "glm"), eff_in(dG, lab_mG, "glm"))),
    sprintf("OR(GRF)=%.3f OR(label)=%.3f",
            eff_in(dG, grf_mG, "glm"), eff_in(dG, lab_mG, "glm")))
.tc("glm: multi-leaf subgroup renders as a disjunction",
    isTRUE(resG$sg_def$is_disjunction),
    resG$sg_def$definition)

# ---- summary ---------------------------------------------------------------
summary_tab <- do.call(rbind, .tc_results)
cat("\n================= GRF label smoke test =================\n")
print(summary_tab, row.names = FALSE)
n_fail <- sum(summary_tab$status == "FAIL")
cat(sprintf("\n%d checks, %d passed, %d failed.\n",
            nrow(summary_tab), sum(summary_tab$status == "PASS"), n_fail))
if (n_fail > 0L) {
  cat("\n** GRF label regression detected. **\n")
} else {
  cat("\nAll GRF label checks passed.\n")
}
