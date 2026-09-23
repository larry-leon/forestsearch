# Usage: Rscript postcond_fits.R <out.rds>
# Fixed-seed fits on every path the declaration-calibration change can reach,
# with timing fields stripped; digests via base serialize + tools::md5sum.
args <- commandArgs(TRUE); out <- args[1]
suppressMessages(devtools::load_all(quiet = TRUE))
strip <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("timing_seconds", "minutes_all", "time_search")
    if (any(drop)) x <- x[!drop]
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- strip(x[[i]])
  }
  x
}
md5 <- function(x) { f <- tempfile(); writeBin(serialize(x, NULL, version = 3), f)
  unname(tools::md5sum(f)) }
gb <- survival::gbsg; gb$id <- seq_len(nrow(gb)); gb$time_months <- gb$rfstime / 30.4375
gb$grade3 <- ifelse(gb$grade == "3", 1, 0)
gargs <- list(df.analysis = gb, outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", maxk = 2,
  hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  n.min = 60, d0.min = 12, d1.min = 12, use_twostage = FALSE,
  seedit = 8316951, quiet = TRUE, details = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L))
res <- list()
res$gbsg_mr_off <- suppressWarnings(do.call(forestsearch, gargs))
res$gbsg_mr_on  <- suppressWarnings(do.call(forestsearch, c(gargs, list(
  mr_inference = TRUE, mr_inference_args = list(draws = 500L)))))
# GLM (continuous) path through forestsearch()
source("tests/testthat/helper-synthetic-dgm.R")
cd <- .make_continuous_data(N = 400L, MD_harm = 2)
cargs <- .fs_args_for("continuous", confounders = c("age", "biomarker", "biomarker_hi", "sex"),
  extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE, hr.threshold = 0.5, hr.consistency = 0.25,
               mr_inference = TRUE, mr_inference_args = list(draws = 300L)))
res$cont_mr_on <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = cd), cargs)))
# fs_mr_inference() direct, supplement S1.7 configuration B (OLS, nested)
set.seed(20260922L); n <- 400L
A <- rbinom(n, 1L, 0.5); X <- sample.int(4L, n, replace = TRUE)
Y <- c(0, 0, -1, -1)[X] * A + rnorm(n)
dB <- data.frame(id = seq_len(n), Y = Y, A = A, X = X)
cB <- list(g1 = which(X == 1L), g2 = which(X <= 2L))
spec <- list(outcome_type = "continuous", effect_measure = "MD", treat.name = "A",
  outcome.name = "Y", event.name = NULL, offset.name = NULL, adjust_covariates = NULL,
  adverse_outcome = TRUE)
res$mr_direct_B <- fs_mr_inference(dB, cB, spec, selected_members = cB$g1,
  admission = list(effect_floor = NULL, consistency = list(c_cons = 0, p_star = 0.9)),
  reselection = "maxeff", draws = 2000L, multiplier = "gaussian",
  ci_method = "ij", seed = 7L)
res <- lapply(res, strip)
dig <- vapply(res, md5, character(1))
# formals of every function defined in the touched files
fns <- unlist(lapply(c("R/fs_mr_inference.R", "R/forestsearch_main.R"), function(f) {
  ex <- parse(f); nm <- vapply(ex, function(e)
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) as.character(e[[2]]) else NA_character_,
    character(1)); nm[!is.na(nm)] }))
fml <- lapply(setNames(fns, fns), function(f) as.list(formals(get(f, asNamespace("forestsearch")))))
saveRDS(list(digest = dig, formals = fml, objects = res, head = system("git rev-parse HEAD", intern = TRUE)), out)
print(dig); cat(length(fml), "functions captured\n")
