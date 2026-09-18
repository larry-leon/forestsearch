# ============================================================================
# accept_directive_b_fits.R  --  Directive B fit-level baseline (5 seeded fits)
#
# Bootstrap and CV are separate entry points and are simply not invoked.  Thresholds are deliberately NOT passed, so each fit
# exercises the inherited defaults -- the thing Directive B changes.
#
# Usage:  Rscript dev/accept_directive_b_fits.R <out_prefix>
# ============================================================================

suppressMessages(devtools::load_all(".", quiet = TRUE))

out_prefix <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_prefix)) out_prefix <- "dev/directive_b_fits"

sys.source(file.path("tests", "testthat", "helper-synthetic-dgm.R"),
           envir = environment())

# N and effect sizes are raised above the factories' defaults so each fit
# actually selects a subgroup -- a "<none>" selection is a weak digest.
# The continuous cell is the exception and stays "<none>" by construction:
# adverse_outcome defaults to FALSE for continuous, and the factory's
# subgroup effect raises y under treatment, i.e. benefit, so nothing is ever
# screened as harm.  Its digest is still a valid parity check.
dat <- list(
  survival   = .make_survival_data(N = 600L, HR_harm = 2.5, seed = 42L),
  binary     = .make_binary_data(N = 600L, OR_harm = 4.0, seed = 42L),
  continuous = .make_continuous_data(N = 600L, MD_harm = 1.5, seed = 42L)
)
CONF <- list(survival   = c("age", "stage", "sex"),
             binary     = c("age", "biomarker_hi", "sex"),
             continuous = c("age", "biomarker_hi", "sex"))

.fit_args <- function(outcome_type, effect_measure) {
  a <- .fs_args_for(outcome_type, confounders = CONF[[outcome_type]])
  a$df.analysis      <- dat[[outcome_type]]
  a$confounders.name <- CONF[[outcome_type]]
  # Inherited-default path: remove the explicit thresholds the test helper sets.
  a$hr.threshold     <- NULL
  a$hr.consistency   <- NULL
  a$effect_measure   <- NULL
  if (!is.na(effect_measure)) a$effect_measure <- effect_measure
  a$mr_inference     <- FALSE
  a$quiet            <- TRUE
  a$details          <- FALSE
  a[!vapply(a, is.null, logical(1)) | names(a) %in% character(0)]
}

CELLS <- list(
  list(id = "binary_unset",    ot = "binary",     em = NA),
  list(id = "binary_OR",       ot = "binary",     em = "OR"),
  list(id = "binary_RD",       ot = "binary",     em = "RD"),
  list(id = "survival_default", ot = "survival",  em = NA),
  list(id = "continuous_MD",   ot = "continuous", em = "MD")
)

rows <- lapply(CELLS, function(cl) {
  args <- .fit_args(cl$ot, cl$em)
  t0 <- unname(proc.time()[3])
  set.seed(20260918L)
  fit <- tryCatch(suppressWarnings(do.call(forestsearch, args)),
                  error = function(e) structure(conditionMessage(e),
                                                class = "fs_probe_error"))
  secs <- round(unname(proc.time()[3]) - t0, 2)
  if (inherits(fit, "fs_probe_error"))
    return(data.frame(cell = cl$id, resolved = "<error>", sg_harm = as.character(fit),
                      n_selected = NA_integer_, digest = NA_character_,
                      secs = secs, stringsAsFactors = FALSE))
  id  <- fit$grp.consistency$sg.harm.id
  key <- list(sg.harm          = fit$sg.harm,
              sg.harm.id       = id,
              effect_measure   = fit$effect_measure,
              threshold_config = fit$threshold_config[c(
                "effect_measure", "screening", "consistency",
                "screening_natural", "consistency_natural", "scale")])
  data.frame(
    cell       = cl$id,
    resolved   = if (is.null(fit$effect_measure)) "HR" else fit$effect_measure,
    sg_harm    = if (is.null(fit$sg.harm)) "<none>"
                 else paste(fit$sg.harm, collapse = " & "),
    n_selected = if (is.null(id)) NA_integer_ else sum(id == 1L),
    digest     = digest::digest(key, algo = "md5"),
    secs       = secs,
    stringsAsFactors = FALSE
  )
})
fits <- do.call(rbind, rows)

sha <- system("git rev-parse HEAD", intern = TRUE)
dirty <- length(system("git status --porcelain --untracked-files=no",
                       intern = TRUE)) > 0L
saveRDS(fits, paste0(out_prefix, ".rds"))
con <- file(paste0(out_prefix, ".txt"), open = "wt")
writeLines(c("# Directive B fit digests (bootstrap off, mr off, seed 20260918)",
             sprintf("# SHA: %s%s", sha, if (dirty) " (tracked tree DIRTY)" else ""),
             sprintf("# total wall clock: %.1f s", sum(fits$secs)), ""), con)
utils::capture.output(print(fits, right = FALSE), file = con)
close(con)
cat(readLines(paste0(out_prefix, ".txt")), sep = "\n")
