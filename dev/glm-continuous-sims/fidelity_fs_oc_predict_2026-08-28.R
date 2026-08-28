# =============================================================================
# fidelity_fs_oc_predict_2026-08-28.R -- the extraction gate for fs_oc_predict()
# -----------------------------------------------------------------------------
# An extraction is only correct if it reproduces what it was extracted from.
# This script
#   1. runs the prediction document's own `anchor`, `worked-scenario` and
#      `worked-predictions` chunks (extracted from the .qmd by chunk label --
#      the code is the document's, not a transcription) to rebuild the M = 16
#      scenario family and the document's Step-4 quantities at Rdraw = 2e5,
#      seed = 20260825;
#   2. wraps the chunk's nine interface fields as an fs_oc_family object and
#      calls fs_oc_predict(family = ., n = 500, c1 = 30, c2 = 10,
#      consistency_method = "split", draws = 2e5, seed = 20260825);
#   3. asserts identical() -- not all.equal() -- on every returned quantity.
#
# Run from the repository root:
#   Rscript dev/glm-continuous-sims/fidelity_fs_oc_predict_2026-08-28.R
# Exit status 0 iff every quantity is bit-identical.
# =============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(forestsearch)
  }
  library(mvtnorm)
})

DOC <- file.path("quarto", "simulations", "actg175", "continuous",
                 "analytic_verification_and_prediction_md_harm.qmd")
stopifnot(file.exists(DOC))

# ---- chunk extraction by label ---------------------------------------------
.chunk_code <- function(lines, label) {
  open  <- grep(sprintf("^```\\{r %s[ ,}]", label), lines)
  if (length(open) != 1L) stop("chunk '", label, "' not found exactly once")
  close <- grep("^```\\s*$", lines)
  close <- close[close > open][1L]
  lines[(open + 1L):(close - 1L)]
}
doc_lines <- readLines(DOC, warn = FALSE)

# The chunks read the payload relative to the document's directory.
old_wd <- setwd(dirname(DOC)); on.exit(setwd(old_wd), add = TRUE)
env <- new.env(parent = globalenv())
invisible(capture.output({
  eval(parse(text = .chunk_code(doc_lines, "anchor")),            envir = env)
  eval(parse(text = .chunk_code(doc_lines, "worked-scenario")),   envir = env)
  eval(parse(text = .chunk_code(doc_lines, "worked-predictions")), envir = env)
}))
setwd(old_wd)

cat("Document chunks evaluated: M =", env$M, "; Rdraw =", env$Rdraw, "\n")
stopifnot(env$M == 16L, env$Rdraw == 2e5, env$n_rep == 500,
          env$c1 == 30, env$c2 == 10)

# ---- the family, exactly the chunk's nine fields ----------------------------
fam <- structure(list(
  lab = env$lab, Pg = env$Pg, PQg = env$PQg, beta_g = env$beta_g,
  se_g = env$se_g, sens_g = env$sens_g, spec_g = env$spec_g,
  ovl = env$ovl, M = env$M, PQ = env$PQ),
  class = c("fs_oc_family", "list"))

pred <- fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                      consistency_method = "split", draws = 2e5,
                      seed = 20260825)

# ---- the gate: identical(), quantity by quantity ----------------------------
ref <- list(
  det_rate    = env$det_rate,
  EnH         = env$EnH,
  Esens       = env$Esens,
  Espec       = env$Espec,
  Eppv        = env$Eppv,
  Enpv        = env$Enpv,
  EbetaH      = env$EbetaH,
  Enaive_bias = env$Enaive_bias,
  mass_below  = env$mass_below,
  P1          = env$P1,
  p_sel       = env$p_sel,
  sel_c       = env$sel_c,
  det_rate_se = env$.mc_se_prop(env$det_rate, env$Rdraw)
)
got <- pred[names(ref)]

ok <- TRUE
cat(sprintf("%-12s %-9s %25s %25s\n", "quantity", "identical", "document", "fs_oc_predict"))
for (nm in names(ref)) {
  same <- identical(unname(ref[[nm]]), unname(got[[nm]]))
  ok <- ok && same
  show <- function(v) if (length(v) == 1L) format(v, digits = 17) else
    sprintf("[%d] %s ...", length(v), format(v[1L], digits = 17))
  cat(sprintf("%-12s %-9s %25s %25s\n", nm, if (same) "TRUE" else "FALSE",
              show(ref[[nm]]), show(got[[nm]])))
  if (!same) {
    cat(sprintf("    max |diff| = %.3e\n", max(abs(ref[[nm]] - got[[nm]]))))
  }
}
vec_ref <- unlist(ref[c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv",
                        "EbetaH", "Enaive_bias", "mass_below")])
vec_got <- unlist(got[c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv",
                        "EbetaH", "Enaive_bias", "mass_below")])
cat("\nidentical() on the nine-quantity numeric vector:",
    identical(unname(vec_ref), unname(vec_got)), "\n")
cat("FIDELITY GATE:", if (ok) "PASS (bit-identical)" else "FAIL", "\n")
quit(status = if (ok) 0L else 1L)
