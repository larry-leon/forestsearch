# diag_phase41_serialization.R --------------------------------------------
# SENTINEL: P41-DIAG-v1-2026-09-02
#
# READ-ONLY diagnostic for the G7-F failure: subgroup_glm() serializes to
# ~3.3 MB on the installed Phase-4.1 build (gate ceiling 512 KB), while a
# minimal package with the identical fitter + estimator code serializes at
# 8.4 KB under both R CMD INSTALL and devtools::install(quick = TRUE) in
# the container (R 4.3.3), with and without --with-keep.source (keep-source
# adds only per-file srcrefs, ~30 KB). This script names the carrier by
# measuring every layer of the object on the machine where the bloat
# reproduces. It modifies nothing and always exits 0.

suppressPackageStartupMessages(library(forestsearch))
S  <- function(x) tryCatch(length(serialize(x, NULL)),
                           error = function(e) NA_integer_)
KB <- function(b) sprintf("%10.1f KB", b / 1024)
hr <- function(t) cat("\n--", t, paste(rep("-", max(1, 60 - nchar(t))),
                                       collapse = ""), "\n")
say <- function(lbl, b) cat(sprintf("%-46s %s\n", lbl, KB(b)))

hr("session")
cat(R.version.string, "\n")
cat("forestsearch:", as.character(utils::packageVersion("forestsearch")),
    "| subgroup_glm exported:",
    "subgroup_glm" %in% getNamespaceExports("forestsearch"), "\n")
cat("keep.source:", getOption("keep.source"),
    "| keep.source.pkgs:", getOption("keep.source.pkgs"), "\n")

hr("headline sizes")
fg <- subgroup_glm()
say("subgroup_glm() whole fitter", S(fg))
cox <- subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim)
say("subgroup_cox() fitter (comparison)", S(cox))
ef_direct <- forestsearch:::.make_lm_estimator_fast("treat_sim", "y_sim",
                                                    adverse_outcome = TRUE)
say("fast estimator, direct from namespace", S(ef_direct))
sl_direct <- make_effect_estimator("continuous", "treat_sim", "y_sim")
say("slow estimator, direct from namespace", S(sl_direct))

hr("attributes vs function vs environment")
say("attr 'formula'", S(attr(fg, "formula")))
say("attr 'effect'",  S(attr(fg, "effect")))
g <- fg; attributes(g) <- NULL
say("fitter with ALL attributes stripped", S(g))
g2 <- g; environment(g2) <- globalenv()
say("fitter, attrs stripped, env -> globalenv", S(g2))
env_f <- environment(fg)
say("env_f serialized alone", S(env_f))
for (nm in ls(env_f, all.names = TRUE)) {
  say(paste0("  env_f$", nm), S(get(nm, envir = env_f)))
}

hr("est_fn dissection")
ef <- get("est_fn", envir = env_f)
cat("est_fn has srcref:", !is.null(attr(ef, "srcref")), "\n")
say("est_fn as stored", S(ef))
h <- ef; attr(h, "srcref") <- NULL
say("est_fn, srcref attr nulled", S(h))
h2 <- tryCatch(utils::removeSource(ef), error = function(e) ef)
say("est_fn, removeSource()", S(h2))
h3 <- ef; environment(h3) <- globalenv()
say("est_fn, env -> globalenv", S(h3))
h4 <- ef; environment(h4) <- globalenv(); h4 <- utils::removeSource(h4)
say("est_fn, env -> globalenv + removeSource", S(h4))
if (!is.null(attr(ef, "srcref"))) {
  sf <- attr(attr(ef, "srcref"), "srcfile")
  cat("srcfile:", tryCatch(sf$filename, error = function(e) "?"),
      "| lines held:",
      tryCatch(length(sf$lines), error = function(e) NA), "\n")
}

hr("est_fn execution frame")
E <- environment(ef)
cat("ls(E):", paste(ls(E, all.names = TRUE), collapse = ", "), "\n")
for (nm in ls(E, all.names = TRUE)) {
  v <- tryCatch(get(nm, envir = E), error = function(e) "<unforced/missing>")
  if (identical(v, "<unforced/missing>")) {
    cat(sprintf("%-46s %s\n", paste0("  E$", nm), "  <unforced/missing>"))
  } else say(paste0("  E$", nm), S(v))
}

hr("parent chain from E (non-special envs get inventoried)")
p <- E
for (i in 1:8) {
  p <- parent.env(p)
  nmE <- environmentName(p)
  special <- isNamespace(p) || identical(p, globalenv()) ||
    identical(p, baseenv()) || identical(p, emptyenv()) ||
    grepl("^package:", nmE)
  cat(sprintf("parent %d: '%s' | namespace=%s | special=%s\n",
              i, nmE, isNamespace(p), special))
  if (!special) {
    b <- ls(p, all.names = TRUE)
    cat("  ORDINARY ENV in chain --", length(b), "bindings:\n")
    for (nm in utils::head(b, 40)) {
      v <- tryCatch(get(nm, envir = p), error = function(e) NULL)
      say(paste0("    $", nm), if (is.null(v)) NA_integer_ else S(v))
    }
    say("  << this env serialized whole >>", S(p))
  }
  if (identical(p, emptyenv())) break
}

hr("bytecode / srcref on the outer fitter")
cat("fitter has srcref:", !is.null(attr(fg, "srcref")),
    "| body class:", class(body(fg)), "\n")
f5 <- tryCatch(utils::removeSource(fg), error = function(e) fg)
say("fitter, removeSource() applied again", S(f5))

hr("verdict hints")
cat("Carrier is wherever the KB column stays ~3000 while siblings are\n")
cat("small. Send the full output back; no fix will be applied here.\n")
