# dev/render_vignettes.R
#
# Test-render a single core vignette WITHOUT invoking devtools::check() or
# pkgdown. Fast iteration loop for editing vignettes during development.
#
# IMPORTANT: forestsearch must be INSTALLED (not just load_all()'d). The
# vignettes use doFuture-based parallel chunks; each worker spawns a fresh
# R session and calls library(forestsearch). Workers cannot see anything
# attached via devtools::load_all(). Install with:
#
#   devtools::install()        # from the package root
#
# Usage (interactive, from the package root):
#
#   source("dev/render_vignettes.R")
#   render_vignette("forestsearch")           # short name
#   render_vignette("forestsearch.Rmd")       # with extension
#   render_vignette("vignettes/methodology.Rmd")  # full path
#
# Output (HTML and intermediate knit artifacts) goes to _vignette_test/ at
# the package root. That directory should be in .gitignore and
# .Rbuildignore so it does not pollute the package build.

render_vignette <- function(name,
                            vignettes_dir = "vignettes",
                            output_dir    = "_vignette_test",
                            quiet         = FALSE) {

  # ── Sanity check: forestsearch must be installed (load_all is not enough).
  if (!requireNamespace("forestsearch", quietly = TRUE)) {
    stop(
      "forestsearch is not installed in the active library.\n",
      "Install via `devtools::install()` from the package root before ",
      "running this script.\n",
      "Note: `devtools::load_all()` will NOT work here because doFuture ",
      "worker sessions cannot see attached-only packages.",
      call. = FALSE
    )
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown is required. Install via install.packages('rmarkdown').",
         call. = FALSE)
  }

  # ── Resolve the source path. Accept: base name, base name with extension,
  #    or a full/relative path.
  if (file.exists(name)) {
    src <- name
  } else {
    base <- if (grepl("\\.[Rr]md$", name)) name else paste0(name, ".Rmd")
    src  <- file.path(vignettes_dir, base)
  }
  if (!file.exists(src)) {
    available <- list.files(vignettes_dir, pattern = "\\.Rmd$",
                            full.names = FALSE)
    stop(
      sprintf("Vignette '%s' not found at '%s'.\nAvailable in %s/: %s",
              name, src, vignettes_dir,
              paste(available, collapse = ", ")),
      call. = FALSE
    )
  }

  # ── Create scratch output and intermediates dirs (both at package root).
  #    Keeping intermediates out of vignettes/ avoids polluting that
  #    directory with knit artifacts that R CMD check would flag.
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  intermediates <- file.path(output_dir, "_intermediates")
  dir.create(intermediates, showWarnings = FALSE, recursive = TRUE)

  # ── Render in a fresh environment so leftover objects from a prior
  #    interactive render do not leak into this one.
  message(sprintf(">>> Rendering %s ...", src))
  t0 <- proc.time()

  out <- rmarkdown::render(
    input             = src,
    output_dir        = output_dir,
    intermediates_dir = intermediates,
    quiet             = quiet,
    envir             = new.env(parent = globalenv())
  )

  dt <- (proc.time() - t0)["elapsed"]
  message(sprintf("<<< Done in %.1f s\n    Output: %s", dt, out))

  invisible(out)
}
