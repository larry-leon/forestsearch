# render_quarto_docs.R
#
# Convenience script to render all forestsearch extended Quarto documentation.
# Run from the repo root (where DESCRIPTION lives).
#
# From RStudio console:
#   source("render_quarto_docs.R")
#
# From terminal:
#   Rscript render_quarto_docs.R
#
# Requirements: Quarto CLI installed and on PATH.
#   Install: https://quarto.org/docs/get-started/
#   Check:   quarto::quarto_version()   # or system("quarto --version")

quarto_dir <- "quarto"

if (!dir.exists(quarto_dir)) {
  stop("Directory '", quarto_dir, "' not found. Run from the repo root.")
}

if (nchar(Sys.which("quarto")) == 0L) {
  stop(
    "Quarto CLI not found on PATH.\n",
    "Install from https://quarto.org/docs/get-started/"
  )
}

message("Rendering all .qmd files in ", quarto_dir, "/ ...")
ret <- system(paste("quarto render", quarto_dir))

if (ret == 0L) {
  output_dir <- file.path(quarto_dir, "_output")
  htmls <- list.files(output_dir, pattern = "\\.html$", full.names = TRUE)
  message("\nDone. Output files:")
  lapply(htmls, function(f) message("  ", f))
  invisible(htmls)
} else {
  stop("quarto render exited with status ", ret)
}
