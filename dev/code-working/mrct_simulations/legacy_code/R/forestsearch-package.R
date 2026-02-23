#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @importFrom graphics par plot plot.new plot.window lines points abline
#' @importFrom graphics legend axis box text mtext title hist barplot rug rect
#' @importFrom stats as.formula coef na.exclude na.omit predict setNames
#' @importFrom utils hasName head tail
#' @importFrom data.table := .SD .GRP .I .N rbindlist setnames setcolorder is.data.table
#' @importFrom future multisession multicore sequential
#' @name forestsearch-imports
#' @rdname forestsearch-package
#' @keywords internal
NULL

# Note: All globalVariables() declarations are consolidated in globals.R
# Note: install.packages is intentionally NOT imported as it's generally not
# appropriate for package code. Use requireNamespace() instead.
