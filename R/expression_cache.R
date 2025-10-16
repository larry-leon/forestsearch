#' Expression Cache Manager for ForestSearch
#'
#' Manages cached parsed expressions and their evaluated results
#' to avoid repeated parsing and evaluation overhead.
#'
#' @export
create_expression_cache <- function() {
  # Create environment for storing cached expressions
  cache <- new.env(hash = TRUE, parent = emptyenv())

  # Track cache statistics (optional, for debugging)
  cache$.stats <- list(
    hits = 0L,
    misses = 0L,
    entries = 0L
  )

  structure(cache, class = "expression_cache")
}

#' Get or Create Cached Expression
#'
#' @param cache Expression cache environment
#' @param expr_text Character string of expression
#' @param parse_only Logical; if TRUE, only cache parsed expression
#' @return Parsed expression or list with parsed and evaluated result
#' @export
get_cached_expression <- function(cache, expr_text, parse_only = TRUE) {
  if (!exists(expr_text, envir = cache)) {
    # Cache miss
    cache$.stats$misses <- cache$.stats$misses + 1L

    parsed_expr <- parse(text = expr_text)

    if (parse_only) {
      cache[[expr_text]] <- parsed_expr
    } else {
      cache[[expr_text]] <- list(
        parsed = parsed_expr,
        text = expr_text
      )
    }
    cache$.stats$entries <- cache$.stats$entries + 1L
    return(cache[[expr_text]])
  } else {
    # Cache hit
    cache$.stats$hits <- cache$.stats$hits + 1L
    return(cache[[expr_text]])
  }
}

#' Evaluate Cached Expression with Data
#'
#' @param cache Expression cache environment
#' @param expr_text Character string of expression
#' @param data Data frame or environment for evaluation
#' @param subset_only Logical; if TRUE, return logical vector, else subset data
#' @return Logical vector or subsetted data
#' @export
eval_cached_expression <- function(cache, expr_text, data, subset_only = FALSE) {
  parsed_expr <- get_cached_expression(cache, expr_text, parse_only = TRUE)

  # Create a hash for this specific data+expression combination
  # This allows caching of evaluation results for repeated data
  data_hash <- paste0(expr_text, "_", digest::digest(data, algo = "xxhash32"))

  if (!exists(data_hash, envir = cache)) {
    result <- eval(parsed_expr, envir = data)

    # Cache the logical vector result (much smaller than full subset)
    if (is.logical(result)) {
      cache[[data_hash]] <- result
    }
  } else {
    result <- cache[[data_hash]]
  }

  if (subset_only) {
    return(result)
  } else {
    return(data[result, , drop = FALSE])
  }
}

#' Clear Expression Cache
#'
#' @param cache Expression cache environment
#' @param pattern Optional pattern to clear specific entries
#' @export
clear_expression_cache <- function(cache, pattern = NULL) {
  if (is.null(pattern)) {
    # Clear everything except stats
    to_remove <- setdiff(ls(envir = cache), ".stats")
    rm(list = to_remove, envir = cache)
    if (exists(".stats", envir = cache)) {
      cache$.stats$entries <- 0L
    }
  } else {
    # Clear entries matching pattern
    to_remove <- grep(pattern, ls(envir = cache), value = TRUE)
    if (length(to_remove) > 0) {
      rm(list = to_remove, envir = cache)
      if (exists(".stats", envir = cache)) {
        cache$.stats$entries <- cache$.stats$entries - length(to_remove)
      }
    }
  }
}

#' Get Cache Statistics
#'
#' @param cache Expression cache environment
#' @return List with cache statistics
#' @export
cache_stats <- function(cache) {
  if (exists(".stats", envir = cache)) {
    stats <- cache$.stats
    stats$hit_rate <- if (stats$hits + stats$misses > 0) {
      stats$hits / (stats$hits + stats$misses)
    } else 0
    return(stats)
  }
  return(NULL)
}
