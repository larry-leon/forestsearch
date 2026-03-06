#' Generate Synthetic Data using Bootstrap with Perturbation
#'
#' @param data Original dataset to bootstrap from
#' @param continuous_vars Character vector of continuous variable names
#' @param cat_vars Character vector of categorical variable names
#' @param n Number of synthetic observations to generate (default: same as original)
#' @param seed Random seed for reproducibility
#' @param noise_level Noise level for perturbation (0 to 1, default 0.1)
#' @param id_var Optional name of ID variable to regenerate (will be numbered 1:n)
#' @param cat_flip_prob Probability of flipping categorical values (default: noise_level/2)
#' @param preserve_bounds Logical: should continuous variables stay within original bounds? (default: TRUE)
#' @param ordinal_vars Optional character vector of ordinal categorical variables
#'                     (these will be perturbed to adjacent values rather than randomly flipped)
#'
#' @return A data frame with synthetic data
#'
#' @examples
#' \donttest{
#' # Example 1: Using with GBSG dataset
#' synth_gbsg <- generate_bootstrap_synthetic(
#'   data = survival::gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
#'   cat_vars = c("meno", "hormon", "status"),
#'   ordinal_vars = c("grade"),
#'   id_var = "pid",
#'   n = 1000,
#'   seed = 123,
#'   noise_level = 0.15
#' )
#'
#' # Example 2: Using with any dataset
#' my_data <- data.frame(
#'   id = 1:100,
#'   height = rnorm(100, 170, 10),
#'   weight = rnorm(100, 70, 15),
#'   age = sample(20:80, 100, replace = TRUE),
#'   gender = sample(c("M", "F"), 100, replace = TRUE),
#'   education = sample(1:5, 100, replace = TRUE),
#'   smoker = sample(0:1, 100, replace = TRUE)
#' )
#'
#' synth_data <- generate_bootstrap_synthetic(
#'   data = my_data,
#'   continuous_vars = c("height", "weight", "age"),
#'   cat_vars = c("gender", "smoker"),
#'   ordinal_vars = c("education"),
#'   id_var = "id",
#'   n = 150,
#'   seed = 456
#' )
#' }
#'@export
generate_bootstrap_synthetic <- function(data,
                                         continuous_vars,
                                         cat_vars,
                                         n = NULL,
                                         seed = 123,
                                         noise_level = 0.1,
                                         id_var = NULL,
                                         cat_flip_prob = NULL,
                                         preserve_bounds = TRUE,
                                         ordinal_vars = NULL) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!all(continuous_vars %in% names(data))) {
    missing <- continuous_vars[!continuous_vars %in% names(data)]
    stop("Continuous variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (!all(cat_vars %in% names(data))) {
    missing <- cat_vars[!cat_vars %in% names(data)]
    stop("Categorical variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (!is.null(ordinal_vars) && !all(ordinal_vars %in% names(data))) {
    missing <- ordinal_vars[!ordinal_vars %in% names(data)]
    stop("Ordinal variables not found in data: ", paste(missing, collapse = ", "))
  }

  if (noise_level < 0 || noise_level > 1) {
    stop("'noise_level' must be between 0 and 1")
  }

  # Set defaults
  if (is.null(n)) {
    n <- nrow(data)
  }

  if (is.null(cat_flip_prob)) {
    cat_flip_prob <- noise_level / 2
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Bootstrap sample from original data
  boot_indices <- sample(nrow(data), n, replace = TRUE)
  synthetic_data <- data[boot_indices, ]

  # Reset row names
  rownames(synthetic_data) <- NULL

  # ================================================================================
  # Perturb continuous variables
  # ================================================================================

  for (var in continuous_vars) {
    # Check if variable is numeric
    if (!is.numeric(data[[var]])) {
      warning(paste("Variable", var, "is not numeric. Skipping continuous perturbation."))
      next
    }

    # Calculate noise based on variable's standard deviation
    var_sd <- sd(data[[var]], na.rm = TRUE)
    noise_sd <- var_sd * noise_level

    # Generate noise
    noise <- rnorm(n, mean = 0, sd = noise_sd)

    # Add noise to the variable
    synthetic_data[[var]] <- synthetic_data[[var]] + noise

    # Preserve bounds if requested
    if (preserve_bounds) {
      var_min <- min(data[[var]], na.rm = TRUE)
      var_max <- max(data[[var]], na.rm = TRUE)
      synthetic_data[[var]] <- pmax(var_min, pmin(var_max, synthetic_data[[var]]))
    }

    # Round if original variable appears to be integer
    if (all(data[[var]] == round(data[[var]]), na.rm = TRUE)) {
      synthetic_data[[var]] <- round(synthetic_data[[var]])
    }
  }

  # ================================================================================
  # Perturb categorical variables
  # ================================================================================

  for (var in cat_vars) {
    # Determine which indices to flip
    flip_indices <- which(runif(n) < cat_flip_prob)

    if (length(flip_indices) == 0) next

    # Check if this is an ordinal variable
    if (!is.null(ordinal_vars) && var %in% ordinal_vars) {
      # Ordinal perturbation: change to adjacent values
      unique_vals <- sort(unique(data[[var]]))

      for (i in flip_indices) {
        current_val <- synthetic_data[[var]][i]
        current_pos <- which(unique_vals == current_val)

        if (length(current_pos) == 0) next

        # Determine possible new values (adjacent only)
        possible_new <- c()
        if (current_pos > 1) possible_new <- c(possible_new, unique_vals[current_pos - 1])
        if (current_pos < length(unique_vals)) possible_new <- c(possible_new, unique_vals[current_pos + 1])

        if (length(possible_new) > 0) {
          synthetic_data[[var]][i] <- sample(possible_new, 1)
        }
      }

    } else {
      # Non-ordinal categorical perturbation
      unique_vals <- unique(data[[var]])

      if (length(unique_vals) == 2) {
        # Binary variable: flip the value
        if (is.numeric(synthetic_data[[var]])) {
          # For 0/1 binary
          synthetic_data[[var]][flip_indices] <- 1 - synthetic_data[[var]][flip_indices]
        } else if (is.logical(synthetic_data[[var]])) {
          # For TRUE/FALSE
          synthetic_data[[var]][flip_indices] <- !synthetic_data[[var]][flip_indices]
        } else {
          # For factor or character binary
          for (i in flip_indices) {
            current_val <- synthetic_data[[var]][i]
            other_val <- unique_vals[unique_vals != current_val][1]
            synthetic_data[[var]][i] <- other_val
          }
        }
      } else {
        # Multi-category variable: randomly assign different value
        for (i in flip_indices) {
          current_val <- synthetic_data[[var]][i]
          other_vals <- unique_vals[unique_vals != current_val]
          if (length(other_vals) > 0) {
            synthetic_data[[var]][i] <- sample(other_vals, 1)
          }
        }
      }
    }
  }

  # ================================================================================
  # Regenerate ID variable if specified
  # ================================================================================

  if (!is.null(id_var) && id_var %in% names(synthetic_data)) {
    synthetic_data[[id_var]] <- 1:n
  }

  # ================================================================================
  # Handle remaining variables (neither continuous nor categorical specified)
  # ================================================================================

  all_specified_vars <- c(continuous_vars, cat_vars)
  if (!is.null(id_var)) {
    all_specified_vars <- c(all_specified_vars, id_var)
  }

  remaining_vars <- setdiff(names(data), all_specified_vars)

  if (length(remaining_vars) > 0) {
    message(paste("Note: The following variables were not specified as continuous or categorical",
                  "and will be kept as-is from bootstrap sample:"))
    message(paste("  ", paste(remaining_vars, collapse = ", ")))
  }

  return(synthetic_data)
}

# ================================================================================
# Helper function to automatically detect variable types
# ================================================================================

#' Automatically Detect Variable Types in a Dataset
#'
#' Analyzes a data frame to automatically classify variables as continuous or categorical,
#' and returns a subset of the data with specified variables excluded.
#'
#' @param data A data frame to analyze
#' @param max_unique_for_cat Integer. Maximum number of unique values for a numeric
#'   variable to be considered categorical. Default is 10.
#' @param exclude_vars Character vector of variable names to exclude from both
#'   classification and the returned dataset (e.g., ID variables, timestamps).
#'   Default is NULL.
#'
#' @return A list containing:
#'   \item{continuous_vars}{Character vector of variable names classified as continuous}
#'   \item{cat_vars}{Character vector of variable names classified as categorical}
#'   \item{data_subset}{Data frame with exclude_vars columns removed}
#'
#' @details
#' The function classifies variables using the following rules:
#' \itemize{
#'   \item Numeric variables with more than \code{max_unique_for_cat} unique values
#'     are classified as continuous
#'   \item Numeric variables with \code{max_unique_for_cat} or fewer unique values
#'     are classified as categorical
#'   \item Factor, character, and logical variables are always classified as categorical
#'   \item Variables listed in \code{exclude_vars} are omitted from classification
#'     and removed from the returned dataset
#' }
#'
#' @examples
#' \dontrun{
#' example_data <- data.frame(
#'   id = 1:100,
#'   age = rnorm(100, 50, 10),
#'   grade = sample(1:3, 100, replace = TRUE),
#'   status = sample(c("Active", "Inactive"), 100, replace = TRUE),
#'   score = runif(100, 0, 100)
#' )
#' result <- detect_variable_types(example_data,
#'                                  max_unique_for_cat = 10,
#'                                  exclude_vars = "id")
#' result$continuous_vars  # c("age", "score")
#' result$cat_vars         # c("grade", "status")
#' names(result$data_subset)  # c("age", "grade", "status", "score")
#' }
#'
#' @keywords internal

detect_variable_types <- function(data, max_unique_for_cat = 10, exclude_vars = NULL) {

  # Input validation
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.numeric(max_unique_for_cat) || max_unique_for_cat < 1) {
    stop("'max_unique_for_cat' must be a positive integer")
  }

  if (!is.null(exclude_vars) && !is.character(exclude_vars)) {
    stop("'exclude_vars' must be a character vector or NULL")
  }

  # Check if exclude_vars exist in data
  if (!is.null(exclude_vars)) {
    missing_vars <- setdiff(exclude_vars, names(data))
    if (length(missing_vars) > 0) {
      warning("Variables not found in data and will be ignored: ",
              paste(missing_vars, collapse = ", "))
    }
  }

  continuous_vars <- c()
  cat_vars <- c()

  for (var in names(data)) {
    if (!is.null(exclude_vars) && var %in% exclude_vars) {
      next
    }

    if (is.numeric(data[[var]])) {
      n_unique <- length(unique(data[[var]]))
      # If numeric with few unique values, treat as categorical
      if (n_unique <= max_unique_for_cat) {
        cat_vars <- c(cat_vars, var)
      } else {
        continuous_vars <- c(continuous_vars, var)
      }
    } else {
      # Factor, character, or logical variables
      cat_vars <- c(cat_vars, var)
    }
  }

  # Create subset of data without excluded variables
  if (!is.null(exclude_vars)) {
    vars_to_keep <- setdiff(names(data), exclude_vars)
    data_subset <- data[, vars_to_keep, drop = FALSE]
  } else {
    data_subset <- data
  }

  return(list(
    continuous_vars = continuous_vars,
    cat_vars = cat_vars,
    data_subset = data_subset
  ))
}


#' Generate Synthetic GBSG Data using Generalized Bootstrap
#'
#' @param n Number of observations
#' @param seed Random seed
#' @param noise_level Noise level for perturbation
#'
#' @return Synthetic GBSG dataset
#' @export

generate_gbsg_bootstrap_general <- function(n = 686, seed = 123, noise_level = 0.1) {

  # Load gbsg data into local environment
  gbsg <- NULL  # Avoid R CMD check NOTE
  data("gbsg", package = "survival", envir = environment())

  synthetic_data <- generate_bootstrap_synthetic(
    data = gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
    cat_vars = c("meno", "hormon", "status"),
    ordinal_vars = c("grade"),
    id_var = "pid",
    n = n,
    seed = seed,
    noise_level = noise_level,
    preserve_bounds = TRUE
  )

  return(synthetic_data)
}

#' Generate Bootstrap Sample with Added Noise
#'
#' Creates a bootstrap sample from a dataset with controlled noise added to both
#' continuous and categorical variables. This function is useful for generating
#' synthetic datasets that maintain the general structure of the original data
#' while introducing controlled variation.
#'
#' @param data A data frame containing the original dataset to bootstrap from.
#' @param n Integer. Number of observations in the output dataset. If NULL (default),
#'   uses the same number of rows as the input data.
#' @param continuous_vars Character vector of column names to treat as continuous
#'   variables. If NULL (default), automatically detects numeric columns.
#' @param cat_vars Character vector of column names to treat as categorical variables.
#'   If NULL (default), automatically detects factors, logical columns, and numeric
#'   columns with 10 or fewer unique values.
#' @param id_var Character string specifying the name of the ID variable column.
#'   This column will be reset to sequential values (1:n) in the output.
#'   Default is "pid".
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#' @param noise_level Numeric between 0 and 1. Controls the amount of noise added.
#'   For continuous variables, this is multiplied by the standard deviation to
#'   determine noise magnitude. For categorical variables, this is divided by 2
#'   to determine the probability of value changes. Default is 0.1.
#'
#' @return A data frame with the same structure as the input data, containing
#'   bootstrap sampled observations with added noise.
#'
#' @details
#' The function performs the following operations:
#'
#' \subsection{Bootstrap Sampling}{
#'   Samples n observations with replacement from the original dataset.
#' }
#'
#' \subsection{Continuous Variable Noise}{
#'   \itemize{
#'     \item Adds Gaussian noise with standard deviation = original SD × noise_level
#'     \item Constrains values to remain within original variable bounds
#'     \item Preserves integer type for variables that appear to be integers
#'   }
#' }
#'
#' \subsection{Categorical Variable Perturbation}{
#'   \itemize{
#'     \item Changes values with probability = noise_level / 2
#'     \item Binary variables: flips to opposite value
#'     \item Multi-level unordered: randomly selects from other levels
#'     \item Ordered factors: weights selection toward adjacent levels
#'     \item Preserves factor levels and ordering from original data
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Load example dataset
#' data(gbsg, package = "survival")
#'
#' # Basic usage with automatic variable detection
#' synthetic_data <- generate_bootstrap_with_noise(
#'   data = gbsg,
#'   seed = 123
#' )
#'
#' # Specify variables explicitly
#' synthetic_data <- generate_bootstrap_with_noise(
#'   data = gbsg,
#'   n = 1000,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
#'   cat_vars = c("meno", "grade", "hormon", "status"),
#'   id_var = "pid",
#'   seed = 456,
#'   noise_level = 0.15
#' )
#'
#' # Create multiple synthetic datasets
#' synthetic_list <- lapply(1:10, function(i) {
#'   generate_bootstrap_with_noise(data = gbsg, seed = i)
#' })
#' }
#'
#' @note
#' \itemize{
#'   \item The function assumes that categorical variables with numeric encoding
#'     should maintain their numeric type unless they are factors in the input
#'   \item Missing values (NA) are handled appropriately in calculations but
#'     are not imputed
#'   \item For ordered factors or variables named "grade", the perturbation
#'     favors transitions to adjacent levels over distant levels
#' }
#'
#' @seealso
#' \code{\link[base]{sample}} for bootstrap sampling,
#' \code{\link[stats]{rnorm}} for noise generation
#'
#' @importFrom stats rnorm sd
#' @keywords datagen bootstrap simulation
#' @export

generate_bootstrap_with_noise <- function(data,
                                          n = NULL,
                                          continuous_vars = NULL,
                                          cat_vars = NULL,
                                          id_var = "pid",
                                          seed = 123,
                                          noise_level = 0.1) {

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }

  if (!is.null(n) && (!is.numeric(n) || n <= 0 || n != round(n))) {
    stop("'n' must be a positive integer")
  }

  if (!is.numeric(noise_level) || noise_level < 0 || noise_level > 1) {
    stop("'noise_level' must be a numeric value between 0 and 1")
  }

  if (!is.numeric(seed)) {
    stop("'seed' must be numeric")
  }

  set.seed(seed)

  # Default n to number of rows in original data
  if (is.null(n)) {
    n <- nrow(data)
  }

  # Auto-detect variable types if not specified
  if (is.null(continuous_vars)) {
    continuous_vars <- names(data)[vapply(data, is.numeric, logical(1))]
    # Remove id variable if present
    continuous_vars <- setdiff(continuous_vars, id_var)
  }

  if (is.null(cat_vars)) {

    cat_vars <- names(data)[vapply(data, function(x) {
      is.factor(x) || is.logical(x) ||
        (is.numeric(x) && length(unique(x)) <= 10L)
    }, logical(1))]

    # Remove any vars already in continuous_vars and id_var
    cat_vars <- setdiff(cat_vars, c(continuous_vars, id_var))
  }

  # Check that specified variables exist in data
  missing_cont <- setdiff(continuous_vars, names(data))
  missing_cat <- setdiff(cat_vars, names(data))
  if (length(missing_cont) > 0) {
    stop("Continuous variables not found in data: ", paste(missing_cont, collapse = ", "))
  }
  if (length(missing_cat) > 0) {
    stop("Categorical variables not found in data: ", paste(missing_cat, collapse = ", "))
  }

  # Bootstrap sample from original data
  boot_indices <- sample(nrow(data), n, replace = TRUE)
  synthetic_data <- data[boot_indices, ]

  # Add noise to continuous variables
  for (var in continuous_vars) {
    if (var %in% names(synthetic_data)) {
      # Calculate noise based on variable's standard deviation
      noise_sd <- sd(data[[var]], na.rm = TRUE) * noise_level
      noise <- rnorm(n, 0, noise_sd)

      # Add noise and ensure values stay within original bounds
      synthetic_data[[var]] <- synthetic_data[[var]] + noise
      synthetic_data[[var]] <- pmax(min(data[[var]], na.rm = TRUE),
                                    pmin(max(data[[var]], na.rm = TRUE),
                                         synthetic_data[[var]]))

      # Round if original variable appears to be integer
      if (all(data[[var]] == round(data[[var]]), na.rm = TRUE)) {
        synthetic_data[[var]] <- round(synthetic_data[[var]])
      }
    }
  }

  # Perturb categorical variables
  for (var in cat_vars) {
    if (var %in% names(synthetic_data)) {
      # Get unique levels
      if (is.factor(synthetic_data[[var]])) {
        levels_vec <- levels(synthetic_data[[var]])
        current_values <- as.character(synthetic_data[[var]])
      } else {
        levels_vec <- sort(unique(data[[var]]))
        current_values <- synthetic_data[[var]]
      }

      n_levels <- length(levels_vec)
      flip_prob <- noise_level / 2
      flip_indices <- which(runif(n) < flip_prob)

      if (length(flip_indices) > 0) {
        if (n_levels == 2) {
          # Binary variable - flip to the other value
          for (i in flip_indices) {
            current <- current_values[i]
            other_level <- setdiff(levels_vec, current)
            synthetic_data[[var]][i] <- other_level[1]
          }
        } else if (n_levels > 2) {
          # Multi-level categorical - change to a different random value
          for (i in flip_indices) {
            current <- current_values[i]
            # Choose from other possible values
            other_levels <- setdiff(levels_vec, current)

            # If ordered/ordinal, prefer adjacent values with higher probability
            if (is.ordered(synthetic_data[[var]]) || var == "grade") {
              current_idx <- which(levels_vec == current)
              # Create weights favoring adjacent values
              weights <- rep(1, length(other_levels))
              for (j in seq_along(other_levels)) {
                other_idx <- which(levels_vec == other_levels[j])
                distance <- abs(current_idx - other_idx)
                weights[j] <- 1 / (distance + 0.5)  # Closer values get higher weight
              }
              new_value <- sample(other_levels, 1, prob = weights)
            } else {
              # For unordered categorical, uniform random selection
              new_value <- sample(other_levels, 1)
            }

            synthetic_data[[var]][i] <- new_value
          }
        }
      }

      # Preserve factor structure if original was factor
      if (is.factor(data[[var]])) {
        synthetic_data[[var]] <- factor(synthetic_data[[var]], levels = levels_vec)
        if (is.ordered(data[[var]])) {
          synthetic_data[[var]] <- ordered(synthetic_data[[var]], levels = levels_vec)
        }
      }
    }
  }

  # Reset row names
  rownames(synthetic_data) <- NULL

  # Reset ID variable if it exists
  if (id_var %in% names(synthetic_data)) {
    synthetic_data[[id_var]] <- 1:n
  }

  return(synthetic_data)
}
