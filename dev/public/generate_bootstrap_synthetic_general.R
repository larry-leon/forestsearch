# Generalized Bootstrap with Perturbation Function for Synthetic Data Generation
# This function can work with any dataset, not just GBSG

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
#' # Example 1: Using with GBSG dataset
#' library(survival)
#' data(cancer)
#' 
#' synth_gbsg <- generate_bootstrap_synthetic(
#'   data = gbsg,
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
#' @param data A data frame
#' @param max_unique_for_cat Maximum unique values to consider a numeric variable as categorical
#' @param exclude_vars Variables to exclude from classification (e.g., ID variables)
#' 
#' @return A list with continuous_vars and cat_vars
detect_variable_types <- function(data, max_unique_for_cat = 10, exclude_vars = NULL) {
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
  
  return(list(
    continuous_vars = continuous_vars,
    cat_vars = cat_vars
  ))
}

# ================================================================================
# Wrapper function specifically for GBSG dataset
# ================================================================================

#' Generate Synthetic GBSG Data using Generalized Bootstrap
#' 
#' @param n Number of observations
#' @param seed Random seed
#' @param noise_level Noise level for perturbation
#' 
#' @return Synthetic GBSG dataset
generate_gbsg_bootstrap_general <- function(n = 686, seed = 123, noise_level = 0.1) {
  library(survival)
  data(cancer)
  
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

# ================================================================================
# Example usage and validation
# ================================================================================

if (FALSE) {  # Set to TRUE to run examples
  
  # Example 1: GBSG dataset with manual specification
  library(survival)
  data(cancer)
  
  synth_gbsg_manual <- generate_bootstrap_synthetic(
    data = gbsg,
    continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"),
    cat_vars = c("meno", "hormon", "status"),
    ordinal_vars = c("grade"),
    id_var = "pid",
    n = 1000,
    seed = 123,
    noise_level = 0.2
  )
  
  # Example 2: GBSG dataset with automatic detection
  var_types <- detect_variable_types(gbsg, max_unique_for_cat = 3, exclude_vars = "pid")
  
  synth_gbsg_auto <- generate_bootstrap_synthetic(
    data = gbsg,
    continuous_vars = var_types$continuous_vars,
    cat_vars = var_types$cat_vars,
    id_var = "pid",
    n = 500,
    seed = 456
  )
  
  # Example 3: Custom dataset
  set.seed(999)
  custom_data <- data.frame(
    patient_id = 1:200,
    age = round(rnorm(200, 50, 15)),
    bmi = rnorm(200, 25, 5),
    blood_pressure = rnorm(200, 120, 20),
    treatment = sample(c("A", "B", "C"), 200, replace = TRUE),
    severity = sample(1:4, 200, replace = TRUE),
    recovered = sample(c(TRUE, FALSE), 200, replace = TRUE)
  )
  
  synth_custom <- generate_bootstrap_synthetic(
    data = custom_data,
    continuous_vars = c("age", "bmi", "blood_pressure"),
    cat_vars = c("treatment", "recovered"),
    ordinal_vars = c("severity"),
    id_var = "patient_id",
    n = 300,
    seed = 789,
    noise_level = 0.15
  )
  
  # Compare original and synthetic
  print("Original data summary:")
  print(summary(custom_data))
  print("\nSynthetic data summary:")
  print(summary(synth_custom))
}

# ================================================================================
# Validation function to compare original and synthetic datasets
# ================================================================================

#' Compare Original and Synthetic Datasets
#' 
#' @param original Original dataset
#' @param synthetic Synthetic dataset
#' @param continuous_vars Character vector of continuous variable names
#' @param cat_vars Character vector of categorical variable names
#' 
#' @return Prints comparison statistics
compare_datasets_general <- function(original, synthetic, continuous_vars, cat_vars) {
  cat("\n=== Dataset Comparison ===\n")
  cat("Original size:", nrow(original), "x", ncol(original), "\n")
  cat("Synthetic size:", nrow(synthetic), "x", ncol(synthetic), "\n\n")
  
  # Compare continuous variables
  if (length(continuous_vars) > 0) {
    cat("Continuous Variables:\n")
    cat(sprintf("%-15s %10s %10s %10s %10s %10s\n", 
                "Variable", "Orig_Mean", "Synth_Mean", "Orig_SD", "Synth_SD", "Mean_Diff%"))
    cat(rep("-", 75), "\n", sep = "")
    
    for (var in continuous_vars) {
      if (var %in% names(original) && var %in% names(synthetic)) {
        orig_mean <- mean(original[[var]], na.rm = TRUE)
        synth_mean <- mean(synthetic[[var]], na.rm = TRUE)
        orig_sd <- sd(original[[var]], na.rm = TRUE)
        synth_sd <- sd(synthetic[[var]], na.rm = TRUE)
        mean_diff_pct <- abs((synth_mean - orig_mean) / orig_mean * 100)
        
        cat(sprintf("%-15s %10.2f %10.2f %10.2f %10.2f %10.1f\n",
                    var, orig_mean, synth_mean, orig_sd, synth_sd, mean_diff_pct))
      }
    }
  }
  
  # Compare categorical variables
  if (length(cat_vars) > 0) {
    cat("\nCategorical Variables:\n")
    
    for (var in cat_vars) {
      if (var %in% names(original) && var %in% names(synthetic)) {
        cat("\n", var, ":\n", sep = "")
        orig_tab <- table(original[[var]])
        synth_tab <- table(synthetic[[var]])
        
        # Get all unique levels
        all_levels <- union(names(orig_tab), names(synth_tab))
        
        cat(sprintf("%-10s %15s %15s\n", "Level", "Original_Prop", "Synthetic_Prop"))
        cat(rep("-", 45), "\n", sep = "")
        
        for (level in all_levels) {
          orig_prop <- ifelse(level %in% names(orig_tab), 
                             prop.table(orig_tab)[level], 0)
          synth_prop <- ifelse(level %in% names(synth_tab), 
                              prop.table(synth_tab)[level], 0)
          cat(sprintf("%-10s %15.3f %15.3f\n", level, orig_prop, synth_prop))
        }
      }
    }
  }
}
