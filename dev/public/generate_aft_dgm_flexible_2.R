
# Improved Factor Variable Processing with Largest Value as Reference
# This function demonstrates the corrected approach for handling factor variables

#' Process Factor Variables with Largest Value as Reference
#'
#' @param data The input dataset
#' @param factor_vars Character vector of factor/categorical variable names
#' @param df_work Working dataframe to add processed variables to
#' @param verbose Print processing details
#'
#' @return df_work with added dummy variables

process_factor_variables <- function(data, factor_vars, df_work, verbose = FALSE) {

  # Process factor variables
  for (var in factor_vars) {

    if (verbose) cat("Processing factor variable:", var, "\n")

    # Get unique values and determine the reference level (largest)
    if (is.factor(data[[var]])) {
      # If it's already a factor, get levels
      all_levels <- levels(data[[var]])
    } else {
      # If not a factor, get unique values
      all_levels <- sort(unique(data[[var]]))
    }

    n_levels <- length(all_levels)

    if (verbose) {
      cat("  Levels found:", paste(all_levels, collapse = ", "), "\n")
    }

    if (n_levels == 1) {
      # Only one level - skip this variable
      if (verbose) cat("  Skipping - only one level\n")
      next
    } else if (n_levels == 2) {
      # Binary variable - create single indicator
      # Use largest value as reference (indicator = 1 for smaller value)
      ref_level <- max(all_levels)  # Largest value is reference
      other_level <- min(all_levels)  # Smaller value gets indicator

      if (verbose) {
        cat("  Binary variable - Reference level (omitted):", ref_level, "\n")
        cat("  Creating indicator for:", other_level, "\n")
      }

      # Create indicator: 1 if equal to smaller value, 0 if equal to larger value
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummy variables
      # Use largest value as reference

      # Determine reference level based on type
      if (is.numeric(all_levels)) {
        # For numeric, use maximum value
        ref_level <- max(all_levels)
        other_levels <- setdiff(all_levels, ref_level)
        other_levels <- sort(other_levels)  # Sort remaining levels
      } else {
        # For character/factor, use last in alphabetical order
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
        other_levels <- setdiff(all_levels, ref_level)
        other_levels <- sort(other_levels)  # Sort remaining levels
      }

      if (verbose) {
        cat("  Multi-level variable - Reference level (omitted):", ref_level, "\n")
        cat("  Creating indicators for:", paste(other_levels, collapse = ", "), "\n")
      }

      # Create dummy variables for all levels except reference
      for (i in seq_along(other_levels)) {
        level <- other_levels[i]
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)

        if (verbose) {
          cat("    Created:", dummy_name, "for level", level, "\n")
        }
      }
    }
  }

  return(df_work)
}

# ================================================================================
# Example: Demonstrating the Improved Processing
# ================================================================================

demonstrate_factor_processing <- function() {

  cat("================================================================================\n")
  cat("DEMONSTRATION: Factor Variable Processing with Largest Value as Reference\n")
  cat("================================================================================\n\n")

  # Create example dataset
  set.seed(123)
  example_data <- data.frame(
    # Binary factor (0/1)
    binary_var = sample(0:1, 100, replace = TRUE),

    # Binary factor (character)
    sex = sample(c("F", "M"), 100, replace = TRUE),

    # Multi-level numeric factor
    grade = sample(1:3, 100, replace = TRUE),

    # Multi-level character factor
    stage = sample(c("I", "II", "III", "IV"), 100, replace = TRUE),

    # Ordered factor
    severity = factor(sample(c("mild", "moderate", "severe"), 100, replace = TRUE),
                      levels = c("mild", "moderate", "severe"), ordered = TRUE)
  )

  # Initialize working dataframe
  df_work <- data.frame(id = 1:100)

  # Define factor variables
  factor_vars <- c("binary_var", "sex", "grade", "stage", "severity")

  # Process with verbose output
  df_work <- process_factor_variables(
    data = example_data,
    factor_vars = factor_vars,
    df_work = df_work,
    verbose = TRUE
  )

  cat("\n================================================================================\n")
  cat("RESULTS\n")
  cat("================================================================================\n\n")

  # Show the created variables
  cat("Created dummy variables:\n")
  dummy_cols <- grep("^z_", names(df_work), value = TRUE)
  for (col in dummy_cols) {
    cat(sprintf("  %-20s Mean = %.3f (proportion where indicator = 1)\n",
                col, mean(df_work[[col]])))
  }

  # Verification
  cat("\n================================================================================\n")
  cat("VERIFICATION\n")
  cat("================================================================================\n\n")

  # For grade (1, 2, 3) - should use 3 as reference
  cat("Grade variable (1, 2, 3):\n")
  cat("  Original distribution:\n")
  print(table(example_data$grade))
  cat("\n  Dummy variables created:\n")
  cat("  z_grade_1 (grade==1):", sum(df_work$z_grade_1), "observations\n")
  cat("  z_grade_2 (grade==2):", sum(df_work$z_grade_2), "observations\n")
  cat("  (grade==3 is reference - no dummy)\n")

  # For stage (I, II, III, IV) - should use IV as reference
  cat("\nStage variable (I, II, III, IV):\n")
  cat("  Original distribution:\n")
  print(table(example_data$stage))
  cat("\n  Dummy variables created:\n")
  cat("  z_stage_I (stage==I):", sum(df_work$z_stage_I), "observations\n")
  cat("  z_stage_II (stage==II):", sum(df_work$z_stage_II), "observations\n")
  cat("  z_stage_III (stage==III):", sum(df_work$z_stage_III), "observations\n")
  cat("  (stage==IV is reference - no dummy)\n")

  return(df_work)
}

# ================================================================================
# Integration into generate_aft_dgm_flex: Corrected Section
# ================================================================================

# Here's the corrected section to replace in generate_aft_dgm_flex:

process_covariates_corrected <- function(data, continuous_vars, factor_vars, df_work) {

  # Process continuous variables (unchanged)
  for (var in continuous_vars) {
    df_work[[paste0("z_", var)]] <- scale(data[[var]])[, 1]  # Standardize
  }

  # Process factor variables with LARGEST value as reference
  for (var in factor_vars) {

    # Get unique values
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]]))
    }

    n_levels <- length(all_levels)

    if (n_levels == 1) {
      # Skip variables with only one level
      next

    } else if (n_levels == 2) {
      # Binary variable
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
        other_level <- min(all_levels)
      } else {
        # For character, last alphabetically
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
        other_level <- setdiff(all_levels, ref_level)
      }

      # Create single indicator for non-reference level
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummies with largest as reference
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
      } else {
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
      }

      other_levels <- setdiff(all_levels, ref_level)
      other_levels <- sort(other_levels)

      # Create dummy for each non-reference level
      for (level in other_levels) {
        # Use the actual level value in the variable name for clarity
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)
      }
    }
  }

  return(df_work)
}

# ================================================================================
# Test the implementation
# ================================================================================

if (FALSE) {  # Set to TRUE to run test

  # Run the demonstration
  result <- demonstrate_factor_processing()

  # Test with GBSG data
  library(survival)
  data(cancer)

  cat("\n\n================================================================================\n")
  cat("TEST WITH GBSG DATA\n")
  cat("================================================================================\n\n")

  # For GBSG, grade has values 1, 2, 3
  # Using the corrected approach, grade=3 should be reference

  df_test <- data.frame(
    id = 1:nrow(gbsg),
    y = gbsg$rfstime,
    event = gbsg$status
  )

  df_test <- process_factor_variables(
    data = gbsg,
    factor_vars = c("grade", "meno"),
    df_work = df_test,
    verbose = TRUE
  )

  cat("\nCreated variables for GBSG:\n")
  dummy_cols <- grep("^z_", names(df_test), value = TRUE)
  for (col in dummy_cols) {
    cat(sprintf("  %-15s n=%3d (%.1f%%)\n",
                col, sum(df_test[[col]]), 100*mean(df_test[[col]])))
  }

  # Verify grade encoding
  cat("\nGrade encoding verification:\n")
  cat("  Grade 1:", sum(gbsg$grade == 1), "→ z_grade_1:", sum(df_test$z_grade_1), "\n")
  cat("  Grade 2:", sum(gbsg$grade == 2), "→ z_grade_2:", sum(df_test$z_grade_2), "\n")
  cat("  Grade 3:", sum(gbsg$grade == 3), "→ REFERENCE (no dummy)\n")
}

# ================================================================================
# Summary of Changes
# ================================================================================

cat("\n================================================================================\n")
cat("SUMMARY OF CHANGES\n")
cat("================================================================================\n\n")

cat("The key changes in factor processing:\n\n")

cat("1. BINARY VARIABLES (2 levels):\n")
cat("   - Old: First level is reference, second gets indicator\n")
cat("   - New: LARGEST value is reference, smaller gets indicator\n")
cat("   - Example: For meno (0,1), creates z_meno for meno==0\n\n")

cat("2. MULTI-LEVEL VARIABLES (>2 levels):\n")
cat("   - Old: First level is reference, others get dummies\n")
cat("   - New: LARGEST value is reference, others get dummies\n")
cat("   - Example: For grade (1,2,3):\n")
cat("     * Creates z_grade_1 for grade==1\n")
cat("     * Creates z_grade_2 for grade==2\n")
cat("     * Grade==3 is reference (no dummy)\n\n")

cat("3. CHARACTER VARIABLES:\n")
cat("   - Uses last in alphabetical order as reference\n")
cat("   - Example: For stage (I,II,III,IV), IV is reference\n\n")

cat("4. NAMING CONVENTION:\n")
cat("   - Dummy names now include actual level value\n")
cat("   - Example: z_grade_1 instead of z_grade_1 (was ambiguous)\n")
cat("   - Makes interpretation clearer\n\n")

cat("This approach is more standard in statistical modeling where:\n")
cat("- Higher values often represent more severe/advanced conditions\n")
cat("- Using them as reference makes interpretation easier\n")
cat("- Coefficients represent effects relative to highest level\n")


# Enhanced AFT Data Generating Mechanism with Flexible Subgroup Definitions
# Supports multiple cutpoint specifications: fixed values, quantiles, functions, etc.

#' Generate Synthetic Survival Data using AFT Model with Flexible Subgroups
#'
#' @param data The input dataset (data.frame)
#' @param continuous_vars Character vector of continuous variable names
#' @param factor_vars Character vector of factor/categorical variable names
#' @param outcome_var Name of the outcome/time variable
#' @param event_var Name of the event/status variable
#' @param treatment_var Name of the treatment variable (if NULL, will be simulated)
#' @param subgroup_vars Character vector of variables defining the subgroup (optional)
#' @param subgroup_cuts Named list of cutpoint specifications (see details)
#' @param model Character: "alt" (with subgroup effects) or "null" (no subgroup effects)
#' @param k_treat Numeric: treatment effect modifier (default = 1)
#' @param k_inter Numeric: interaction effect modifier (default = 1)
#' @param n_super Integer: size of super population (default = 5000)
#' @param cens_type Character: "weibull" or "uniform" censoring
#' @param cens_params List: parameters for censoring distribution
#' @param seed Integer: random seed for reproducibility
#' @param verbose Logical: print diagnostic information
#'
#' @details
#' The subgroup_cuts parameter accepts flexible specifications:
#'   - Numeric value: subgroup_cuts = list(er = 20) means er <= 20
#'   - Quantile specification: subgroup_cuts = list(er = list(type = "quantile", value = 0.25))
#'   - Function: subgroup_cuts = list(er = list(type = "function", fun = median))
#'   - Range: subgroup_cuts = list(er = list(type = "range", min = 10, max = 50))
#'   - Multiple conditions: subgroup_cuts = list(er = list(type = "multiple", values = c(10, 20, 30)))
#'   - Custom function: subgroup_cuts = list(er = list(type = "custom", fun = function(x) x <= quantile(x, 0.3)))
#'
#' @return List containing:
#'   - df_super: Super population dataset
#'   - model_params: Model parameters (coefficients, scale, etc.)
#'   - subgroup_info: Information about subgroups
#'   - hazard_ratios: True hazard ratios
#'
#' @examples
#' # Example with various cutpoint specifications
#' dgm <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "pgr", "age"),
#'   subgroup_cuts = list(
#'     er = list(type = "quantile", value = 0.25),        # er <= 25th percentile
#'     pgr = list(type = "function", fun = median),       # pgr <= median
#'     age = list(type = "range", min = 40, max = 60)     # age between 40 and 60
#'   ),
#'   model = "alt",
#'   verbose = TRUE
#' )

generate_aft_dgm_flex <- function(data,
                                 continuous_vars,
                                 factor_vars,
                                 outcome_var,
                                 event_var,
                                 treatment_var = NULL,
                                 subgroup_vars = NULL,
                                 subgroup_cuts = NULL,
                                 model = "alt",
                                 k_treat = 1,
                                 k_inter = 1,
                                 n_super = 5000,
                                 cens_type = "weibull",
                                 cens_params = list(),
                                 seed = 8316951,
                                 verbose = TRUE, standardize = FALSE) {

  # ============================================================================
  # Helper function to process cutpoint specifications
  # ============================================================================

  process_cutpoint <- function(var_data, cut_spec, var_name = "") {
    # If cut_spec is a simple numeric value, treat as fixed cutpoint
    if (is.numeric(cut_spec) && length(cut_spec) == 1) {
      return(var_data <= cut_spec)
    }

    # If it's a list, process based on type
    if (is.list(cut_spec)) {
      cut_type <- cut_spec$type

      if (cut_type == "quantile") {
        # Quantile-based cutpoint
        cutpoint <- quantile(var_data, probs = cut_spec$value, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "function") {
        # Function-based cutpoint (e.g., median, mean)
        cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        return(var_data <= cutpoint)

      } else if (cut_type == "range") {
        # Range-based (between min and max)
        return(var_data >= cut_spec$min & var_data <= cut_spec$max)

      } else if (cut_type == "greater") {
        # Greater than cutpoint
        if (!is.null(cut_spec$value)) {
          cutpoint <- cut_spec$value
        } else if (!is.null(cut_spec$quantile)) {
          cutpoint <- quantile(var_data, probs = cut_spec$quantile, na.rm = TRUE)
        } else if (!is.null(cut_spec$fun)) {
          cutpoint <- cut_spec$fun(var_data, na.rm = TRUE)
        }
        return(var_data > cutpoint)

      } else if (cut_type == "multiple") {
        # Multiple cutpoints (var in specified values)
        return(var_data %in% cut_spec$values)

      } else if (cut_type == "custom") {
        # Custom function that returns logical vector
        return(cut_spec$fun(var_data))

      } else {
        stop(paste("Unknown cutpoint type:", cut_type, "for variable:", var_name))
      }
    }

    # Default: if no specification, use median
    if (is.null(cut_spec)) {
      if (verbose) cat("  Using median as default cutpoint for", var_name, "\n")
      return(var_data <= median(var_data, na.rm = TRUE))
    }

    stop(paste("Invalid cutpoint specification for variable:", var_name))
  }

  # ============================================================================
  # Input Validation
  # ============================================================================

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (!model %in% c("alt", "null")) {
    stop("'model' must be either 'alt' or 'null'")
  }

  if (!cens_type %in% c("weibull", "uniform")) {
    stop("'cens_type' must be either 'weibull' or 'uniform'")
  }

  # Check that required variables exist
  required_vars <- c(outcome_var, event_var)
  if (!is.null(treatment_var)) required_vars <- c(required_vars, treatment_var)

  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Required variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Check continuous and factor variables
  all_covars <- c(continuous_vars, factor_vars)
  missing_covars <- setdiff(all_covars, names(data))
  if (length(missing_covars) > 0) {
    stop("Covariate variables not found in data: ", paste(missing_covars, collapse = ", "))
  }

  # ============================================================================
  # Data Preparation
  # ============================================================================

  set.seed(seed)

  # Create working dataset
  df_work <- data.frame(
    id = 1:nrow(data),
    y = data[[outcome_var]],
    event = ifelse(data[[event_var]] == 1, 1, 0)
  )

  # Add treatment (or simulate if not provided)
  if (!is.null(treatment_var)) {
    df_work$treat <- data[[treatment_var]]
  } else {
    df_work$treat <- rbinom(nrow(data), size = 1, prob = 0.5)
    if (verbose) cat("Treatment variable simulated (50/50 randomization)\n")
  }

  # ============================================================================
  # Process Covariates
  # ============================================================================

  # Process continuous variables
  for (var in continuous_vars) {
 if(standardize){
   df_work[[paste0("z_", var)]] <- scale(data[[var]])[, 1]  # Standardize
 } else {
   df_work[[paste0("z_", var)]] <- data[[var]]
 }
   }

  # # Original version
  # # Process factor variables
  # for (var in factor_vars) {
  #   if (is.factor(data[[var]])) {
  #     # Create dummy variables for factors
  #     dummies <- model.matrix(~ data[[var]] - 1)
  #     # Keep all but first level (reference)
  #     if (ncol(dummies) > 1) {
  #       for (j in 2:ncol(dummies)) {
  #         dummy_name <- paste0("z_", var, "_", j-1)
  #         df_work[[dummy_name]] <- dummies[, j]
  #       }
  #     }
  #   } else {
  #     # Treat as binary
  #     df_work[[paste0("z_", var)]] <- as.numeric(data[[var]])
  #   }
  # }

  # Process factor variables with LARGEST value as reference (UNLESS already binary, then retain)
  for (var in factor_vars) {

    # Get unique values
    if (is.factor(data[[var]])) {
      all_levels <- levels(data[[var]])
    } else {
      all_levels <- sort(unique(data[[var]]))
    }

    n_levels <- length(all_levels)

    if (n_levels == 1) {
      # Skip variables with only one level
      next

    } else if (n_levels == 2) {
      # NOTE: revising to keep as-is (switching max/min roles)
      # Binary variable
      if (is.numeric(all_levels)) {
        ref_level <- min(all_levels)
        other_level <- max(all_levels)
      } else {
        # For character, last alphabetically
        ref_level <- sort(all_levels, decreasing = FALSE)[1]
        other_level <- setdiff(all_levels, ref_level)
      }

      # Create single indicator for non-reference level
      df_work[[paste0("z_", var)]] <- as.numeric(data[[var]] == other_level)

    } else {
      # Multiple levels - create dummies with largest as reference
      if (is.numeric(all_levels)) {
        ref_level <- max(all_levels)
      } else {
        ref_level <- sort(all_levels, decreasing = TRUE)[1]
      }

      other_levels <- setdiff(all_levels, ref_level)
      other_levels <- sort(other_levels)

      # Create dummy for each non-reference level
      for (level in other_levels) {
        # Use the actual level value in the variable name for clarity
        dummy_name <- paste0("z_", var, "_", level)
        df_work[[dummy_name]] <- as.numeric(data[[var]] == level)
      }
    }
  }

  # Add any baseline variables in original data that are not specified in the DGM outcome process
  # Identify processed variables
  processed <- c(continuous_vars, factor_vars, outcome_var, event_var, treatment_var)

  # Find unprocessed variables
  unprocessed <- setdiff(names(data), processed)

  if (length(unprocessed) > 0) {
    # Add them to df_work
    for (var in unprocessed) {
      if (!var %in% names(df_work)) {  # Avoid duplicates
        df_work[[var]] <- data[[var]]
      }
    }

    if (verbose) {
      cat("\nAdded", length(unprocessed), "unprocessed variables:",
          paste(unprocessed, collapse = ", "), "\n")
    }
  }

  # ============================================================================
  # Define Subgroups with Flexible Cutpoints
  # ============================================================================

  df_work$flag_harm <- 0
  interaction_term <- NULL
  subgroup_definitions <- list()

  if (model == "alt" && !is.null(subgroup_vars)) {
    # Create subgroup indicators
    subgroup_indicators <- list()

    if (verbose) {
      cat("\n=== Subgroup Definitions ===\n")
    }

    for (var in subgroup_vars) {
      # Get the cutpoint specification for this variable
      cut_spec <- subgroup_cuts[[var]]

      if (var %in% continuous_vars || is.numeric(data[[var]])) {
        # Process continuous variable with flexible cutpoint
        subgroup_indicators[[var]] <- process_cutpoint(
          var_data = data[[var]],
          cut_spec = cut_spec,
          var_name = var
        )

        # Store the actual cutpoint used for reporting
        if (is.numeric(cut_spec) && length(cut_spec) == 1) {
          actual_cutpoint <- cut_spec
          subgroup_definitions[[var]] <- paste(var, "<=", actual_cutpoint)
        } else if (is.list(cut_spec)) {
          if (cut_spec$type == "quantile") {
            actual_cutpoint <- quantile(data[[var]], probs = cut_spec$value, na.rm = TRUE)
            subgroup_definitions[[var]] <- paste(var, "<=", round(actual_cutpoint, 2),
                                                "(", cut_spec$value*100, "th percentile)")
          } else if (cut_spec$type == "function") {
            actual_cutpoint <- cut_spec$fun(data[[var]], na.rm = TRUE)
            fun_name <- deparse(substitute(cut_spec$fun))
            subgroup_definitions[[var]] <- paste(var, "<=", round(actual_cutpoint, 2),
                                                "(", fun_name, ")")
          } else if (cut_spec$type == "range") {
            subgroup_definitions[[var]] <- paste(cut_spec$min, "<=", var, "<=", cut_spec$max)
          } else if (cut_spec$type == "greater") {
            if (!is.null(cut_spec$value)) {
              actual_cutpoint <- cut_spec$value
            } else if (!is.null(cut_spec$quantile)) {
              actual_cutpoint <- quantile(data[[var]], probs = cut_spec$quantile, na.rm = TRUE)
            }
            subgroup_definitions[[var]] <- paste(var, ">", round(actual_cutpoint, 2))
          } else if (cut_spec$type == "custom") {
            subgroup_definitions[[var]] <- paste(var, "(custom function)")
          }
        }

      } else {
        # For factor variables
        if (!is.null(cut_spec)) {
          if (is.character(cut_spec) || is.factor(cut_spec)) {
            subgroup_indicators[[var]] <- data[[var]] == cut_spec
            subgroup_definitions[[var]] <- paste(var, "==", cut_spec)
          } else if (is.list(cut_spec) && cut_spec$type == "multiple") {
            subgroup_indicators[[var]] <- data[[var]] %in% cut_spec$values
            subgroup_definitions[[var]] <- paste(var, "in",
                                                paste(cut_spec$values, collapse = ", "))
          }
        } else {
          # Default: use first level
          first_level <- levels(as.factor(data[[var]]))[1]
          subgroup_indicators[[var]] <- data[[var]] == first_level
          subgroup_definitions[[var]] <- paste(var, "==", first_level)
        }
      }

      if (verbose) {
        cat("  ", subgroup_definitions[[var]], "\n")
        cat("    Proportion in subgroup:",
            round(mean(subgroup_indicators[[var]], na.rm = TRUE), 3), "\n")
      }
    }

    # Create harm flag (all subgroup conditions met)
    df_work$flag_harm <- as.numeric(Reduce("&", subgroup_indicators))

    # Create interaction term
    if (length(subgroup_indicators) > 0) {
      interaction_term <- df_work$treat * df_work$flag_harm
    }

    if (verbose) {
      cat("\nOverall subgroup (all conditions met):\n")
      cat("  Size:", sum(df_work$flag_harm), "out of", nrow(df_work), "\n")
      cat("  Proportion:", round(mean(df_work$flag_harm), 3), "\n")
    }
  }

  # ============================================================================
  # Fit AFT Model (Weibull)
  # ============================================================================

  # Prepare model matrix
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols)])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Fit Weibull AFT model
  formula_str <- paste("Surv(y, event) ~ ", paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survreg(as.formula(formula_str),
                     data = df_work,
                     dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]  # Intercept
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Apply effect modifiers AFT parameterization
  # gamma["treat"] <- k_treat * gamma["treat"]
  # if ("treat_harm" %in% names(gamma)) {
  #   gamma["treat_harm"] <- k_inter * gamma["treat_harm"]
  # }

  # Weibull parameterization
  b0 <- -gamma / tau

  # Apply effect modifiers Weibull log(hazard-ratio) parameterization
  b0["treat"] <- k_treat * b0["treat"]
  if ("treat_harm" %in% names(b0)) {
    b0["treat_harm"] <- k_inter * b0["treat_harm"]
  }

# Transform to revised gamma

  gamma <- -b0 * tau


  # if (verbose) {
  #   cat("\n=== Model Parameters (AFT, log(T)) ===\n")
  #   cat("Intercept (mu):", round(mu, 3), "\n")
  #   cat("Scale (tau):", round(tau, 3), "\n")
  #   cat("Treatment effect:", round(gamma["treat"], 3), "\n")
  #   if ("treat_harm" %in% names(gamma)) {
  #     cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
  #   }
  #  }

  # ============================================================================
  # Generate Super Population
  # ============================================================================

  # Sample with replacement to create super population
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat

  # Sample indices
  idx_sample <- sample(1:nrow(df_work), size = n_super, replace = TRUE)
  df_super <- df_work[idx_sample, ]


  if (verbose) {
    cat("\nOverall subgroup in super-population (all conditions met):\n")
    cat("  Proportion:", round(mean(df_super$flag_harm), 3), "\n")
  }

  if (verbose) {
    cat("\n=== Model Parameters (AFT, log(T)) ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (tau):", round(tau, 3), "\n")
    cat("Treatment effect:", round(gamma["treat"], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  # Assign treatment
  df_super$treat[1:n_treat] <- 1
  df_super$treat[(n_treat + 1):n_super] <- 0

  # Reset IDs
  df_super$id <- 1:n_super

  # Calculate linear predictors for potential outcomes

  X_super <- as.matrix(df_super[, c("treat", covariate_cols)])
  if (!is.null(interaction_term)) {
    # Recalculate interaction for super population
    df_super$treat_harm <- df_super$treat * df_super$flag_harm
    X_super <- cbind(X_super, treat_harm = df_super$treat_harm)
  }
  # Potential outcomes under treatment
  X_treat <- X_super
  X_treat[, "treat"] <- 1
  if ("treat_harm" %in% colnames(X_treat)) {
    X_treat[, "treat_harm"] <- df_super$flag_harm
  }

  # Potential outcomes under control
  X_control <- X_super
  X_control[, "treat"] <- 0
  if ("treat_harm" %in% colnames(X_control)) {
    X_control[, "treat_harm"] <- 0
  }

  # Linear predictors per AFT
  df_super$lin_pred_1 <- X_treat %*% gamma
  df_super$lin_pred_0 <- X_control %*% gamma
  df_super$lin_pred_obs <- X_super %*% gamma

  ## PO log-hazards excluding baseline
  #  Used for calculating empirical version of CDEs (controlled direct effects)
  #  theta0 = L0'beta
  # If additional confounder "w" is included
  # theta0 <- -c(zmat.0%*%gamma.true + w*gamma_w)/tau.strataO
  # Here:

  df_super$theta_0 <- X_control %*% b0
  df_super$theta_1 <- X_treat %*% b0

  df_super$loghr_po <- with(df_super, theta_1 - theta_0)

  # # theta1 = L1'beta
  # theta1 <- -c(zmat.1%*%gamma.true + w*gamma_w)/tau.strataO

  # Notes:
  # linear predictor (setting treat=1)
  # z.1 <- zmat.1[,"z"]
  # # log(Y) AFT parameterization
  # eta1 <- mu + c(zmat.1%*%gamma.true)+w*gamma_w
  # # Weibull log(hazard ratio) setting treat=1
  # # and excluding mu and w*bw (since taking difference below)
  # phi1 <- (-1)*c(zmat.1%*%gamma.true)/tau.strataO
  # log.Y1 <- eta1 + tau.strataO*epsilon
  #
  #  Setting treat=0
  # z.0 <- zmat.0[,"z"]
  # eta0 <- mu + c(zmat.0%*%gamma.true)+w*gamma_w
  # log.Y0 <- eta0 + tau.strataO*epsilon
  # phi0 <- (-1)*c(zmat.0%*%gamma.true)/tau.strataO
  # Potential outcome log(hr) difference
  #loghr.po <- phi1-phi0
  #ahr_empirical <- with(dfs,exp(mean(loghr.po)))
  #res$AHR <- ahr_empirical

  # ============================================================================
  # Calculate True Hazard Ratios
  # ============================================================================

  # Generate potential outcomes for HR calculation
  epsilon <- log(rexp(n_super))  # Extreme value distribution

  # Under treatment
  logT_1 <- mu + tau * epsilon + df_super$lin_pred_1
  T_1 <- exp(logT_1)

  # Under control
  logT_0 <- mu + tau * epsilon + df_super$lin_pred_0
  T_0 <- exp(logT_0)

  # Calculate empirical hazard ratios
  df_temp <- data.frame(
    time = c(T_1, T_0),
    event = 1,
    treat = c(rep(1, n_super), rep(0, n_super)),
    flag_harm = rep(df_super$flag_harm, 2)
  )

  hr_overall <- exp(coxph(Surv(time, event) ~ treat, data = df_temp)$coefficients)

  hr_results <- list(overall = hr_overall)

  AHR <- with(df_super, exp(mean(loghr_po)))
  hr_results$AHR <- AHR

  AHR_harm <- with(subset(df_super, flag_harm == 1), exp(mean(loghr_po)))
  hr_results$AHR_harm <- AHR_harm


  AHR_no_harm <- with(subset(df_super, flag_harm == 0), exp(mean(loghr_po)))
  hr_results$AHR_no_harm <- AHR_no_harm

  if (model == "alt" && sum(df_super$flag_harm) > 0) {
    hr_harm <- exp(coxph(Surv(time, event) ~ treat,
                        data = subset(df_temp, flag_harm == 1))$coefficients)
    hr_no_harm <- exp(coxph(Surv(time, event) ~ treat,
                           data = subset(df_temp, flag_harm == 0))$coefficients)

    hr_results$harm_subgroup <- hr_harm
    hr_results$no_harm_subgroup <- hr_no_harm

    if (verbose) {
      cat("\n=== Hazard Ratios (super popln)===\n")
      cat("Overall HR:", round(hr_overall, 3), "\n")
      cat("Causal AHR:", round(AHR, 3), "\n")
      cat("Harm subgroup HR:", round(hr_harm, 3), "\n")
      cat("Harm subgroup AHR:", round(AHR_harm, 3), "\n")
      cat("No-harm subgroup HR:", round(hr_no_harm, 3), "\n")
      cat("No-harm subgroup AHR:", round(AHR_no_harm, 3), "\n")
    }
  } else if (verbose) {
    cat("\n=== Hazard Ratios (super popln) ===\n")
    cat("Overall HR:", round(hr_overall, 3), "\n")
    cat("Causal AHR:", round(AHR, 3), "\n")
  }

  # ============================================================================
  # Prepare Censoring Parameters
  # ============================================================================

  cens_model <- NULL

  if (cens_type == "weibull") {
    # Fit censoring model
    X_cens <- as.matrix(df_work[, c("treat", covariate_cols)])

    fit_cens <- survreg(Surv(y, 1 - event) ~ X_cens,
                       data = df_work,
                       dist = "weibull")

    mu_cens <- coef(fit_cens)[1]
    tau_cens <- fit_cens$scale
    gamma_cens <- coef(fit_cens)[-1]

    # Store censoring parameters
    cens_model <- list(
      mu = mu_cens,
      tau = tau_cens,
      gamma = gamma_cens,
      type = "weibull"
    )

    # Calculate censoring linear predictors for super population
    df_super$lin_pred_cens_1 <- X_treat[, 1:ncol(X_cens)] %*% gamma_cens
    df_super$lin_pred_cens_0 <- X_control[, 1:ncol(X_cens)] %*% gamma_cens

  } else if (cens_type == "uniform") {
    # Use provided or default uniform censoring parameters
    if (is.null(cens_params$min) || is.null(cens_params$max)) {
      # Default: use range of observed times
      cens_params$min <- min(df_work$y) * 0.5
      cens_params$max <- max(df_work$y) * 1.5
    }

    cens_model <- list(
      min = cens_params$min,
      max = cens_params$max,
      type = "uniform"
    )
  }

  # ============================================================================
  # Prepare Output
  # ============================================================================

  # Model parameters
  model_params <- list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    censoring = cens_model
  )

  # Subgroup information
  subgroup_info <- list(
    vars = subgroup_vars,
    cuts = subgroup_cuts,
    definitions = subgroup_definitions,
    size = sum(df_super$flag_harm),
    proportion = mean(df_super$flag_harm)
  )

  # Analysis variables (for downstream use)
  analysis_vars <- list(
    continuous = continuous_vars,
    factor = factor_vars,
    covariates = covariate_cols,
    treatment = "treat",
    outcome = "y_sim",
    event = "event_sim"
  )

  # Return comprehensive results
  results <- list(
    df_super = df_super,
    model_params = model_params,
    subgroup_info = subgroup_info,
    hazard_ratios = hr_results,
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )

  class(results) <- c("aft_dgm_flex", "list")

  return(results)
}

# ================================================================================
# Convenience function to find quantile for target subgroup proportion
# ================================================================================

#' Find Quantile for Target Subgroup Proportion
#'
#' @param data Data containing the variable
#' @param var_name Name of the variable
#' @param target_prop Target proportion for the subgroup
#' @param direction "less" for <=, "greater" for >
#' @param tol Tolerance for root finding
#'
#' @return The quantile value that achieves the target proportion

find_quantile_for_proportion <- function(data, var_name, target_prop,
                                        direction = "less", tol = 0.0001) {

  var_data <- data[[var_name]]

  # Objective function
  obj_fun <- function(q) {
    cutpoint <- quantile(var_data, probs = q, na.rm = TRUE)
    if (direction == "less") {
      actual_prop <- mean(var_data <= cutpoint, na.rm = TRUE)
    } else {
      actual_prop <- mean(var_data > cutpoint, na.rm = TRUE)
    }
    return(actual_prop - target_prop)
  }

  # Find root
  result <- uniroot(obj_fun, interval = c(0, 1), tol = tol)

  return(list(
    quantile = result$root,
    cutpoint = quantile(var_data, probs = result$root, na.rm = TRUE),
    actual_proportion = target_prop
  ))
}

# ================================================================================
# The simulation function remains the same as before
# ================================================================================

simulate_from_dgm <- function(dgm,
                             n = NULL,
                             rand_ratio = 1,
                             max_follow = Inf,
                             cens_adjust = 0,
                             seed = NULL) {

  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm"))) {
    stop("dgm must be an object created by generate_aft_dgm_flex()")
  }

  if (!is.null(seed)) set.seed(seed)

  df_super <- dgm$df_super
  params <- dgm$model_params

  # Determine sample size
  if (is.null(n)) {
    df_sim <- df_super
    n <- nrow(df_sim)
  } else {
    # Sample from super population
    n_treat <- round(n * rand_ratio / (1 + rand_ratio))
    n_control <- n - n_treat

    idx_sample <- sample(1:nrow(df_super), size = n, replace = TRUE)
    df_sim <- df_super[idx_sample, ]

    # Reassign treatment
    df_sim$treat_sim <- df_sim$treat
    df_sim$treat_sim[1:n_treat] <- 1
    df_sim$treat_sim[(n_treat + 1):n] <- 0

    # Update linear predictors based on assigned treatment
    df_sim$lin_pred_obs <- ifelse(df_sim$treat_sim == 1,
                                  df_sim$lin_pred_1,
                                  df_sim$lin_pred_0)

    # Reset IDs
    df_sim$id <- 1:n
  }

  # Generate survival times
  epsilon <- log(rexp(n))  # Extreme value distribution
  logT_sim <- params$mu + params$tau * epsilon + df_sim$lin_pred_obs
  T_sim <- exp(logT_sim)

  # Generate censoring times
  if (params$censoring$type == "weibull") {
    # Weibull censoring
    lin_pred_cens <- ifelse(df_sim$treat_sim == 1,
                           df_sim$lin_pred_cens_1,
                           df_sim$lin_pred_cens_0)

    epsilon_cens <- log(rexp(n))
    logC_sim <- params$censoring$mu + cens_adjust +
                params$censoring$tau * epsilon_cens + lin_pred_cens
    C_sim <- exp(logC_sim)

  } else if (params$censoring$type == "uniform") {
    # Uniform censoring
    C_sim <- runif(n,
                  min = params$censoring$min,
                  max = params$censoring$max)
  }

  # Apply administrative censoring
  C_sim <- pmin(C_sim, max_follow)

  # Observed times and events
  df_sim$y_sim <- pmin(T_sim, C_sim)
  df_sim$event_sim <- ifelse(T_sim <= C_sim, 1, 0)
  df_sim$t_true <- T_sim
  df_sim$c_time <- C_sim


  #print(names(df_sim))

  # Add analysis variables
  # analysis_cols <- c("id", "treat", "flag_harm",
  #                   dgm$analysis_vars$covariates,
  #                   "y_sim", "event_sim", "t_true")

  # Keep only necessary columns
  #df_sim <- df_sim[, intersect(analysis_cols, names(df_sim))]

  return(df_sim)
}



# Testing

# Your original request - using quantiles
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),   # er <= 25th percentile
    pgr = list(type = "quantile", value = 0.50)   # pgr <= median
  ),
  model = "alt"
)



# Your original request - using quantiles
dgm <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 0.0
)


result <- find_k_inter_for_target_hr(
  target_hr_harm = 2.0,
  data = gbsg,
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
    continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0
)

# Result: k_inter = 1.00 achieves HR_harm = 2.0



base_params <- list(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),
    meno = 0
  ),
  k_treat = 1.0,
  n_super = 5000  # Using smaller for faster demonstration
)


sensitivity_results <- do.call(sensitivity_analysis_k_inter, c(
  list(
    k_inter_range = c(-1.5, 1.5),
    n_points = 11,
    model = "alt"
  ),
  base_params
))

cat("\nSensitivity results:\n")
print(round(sensitivity_results, 3))



df_gbsg <- gbsg
df_gbsg$tte <- with(gbsg, rfstime/30.4375)

dgm_null <- generate_aft_dgm_flex(
  data = df_gbsg,
  n_super = 5000,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "tte",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 0.0
)

# Testing simulation module

draw_null <- simulate_from_dgm(dgm = dgm_null, n = 700, rand_ratio = 1, max_follow = Inf, seed = 123)

library(weightedsurv)

dfcount <- df_counting(df = draw_null, tte.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim")

par(mfrow=c(1,1))

plot_weighted_km(dfcount, conf.int=TRUE, show.logrank = TRUE, ymax = 1.05)
title(main="Simulated data")

library(forestsearch)

names(draw_null)

confounders.name <- c("z_age","z_er","z_pgr", "z_meno", "z_grade_1", "z_grade_2", "size", "nodes")

library(doFuture)
library(doRNG)

registerDoFuture()
registerDoRNG()

# conf_force 'forces' a specific covariate cut, eg. 'age <= 65'

system.time({fs <- forestsearch(draw_null,  confounders.name=confounders.name,
                                outcome.name = "y_sim", treat.name = "treat_sim", event.name = "event_sim", id.name = "id",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("z_age <= 65", "z_er <= 0", "z_er <= 1", "z_er <= 2","z_er <= 5"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 12, d1.min = 12,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="multisession", workers = 12, show_message = TRUE)
)
})




dgm_alt <- generate_aft_dgm_flex(
  data = df_gbsg,
  n_super = 5000,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "tte",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.2612616),   # er <= 25th percentile
    meno = 0
  ),
  model = "alt",
  k_inter = 1.29
)

# Testing simulation module

draw_alt <- simulate_from_dgm(dgm = dgm_alt, n = 700, rand_ratio = 1, max_follow = Inf, seed = 123)

dfcount <- df_counting(df = draw_alt, tte.name = "y_sim", event.name = "event_sim", treat.name = "treat_sim")

par(mfrow=c(1,1))

plot_weighted_km(dfcount, conf.int=TRUE, show.logrank = TRUE, ymax = 1.05)
title(main="Simulated data")


system.time({fs <- forestsearch(draw_alt,  confounders.name=confounders.name,
                                outcome.name = "y_sim", treat.name = "treat_sim", event.name = "event_sim", id.name = "id",
                                potentialOutcome.name = "loghr_po",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("z_age <= 65", "z_er <= 0", "z_er <= 1", "z_er <= 2","z_er <= 5"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 12, d1.min = 12,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="callr", workers = 12, show_message = TRUE)
)
})


res_tabs <- sg_tables(fs, ndecimals = 3)
res_tabs$sg10_out


res_tabs$tab_estimates

plan("sequential")


NB <- 30

t.start <- proc.time()[3]

fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = FALSE, details = TRUE)

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes (total) for bootstrap (boots,mins)",c(NB,t.min),"\n")
cat("Projected minutes for 1000",c(t.min*(1000/NB)),"\n")


sg_tab <- fs_bc$summary$table
sg_tab

plan("sequential")

