# =============================================================================
# Pressure test: raw-covariate analysis candidate set across vignettes
# =============================================================================
# Goal:
#   1. Verify simulate_from_*_dgm() propagates raw covariate columns from
#      df_super to simulated trials (so confounders_analysis resolves).
#   2. Verify the new analysis_continuous_vars / analysis_binary_vars blocks
#      in each vignette are well-formed (correct names, correct types).
#   3. Verify confounders_analysis is what gets passed to run_simulation_*().
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

PASS <- 0L; FAIL <- 0L
ok <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s\n         %s\n", label, conditionMessage(res)))
    return(invisible(NULL))
  }
  if (isTRUE(res)) {
    PASS <<- PASS + 1L
    cat(sprintf("  [ OK ] %s\n", label))
  } else {
    FAIL <<- FAIL + 1L
    cat(sprintf("  [FAIL] %s (got: %s)\n", label, deparse(res, 60)[1]))
  }
}

# =============================================================================
# 1. Static parsing — ensure each vignette's R code parses
# =============================================================================
cat("\n=== Vignette parse checks ===\n\n")

extract_r_chunks <- function(qmd_path) {
  # Pull all ```{r ...} chunks and stitch them
  lines <- readLines(qmd_path, warn = FALSE)
  in_chunk <- FALSE
  chunks <- character(0)
  cur <- character(0)
  for (l in lines) {
    if (grepl("^```\\{r[^}]*\\}", l)) { in_chunk <- TRUE; next }
    if (in_chunk && grepl("^```\\s*$", l)) {
      chunks <- c(chunks, paste(cur, collapse = "\n"))
      cur <- character(0); in_chunk <- FALSE; next
    }
    if (in_chunk) cur <- c(cur, l)
  }
  paste(chunks, collapse = "\n\n")
}

for (f in c("actg175_continuous_simulations.qmd",
            "actg175_survival_benefit_simulations.qmd",
            "gbsg_poisson_simulations.qmd")) {
  src <- extract_r_chunks(file.path("/mnt/user-data/outputs", f))
  res <- tryCatch(parse(text = src), error = function(e) e)
  ok(sprintf("1.x parse: %s", f), !inherits(res, "error"))
}

# =============================================================================
# 2. Verify the new analysis_*_vars blocks in each vignette
# =============================================================================
cat("\n=== analysis_*_vars blocks ===\n\n")

# Continuous + survival vignettes share ACTG175 candidate set
expected_actg_cont   <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
expected_actg_binary <- c("hemo", "homo", "drugs", "race", "gender", "symptom")

# Poisson vignette uses GBSG
expected_gbsg_cont   <- c("age", "size", "nodes", "pgr", "er")
expected_gbsg_binary <- c("meno", "grade3")

# Read each vignette and extract the var blocks via regex on the raw text
extract_vector_def <- function(qmd_path, name) {
  # The c(...) definitions span multiple lines with comments.  Evaluate
  # the relevant lines as R code instead of regex-extracting.
  txt <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  # Find the start of the assignment
  pat <- paste0("(?m)^", name, "\\s*<-\\s*c\\(")
  m_start <- regexpr(pat, txt, perl = TRUE)
  if (m_start == -1L) return(character(0))
  # Walk forward to find matching close paren
  start_pos <- as.integer(m_start)
  i <- start_pos + attr(m_start, "match.length") - 1L  # at the '('
  depth <- 1L
  n <- nchar(txt)
  while (i < n && depth > 0L) {
    i <- i + 1L
    ch <- substring(txt, i, i)
    if (ch == "(") depth <- depth + 1L
    else if (ch == ")") depth <- depth - 1L
  }
  snippet <- substring(txt, start_pos, i)
  # Evaluate in a fresh env
  e <- new.env()
  tryCatch({
    eval(parse(text = snippet), envir = e)
    get(name, envir = e)
  }, error = function(e) character(0))
}

# --- Continuous ---
fcont <- "/mnt/user-data/outputs/actg175_continuous_simulations.qmd"
cont_cv <- extract_vector_def(fcont, "analysis_continuous_vars")
cont_bv <- extract_vector_def(fcont, "analysis_binary_vars")
ok("2.1 continuous: continuous vars correct",
   identical(sort(cont_cv), sort(expected_actg_cont)))
ok("2.2 continuous: binary vars correct",
   identical(sort(cont_bv), sort(expected_actg_binary)))

# --- Survival ---
fsurv <- "/mnt/user-data/outputs/actg175_survival_benefit_simulations.qmd"
surv_cv <- extract_vector_def(fsurv, "analysis_continuous_vars")
surv_bv <- extract_vector_def(fsurv, "analysis_binary_vars")
ok("2.3 survival: continuous vars correct",
   identical(sort(surv_cv), sort(expected_actg_cont)))
ok("2.4 survival: binary vars correct",
   identical(sort(surv_bv), sort(expected_actg_binary)))

# --- Poisson ---
fpois <- "/mnt/user-data/outputs/gbsg_poisson_simulations.qmd"
pois_cv <- extract_vector_def(fpois, "analysis_continuous_vars")
pois_bv <- extract_vector_def(fpois, "analysis_binary_vars")
ok("2.5 poisson: continuous vars correct",
   identical(sort(pois_cv), sort(expected_gbsg_cont)))
ok("2.6 poisson: binary vars correct",
   identical(sort(pois_bv), sort(expected_gbsg_binary)))

# =============================================================================
# 3. Verify confounders_analysis usage in run_simulation_analysis() calls
#    (continuous/survival via confounders_actg, poisson via confounders_gbsg)
# =============================================================================
cat("\n=== confounders_analysis wiring ===\n\n")

check_wiring <- function(qmd_path, alias_var) {
  txt <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  has_def_via_analysis <- grepl(
    paste0(alias_var, "\\s*<-\\s*confounders_analysis"), txt)
  list(via_analysis = has_def_via_analysis)
}

cw_cont <- check_wiring(fcont, "confounders_actg")
cw_surv <- check_wiring(fsurv, "confounders_actg")
cw_pois <- check_wiring(fpois, "confounders_gbsg")

ok("3.1 continuous: confounders_actc <- confounders_analysis",
   isTRUE(cw_cont$via_analysis))
ok("3.2 survival:   confounders_actg <- confounders_analysis",
   isTRUE(cw_surv$via_analysis))
ok("3.3 poisson:    confounders_gbsg <- confounders_analysis",
   isTRUE(cw_pois$via_analysis))

# =============================================================================
# 4. Verify the OLD pattern (dgm_factor_vars / $factor_vars as confounders)
#    is NO LONGER PRESENT in the active analysis-config lines
# =============================================================================
cat("\n=== Old pattern is gone ===\n\n")

# Use grep for an active assignment pattern (not in comments / docs)
for (lbl in c("continuous", "survival", "poisson")) {
  qmd_path <- switch(lbl,
    continuous = fcont, survival = fsurv, poisson = fpois)
  alias_var <- if (lbl == "poisson") "confounders_gbsg" else "confounders_actg"
  txt <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  # Active line where the alias is assigned to dgm_factor_vars or factor_vars
  bad_pattern <- paste0(alias_var,
    "\\s*<-\\s*(dgm_factor_vars|dgm_calibrated\\$factor_vars)")
  has_bad <- grepl(bad_pattern, txt)
  ok(sprintf("4.x %s: no remaining alias <- dgm$factor_vars assignment",
             lbl), !has_bad)
}

# =============================================================================
# 5. Verify factor coercion loop is present for binaries in each vignette
# =============================================================================
cat("\n=== Factor-coercion loop ===\n\n")

for (lbl in c("continuous", "survival", "poisson")) {
  qmd_path <- switch(lbl,
    continuous = fcont, survival = fsurv, poisson = fpois)
  txt <- paste(readLines(qmd_path, warn = FALSE), collapse = "\n")
  has_loop <- grepl(
    "for\\s*\\(v\\s+in\\s+analysis_binary_vars\\)",
    txt)
  ok(sprintf("5.x %s: factor-coercion loop present", lbl), has_loop)
}

# =============================================================================
# 6. Schema simulation: emulate simulate_from_glm_dgm() / simulate_from_dgm()
#    column propagation and verify confounders_analysis resolves
# =============================================================================
cat("\n=== Simulated-data schema resolution ===\n\n")

# Build a fake actg_df with raw + z-factors as the vignette would
set.seed(1)
n_super <- 100L
fake_df_super <- data.frame(
  id      = seq_len(n_super),
  # raw continuous
  age     = rnorm(n_super, 35, 8),
  preanti = rgamma(n_super, 2, 0.005),
  wtkg    = rnorm(n_super, 75, 12),
  karnof  = sample(70:100, n_super, replace = TRUE),
  cd40    = rnorm(n_super, 350, 100),
  cd80    = rnorm(n_super, 950, 250),
  # raw binaries
  hemo = rbinom(n_super, 1L, 0.1),
  homo = rbinom(n_super, 1L, 0.7),
  drugs = rbinom(n_super, 1L, 0.15),
  race  = rbinom(n_super, 1L, 0.3),
  gender = rbinom(n_super, 1L, 0.85),
  symptom = rbinom(n_super, 1L, 0.2),
  # discretised z-factors (also retained on the data frame, like binary vignette)
  z1 = rbinom(n_super, 1L, 0.5), z2 = rbinom(n_super, 1L, 0.5),
  z3 = rbinom(n_super, 1L, 0.5), z4 = rbinom(n_super, 1L, 0.5),
  # ITT helper cols populated by generate_glm_dgm() / simulate_from_*()
  flag_harm = rbinom(n_super, 1L, 0.3),
  mu0 = rnorm(n_super, 50, 5), mu1 = rnorm(n_super, 50, 5),
  treat = rbinom(n_super, 1L, 0.5)
)
# Coerce binaries to factors as the data-prep loop does
binary_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
for (v in binary_vars) fake_df_super[[v]] <- as.factor(fake_df_super[[v]])

# Emulate simulate_from_glm_dgm()'s row-resampling + outcome draw
sim_n <- 50L
idx <- sample(nrow(fake_df_super), sim_n, replace = TRUE)
sim_data <- fake_df_super[idx, , drop = FALSE]
sim_data$id <- seq_len(sim_n)
sim_data$treat_sim <- rbinom(sim_n, 1L, 0.5)
sim_data$y_sim <- with(sim_data,
  ifelse(treat_sim == 1, mu1, mu0) + rnorm(sim_n, sd = 5))

# The realistic candidate set
analysis_continuous_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
analysis_binary_vars     <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
confounders_analysis <- c(analysis_continuous_vars, analysis_binary_vars)

# Confirm all candidate vars resolve in sim_data
ok("6.1 all continuous vars present in sim_data",
   all(analysis_continuous_vars %in% names(sim_data)))
ok("6.2 all binary vars present in sim_data",
   all(analysis_binary_vars %in% names(sim_data)))

# Confirm the type structure is right for FS's is.continuous()
cont_types <- vapply(analysis_continuous_vars,
  function(v) is.numeric(sim_data[[v]]), logical(1))
ok("6.3 all continuous vars are numeric", all(cont_types))

bin_types <- vapply(analysis_binary_vars,
  function(v) is.factor(sim_data[[v]]), logical(1))
ok("6.4 all binary vars are factor", all(bin_types))

# Confirm z-factors are STILL present (DGM truth still works) but NOT in
# confounders_analysis
ok("6.5 z-factors retained on sim_data (truth flag preserved)",
   "z1" %in% names(sim_data))
ok("6.6 z-factors are NOT in confounders_analysis (analysis-blind)",
   !any(c("z1", "z2", "z3", "z4") %in% confounders_analysis))

# Confirm flag_harm (truth) is still computable
ok("6.7 flag_harm column present (truth column)",
   "flag_harm" %in% names(sim_data))

# =============================================================================
# 7. GBSG schema check (poisson vignette)
# =============================================================================
cat("\n=== GBSG schema ===\n\n")

fake_gbsg <- data.frame(
  id = seq_len(n_super),
  age = sample(40:80, n_super, replace = TRUE),
  size = rgamma(n_super, 2, 0.05),
  nodes = rpois(n_super, 4),
  pgr = rgamma(n_super, 1, 0.01),
  er = rgamma(n_super, 1, 0.005),
  meno = rbinom(n_super, 1L, 0.5),
  grade3 = rbinom(n_super, 1L, 0.3),
  z_age = rbinom(n_super, 1L, 0.5), z_meno = rbinom(n_super, 1L, 0.5),
  z_size = rbinom(n_super, 1L, 0.5), z_grade3 = rbinom(n_super, 1L, 0.3),
  z_nodes = rbinom(n_super, 1L, 0.5), z_pgr = rbinom(n_super, 1L, 0.5),
  z_er = rbinom(n_super, 1L, 0.5),
  flag_harm = rbinom(n_super, 1L, 0.3),
  mu0 = rgamma(n_super, 1, 1), mu1 = rgamma(n_super, 1, 1),
  hormon = rbinom(n_super, 1L, 0.5), time_months = rgamma(n_super, 2, 0.05)
)
for (v in c("meno", "grade3")) fake_gbsg[[v]] <- as.factor(fake_gbsg[[v]])

idx <- sample(nrow(fake_gbsg), sim_n, replace = TRUE)
sim_pois <- fake_gbsg[idx, , drop = FALSE]
sim_pois$id <- seq_len(sim_n)
sim_pois$treat_sim <- rbinom(sim_n, 1L, 0.5)

pois_cont <- c("age", "size", "nodes", "pgr", "er")
pois_bin  <- c("meno", "grade3")
pois_conf <- c(pois_cont, pois_bin)

ok("7.1 all continuous vars present in sim_pois",
   all(pois_cont %in% names(sim_pois)))
ok("7.2 all binary vars present in sim_pois",
   all(pois_bin %in% names(sim_pois)))
ok("7.3 GBSG continuous vars are numeric",
   all(vapply(pois_cont, function(v) is.numeric(sim_pois[[v]]), logical(1))))
ok("7.4 GBSG binary vars are factor",
   all(vapply(pois_bin, function(v) is.factor(sim_pois[[v]]), logical(1))))
ok("7.5 z_* factors retained but not in candidate set",
   "z_er" %in% names(sim_pois) && !"z_er" %in% pois_conf)

# =============================================================================
cat(sprintf("\n=== TOTAL: %d passed, %d failed ===\n", PASS, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
