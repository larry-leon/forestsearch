# Replicate 1 of md40 n500 through the template's own chunks with ci_method = "field":
# keep fs.est for the re-selection-map check.
Sys.setenv(FS_MD_CI = "field", FS_MD_NSIMS = "5", FS_MD_CAMPAIGN = "rep1", FS_MD_QUICKRUN = "TRUE")
qmd <- "sim_fs_maxeffCons_mr_field_md_template.qmd"; txt <- readLines(qmd)
chunk <- function(label) { s <- grep(sprintf("^```\\{r %s[,}]", label), txt)[1]
  e <- s + which(grepl("^```\\s*$", txt[(s+1):length(txt)]))[1]; txt[(s+1):(e-1)] }
for (l in c("setup-knobs", "build-dgm", "machinery")) eval(parse(text = chunk(l)), envir = globalenv())
# Re-run the body of record_replicate() for sim 1 but keep fs.est: call forestsearch() directly
# with the same arguments (copy of the recorder's call).
sd_i <- seed_base + 1L; RNGkind("L'Ecuyer-CMRG"); set.seed(sd_i)
df <- simulate_from_glm_dgm(dgm, n = n_sample, seed = sd_i); df[[id_name]] <- df[[id_name]] %||% seq_len(nrow(df))
t0 <- proc.time()[3]
fs.est <- suppressWarnings(suppressMessages(forestsearch(
  df.analysis = df, confounders.name = confounders_analysis,
  outcome.name = outcome_name, treat.name = treat_name, id.name = id_name,
  outcome_type = "continuous", effect_measure = "MD",
  effect.threshold = md_threshold, consistency.threshold = md_consistency,
  pconsistency.threshold = pconsistency, fs.splits = fs_splits,
  n.min = n_min, d0.min = d0_min, d1.min = d1_min, maxk = maxk,
  vi.grf.min = vi_grf_min, sg_focus = sg_focus, selection_rule = selection_rule,
  effect_neighborhood = effect_neighborhood, stop_threshold = stop_threshold,
  consistency_method = consistency_method, conf.cont_jcuts = fs_conf.cont_jcuts,
  use_lasso = use_lasso, use_dina = use_dina, use_grf = use_grf,
  use_twostage = use_twostage, is.RCT = is_rct, adverse_outcome = adverse_outcome,
  details = FALSE, quiet = TRUE, seedit = sd_i, parallel_args = inner_parallel,
  mr_inference = TRUE, mr_inference_args = mr_inference_args)))
cat(sprintf("forestsearch + field: %.1f s\n", proc.time()[3] - t0))
cat("names(fs.est):", paste(names(fs.est), collapse = " "), "\n")
g <- fs.est$mr_inference
cat("names(g):", paste(names(g), collapse = " "), "\n")
cat(sprintf("selected_index %d | selected_label %s | n_family %d | p_hat[sel] %.4f | selection_rate %.4f | sum(p_hat) %.4f\n",
            g$selected_index, g$selected_label, g$n_family, g$reselection$p_hat[g$selected_index], g$selection_rate, sum(g$reselection$p_hat)))
cat("sg.harm:", paste(fs.est$sg.harm, collapse = " & "), "| |Hhat| =", sum(fs.est$grp.consistency$sg.harm.id == 1L), "\n")
saveRDS(list(fs.est = fs.est, df = df, mr_inference_args = mr_inference_args, dgm = dgm), file.path(Sys.getenv("SP"), "rep1_field.rds"))
