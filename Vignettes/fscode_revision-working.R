

# library(data.table)
# (library(glmnet,quietly=TRUE))
# library(gt)
# suppressMessages(library(randomForest,quietly=TRUE))
# suppressMessages(library(survival,quietly=TRUE))
# suppressMessages(library(grf,quietly=TRUE))
# suppressMessages(library(policytree,quietly=TRUE))
# suppressMessages(library(DiagrammeR,quietly=TRUE))
# suppressMessages(library(data.table,quietly=TRUE))
# suppressMessages(library(dplyr,quietly=TRUE))
# library(future.apply)

library(weightedSurv)

# If working from main directory

#source_fs_functions("R/")

# If from dev subdirectory ("dev/R/")

rm(list=ls())

#dir_path <- "/Users/leolarr2/Library/CloudStorage/OneDrive-MerckSharp&DohmeLLC/documents/GitHub/forestSearchWorking/dev/R_revised"

dir_path <- "/Users/leolarr2/Library/CloudStorage/OneDrive-MerckSharp&DohmeLLC/documents/GitHub/forestsearch/R"


# Source all files in new directory
# be careful not to load files that over-write functions
files <- list.files(dir_path, pattern = "\\.R$", full.names = TRUE)
sapply(files, source)


# List all objects in the Global Environment
objs <- ls(envir = .GlobalEnv)
# If there are objects, check which are functions
if (length(objs) > 0) {
  funcs <- objs[sapply(objs, function(x) is.function(get(x, envir = .GlobalEnv)))]
  count <- length(funcs)
} else {
  count <- 0
}
count
# 76

# List all objects in the 'stats' package environment
objs <- ls('package:forestsearch')
# Filter to keep only functions
funcs <- objs[sapply(objs, function(x) is.function(get(x, envir = asNamespace('forestsearch'))))]
# Count the number of functions
length(funcs)


library(forestsearch)


library(survival)
dfa <- within(gbsg,{
  id <- as.numeric(c(1:nrow(gbsg)))
  # time to months
  tte <- rfstime/30.4375
  grade3 <- ifelse(grade=="3", 1, 0)
  treat <- hormon
  event <- status
})
confounders.name<-c("age","meno","size","grade3","nodes","pgr","er")


library(doFuture)
library(doRNG)
registerDoFuture()
registerDoRNG()

system.time({fs <- forestsearch(df.analysis = dfa,  confounders.name=confounders.name,
                   outcome.name = "tte", treat.name = "treat", event.name = "event",
                   hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                   sg_focus = "hrMaxSG",
                   showten_subgroups = FALSE, details=TRUE,
                   conf_force = NULL,
                   cut_type = "default", use_grf = TRUE, use_lasso = TRUE,
                   maxk = 2,
                   plot.sg = TRUE, by.risk = 12,
                   parallel_args = list(plan="multisession", workers = 12, show_message = TRUE)
                   )
})


res_tabs <- sg_tables(fs, ndecimals = 3)

res_tabs$tab_estimates

res_tabs$sg10_out



#plan(sequential)

temp <- get_FSdata(df.analysis = dfa, confounders.name = confounders.name, cut_type = "default", use_grf = TRUE,
                   outcome.name = "tte", event.name = "event", use_lasso = TRUE, conf_force=c("age <= 65","meno == 0"))



# registerDoRNG() sets up reproducible random number generation for all subsequent foreach parallel operations in the session.
# You do not need to call it before every foreach or dofuture loop.
# If you restart your R session or change the parallel backend, you should call it again.
# Call once per session/backend setup.
# Not needed before every parallel loop.

NB <- 5
t.start <- proc.time()[3]
# Bootstrap bias-correction
# fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs,nb_boots = NB, show_three = TRUE, details = TRUE, reset_paralle_fs = FALSE, boot_workers =3)
fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = TRUE, details = TRUE, reset_parallel_fs = FALSE,
                                         parallel_args = list(plan = "multisession", workers = 1, show_message = TRUE) )
t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes (total) for bootstrap (boots,mins)",c(NB,t.min),"\n")
cat("Projected minutes for 1000",c(t.min*(1000/NB)),"\n")


t.start<-proc.time()[3]
# Kfolds = n (default to n-fold cross-validations)
fs_OOB <- forestsearch_Kfold(fs.est = fs, details = TRUE, parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE))
summary_OOB <- forestsearch_KfoldOut(res=fs_OOB,details=TRUE,outall=TRUE)
table(summary_OOB$SGs_found[,1])
table(summary_OOB$SGs_found[,2])
t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes for N-fold",c(t.min),"\n")


t.start<-proc.time()[3]
fs_ten <- forestsearch_tenfold(fs.est = fs, sims = 5, details = TRUE,
                               parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE) )

print(fs_ten$find_summary)
print(fs_ten$sens_summary)
print(head(fs_ten$sens_out))
print(head(fs_ten$find_out))
t.now<-proc.time()[3]
tall.min<-(t.now-t.start)/60
cat("Minutes (all analyses)",c(tall.min),"\n")


library(foreach)
library(doFuture)
registerDoFuture()
plan(multisession)

# Function to add ID column
add_id_column <- function(df.analysis, id.name = NULL) {
  if (is.null(id.name)) {
    df.analysis$id <- seq_len(nrow(df.analysis))
    id.name <- "id"
  } else if (!(id.name %in% names(df.analysis))) {
    df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
  }
  return(df.analysis)
}

# List of data frames and id names
dfs <- list(data.frame(a = 1:3), data.frame(a = 4:6))
id.names <- list(NULL, "row_id")

# Parallel execution
results <- foreach(i = 1:2, .options.future = list(add = "add_id_column")) %dofuture% {
  add_id_column(dfs[[i]], id.names[[i]])
}
print(results)


library(DiagrammeR)
graph_code <- 'digraph forestsearch_flow {
  node [shape=box, style=filled, fillcolor="#e6f2ff", fontname="Arial"]
  edge [fontname="Arial"]

  forestsearch [label="forestsearch (main entry)", fillcolor="#b3d1ff"]
  get_FSdata [label="get_FSdata"]
  grf_subg [label="grf.subg.harm.survival"]
  subgroup_search [label="subgroup.search"]
  subgroup_consistency [label="subgroup.consistency"]
  get_dfpred [label="get_dfpred"]
  SG_tab_estimates [label="SG_tab_estimates"]
  forestsearch_bootstrap [label="forestsearch_bootstrap_dofuture"]

  # Main flow
  forestsearch -> get_FSdata
  forestsearch -> grf_subg
  forestsearch -> subgroup_search
  forestsearch -> subgroup_consistency
  forestsearch -> get_dfpred
  forestsearch -> SG_tab_estimates
  forestsearch -> forestsearch_bootstrap

  # get_FSdata helpers
  get_FSdata -> lasso_selection
  get_FSdata -> get_conf_force
  get_FSdata -> filter_by_lassokeep
  get_FSdata -> is_continuous
  get_FSdata -> cut_var
  get_FSdata -> process_conf_force_expr
  get_FSdata -> dummy
  get_FSdata -> dummy2

  # grf_subg helpers
  grf_subg -> policy_tree
  grf_subg -> causal_survival_forest
  grf_subg -> double_robust_scores
  grf_subg -> aggregate

  # subgroup_search helpers
  subgroup_search -> get_combinations_info
  subgroup_search -> get_subgroup_membership
  subgroup_search -> get_covs_in
  subgroup_search -> extract_idx_flagredundancy

  # subgroup_consistency helpers
  subgroup_consistency -> sort_subgroups
  subgroup_consistency -> extract_subgroup
  subgroup_consistency -> sg_consistency_out
  subgroup_consistency -> remove_redundant_subgroups
  subgroup_consistency -> get_split_hr
  subgroup_consistency -> setup_parallel_SGcons

  # forestsearch_bootstrap helpers
  forestsearch_bootstrap -> bootstrap_results
  forestsearch_bootstrap -> bootstrap_ystar
  forestsearch_bootstrap -> get_dfRes
  forestsearch_bootstrap -> get_Cox_sg
  forestsearch_bootstrap -> ci_est
  forestsearch_bootstrap -> get_targetEst
  forestsearch_bootstrap -> fit_cox_models
  forestsearch_bootstrap -> build_cox_formula
  forestsearch_bootstrap -> ensure_packages
  forestsearch_bootstrap -> setup_parallel_SGcons

  # SG_tab_estimates helpers
  SG_tab_estimates -> analyze_subgroup
  analyze_subgroup [label="analyze_subgroup"]
  analyze_subgroup -> cox_summary
  analyze_subgroup -> km_summary
  analyze_subgroup -> rmst_calculation
  analyze_subgroup -> calculate_counts
  analyze_subgroup -> calculate_potential_hr
  SG_tab_estimates -> prepare_subgroup_data
  SG_tab_estimates -> format_results
  SG_tab_estimates -> format_CI
}
'
DiagrammeR::grViz(graph_code)



library(DiagrammeR)

flowchart <- grViz("
digraph flowchart {
  graph [rankdir = LR]
  node [shape=box, style=filled, fillcolor=\"#e6f2ff\", fontname=\"Arial\"]
  edge [fontname=\"Arial\"]

  forestsearch [label=\"forestsearch (main entry)\", fillcolor=\"#b3d1ff\"]
  get_FSdata [label=\"get_FSdata\"]
  grf_subg [label=\"grf.subg.harm.survival\"]
  subgroup_search [label=\"subgroup.search\"]
  subgroup_consistency [label=\"subgroup.consistency\"]
  get_dfpred [label=\"get_dfpred\"]
  SG_tab_estimates [label=\"SG_tab_estimates\"]
  forestsearch_bootstrap [label=\"forestsearch_bootstrap_dofuture\"]

  # Main flow
  forestsearch -> get_FSdata
  forestsearch -> grf_subg
  forestsearch -> subgroup_search
  forestsearch -> subgroup_consistency
  forestsearch -> get_dfpred
  forestsearch -> SG_tab_estimates
  forestsearch -> forestsearch_bootstrap

  # get_FSdata helpers
  get_FSdata -> lasso_selection
  get_FSdata -> get_conf_force
  get_FSdata -> filter_by_lassokeep
  get_FSdata -> is_continuous
  get_FSdata -> cut_var
  get_FSdata -> process_conf_force_expr
  get_FSdata -> dummy
  get_FSdata -> dummy2

  # grf_subg helpers
  grf_subg -> policy_tree
  grf_subg -> causal_survival_forest
  grf_subg -> double_robust_scores
  grf_subg -> aggregate

  # subgroup_search helpers
  subgroup_search -> get_combinations_info
  subgroup_search -> get_subgroup_membership
  subgroup_search -> get_covs_in
  subgroup_search -> extract_idx_flagredundancy

  # subgroup_consistency helpers
  subgroup_consistency -> sort_subgroups
  subgroup_consistency -> extract_subgroup
  subgroup_consistency -> sg_consistency_out
  subgroup_consistency -> remove_redundant_subgroups
  subgroup_consistency -> get_split_hr
  subgroup_consistency -> setup_parallel_SGcons

  # forestsearch_bootstrap helpers
  forestsearch_bootstrap -> bootstrap_results
  forestsearch_bootstrap -> bootstrap_ystar
  forestsearch_bootstrap -> get_dfRes
  forestsearch_bootstrap -> get_Cox_sg
  forestsearch_bootstrap -> ci_est
  forestsearch_bootstrap -> get_targetEst
  forestsearch_bootstrap -> fit_cox_models
  forestsearch_bootstrap -> build_cox_formula
  forestsearch_bootstrap -> ensure_packages
  forestsearch_bootstrap -> setup_parallel_SGcons

  # SG_tab_estimates helpers
  SG_tab_estimates -> analyze_subgroup
  analyze_subgroup [label=\"analyze_subgroup\"]
  analyze_subgroup -> cox_summary
  analyze_subgroup -> km_summary
  analyze_subgroup -> rmst_calculation
  analyze_subgroup -> calculate_counts
  analyze_subgroup -> calculate_potential_hr
  SG_tab_estimates -> prepare_subgroup_data
  SG_tab_estimates -> format_results
  SG_tab_estimates -> format_CI
}
", width = 800, height = 400)





fun_list <- c("dfnew_boot","dfnew","calc_cov",  "calculate_counts", "analyze_subgroups", "calculate_potential_hr","ci.est","count.id","CV_sgs",
  "cox_summary","df_counting","double_robust_scores", "extract_subgroup","format_results", "get_targetEst","getci_Cox",
  "getCIs","grf.estimates.out","hrCI_format","km_summary","n_pcnt","plot_subgroup","plot_weighted_km",
  "prepare_subgroup_data","quiet","rmst_calculation","sg_tables","sort_subgroups","SummaryStat","var_summary",
  "get_FSdata", "dummy","run_bootstrap",
  "forestsearch", "forestsearch_bootstrap_dofuture","get_combinations_info",
  "get_dfpred",
  "grf.subg.harm.survival",
  "subgroup.search",
  "subgroup.consistency",
  "lasso_selection",
  "get_Cox_sg",
  "get_conf_force",
  "filter_by_lassokeep",
  "is.continuous",
  "process_conf_force_expr",
  "is_flag_continuous",
  "is_flag_drop", "acm.disjctif",  "acm.util.df2", "acm.util.df", "dummy2","ztrail","one.zero",
  "get_dfRes", "get_subgroup_membership",
  "SG_tab_estimates",
  "prepare_data",
  "run_grf",
  "evaluate_subgroups",
  "summarize_results",
  "clean_data", "qlow", "qhigh","FS_labels","thiscut","get_cut_name",
  "bootstrap_results", "remove_redundant_subgroups", "sg_consistency_out","get_split_hr","cut_var",
  "bootstrap_ystar", "ensure_packages", "fit_cox_models", "build_cox_formula", "cox.formula.boot",
  "format_CI","setup_parallel_SGcons", "get_covs_in", "extract_idx_flagredundancy"
)
