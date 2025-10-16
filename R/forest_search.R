
# List of required packages for ForestSearch analysis

required_packages <- c("grf","policytree","data.table","randomForest","survival","weightedSurv","future.apply")
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if(length(missing) > 0) stop("Missing required packages: ", paste(missing, collapse = ", "))


#' Add ID Column to Data Frame
#'
#' Ensures that a data frame has a unique ID column. If \code{id.name} is not provided,
#' a column named "id" is added. If \code{id.name} is provided but does not exist in the data frame,
#' it is created with unique integer values.
#'
#' @param df.analysis Data frame to which the ID column will be added.
#' @param id.name Character. Name of the ID column to add (default is \code{NULL}, which uses "id").
#'
#' @return Data frame with the ID column added if necessary.
#' @export

add_id_column <- function(df.analysis, id.name = NULL) {
  if (is.null(id.name)) {
    df.analysis$id <- seq_len(nrow(df.analysis))
    id.name <- "id"
  } else if (!(id.name %in% names(df.analysis))) {
    df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
  }
  return(df.analysis)
}


#' Generate Prediction Dataset with Subgroup Treatment Recommendation
#'
#' Creates a prediction dataset with a treatment recommendation flag based on subgroup definition.
#'
#' @param df.predict Data frame for prediction (test or validation set).
#' @param sg.harm Character vector of subgroup-defining covariate names.
#' @param version Integer; 1 uses \code{dummy()}, 2 uses \code{dummy2()} for factor encoding.
#'
#' @return Data frame with treatment recommendation flag (\code{treat.recommend}).
#' @export

get_dfpred <- function(df.predict, sg.harm, version = 1) {
  if (version == 1) df.pred <- dummy(df.predict)
  if (version == 2) df.pred <- dummy2(df.predict)
  df.pred$treat.recommend <- NA
  id.harm <- paste(sg.harm, collapse = "==1 & ")
  id.harm <- paste(id.harm, "==1")
  df.pred.0 <- subset(df.pred, eval(parse(text = id.harm)))
  if (nrow(df.pred.0) > 0) {
    df.pred.0$treat.recommend <- 0
  }
  id.noharm <- paste(sg.harm, collapse = "!=1 | ")
  id.noharm <- paste(id.noharm, "!=1")
  df.pred.1 <- subset(df.pred, eval(parse(text = id.noharm)))
  if (nrow(df.pred.1) > 0) {
    df.pred.1$treat.recommend <- 1
  }
  if (nrow(df.pred.0) > 0 && nrow(df.pred.1) > 0) df.pred.out <- data.frame(rbind(df.pred.0, df.pred.1))
  if (nrow(df.pred.0) == 0 && nrow(df.pred.1) > 0) df.pred.out <- data.frame(df.pred.1)
  if (nrow(df.pred.0) > 0 && nrow(df.pred.1) == 0) df.pred.out <- data.frame(df.pred.0)
  return(df.pred.out)
}

#' ForestSearch: Subgroup Identification and Consistency Analysis
#'
#' Performs subgroup identification and consistency analysis for treatment effect heterogeneity using ForestSearch.
#' Supports LASSO-based dimension reduction, GRF-based variable selection, and flexible cut strategies.
#' Returns subgroup definitions, candidate/evaluated confounders, and prediction datasets.
#'
#' @param df.analysis Data frame for analysis.
#' @param outcome.name Character. Name of outcome variable (e.g., time-to-event).
#' @param event.name Character. Name of event indicator variable (0/1).
#' @param treat.name Character. Name of treatment group variable (0/1).
#' @param id.name Character. Name of ID variable.
#' @param potentialOutcome.name Character. Name of potential outcome variable (optional).
#' @param confounders.name Character vector of confounder variable names.
#' @param parallel_args List. Parallelization arguments (plan, workers).
#' @param df.predict Data frame for prediction (optional).
#' @param df.test Data frame for test set (optional).
#' @param is.RCT Logical. Is the data from a randomized controlled trial?
#' @param seedit Integer. Random seed.
#' @param est.scale Character. Effect scale (e.g., \"hr\").
#' @param use_lasso Logical. Use LASSO for dimension reduction.
#' @param use_grf Logical. Use GRF for variable selection.
#' @param plot.grf Logical.  Flag to return GRF plot
#' @param grf_res List. Precomputed GRF results (optional).
#' @param grf_cuts Character vector of GRF cut expressions (optional).
#' @param max_n_confounders Integer. Maximum number of confounders to evaluate.
#' @param grf_depth Integer. Depth for GRF tree.
#' @param dmin.grf Integer. Minimum subgroup size for GRF.
#' @param frac.tau Numeric. Fraction of tau for GRF horizon.
#' @param conf_force Character vector of forced cut expressions.
#' @param defaultcut_names Character vector of confounders to force default cuts.
#' @param cut_type Character. \"default\" or \"median\" for cut strategy.
#' @param exclude_cuts Character vector of cut expressions to exclude.
#' @param replace_med_grf Logical. Remove median cuts that overlap with GRF cuts.
#' @param cont.cutoff Integer. Cutoff for continuous variable determination.
#' @param conf.cont_medians Character vector of continuous confounders to cut at median.
#' @param conf.cont_medians_force Character vector of additional continuous confounders to force median cut.
#' @param n.min Integer. Minimum subgroup size.
#' @param hr.threshold Numeric. Hazard ratio threshold for subgroup selection.
#' @param hr.consistency Numeric. Hazard ratio threshold for consistency.
#' @param sg_focus Character. Subgroup focus criterion (\"hr\", \"hrMaxSG\", \"hrMinSG\", \"maxSG\", \"minSG\").
#' @param fs.splits Integer. Number of splits for consistency analysis.
#' @param m1.threshold Numeric. Threshold for m1 (default: Inf).
#' @param stop.threshold Numeric. Stopping threshold for subgroup search.
#' @param pconsistency.threshold Numeric. Consistency threshold for subgroup search.
#' @param showten_subgroups Logical. Show top ten subgroups.
#' @param d0.min Integer. Minimum number of events in control.
#' @param d1.min Integer. Minimum number of events in treatment.
#' @param max.minutes Numeric. Maximum minutes for search.
#' @param minp Numeric. Minimum proportion for subgroup.
#' @param details Logical. Print details during execution.
#' @param maxk Integer. Maximum number of subgroup factors.
#' @param by.risk Numeric. Interval for risk table time points.
#' @param plot.sg Logical. Plot subgroups.
#' @param max_subgroups_search Integer. Maximum number of subgroups to search.
#' @param vi.grf.min Numeric. Minimum variable importance for GRF screening.
#' @importFrom stats complete.cases median quantile
#' @importFrom grf causal_survival_forest variable_importance
#' @importFrom data.table data.table
#' @importFrom future.apply future_lapply
#' @importFrom randomForest randomForest
#' @importFrom survival Surv
#' @importFrom weightedSurv df_counting
#' @export

forestsearch <- function(df.analysis,
outcome.name = "tte",
event.name = "event",
treat.name = "treat",
id.name = "id",
potentialOutcome.name=NULL,
confounders.name = NULL,
parallel_args = list(plan = "multisession", workers = 6),
df.predict = NULL,
df.test = NULL,
is.RCT = TRUE, seedit = 8316951,
est.scale="hr",
use_lasso = TRUE,
use_grf = TRUE,
grf_res = NULL,
grf_cuts = NULL,
max_n_confounders = 1000,
grf_depth = 2,
dmin.grf = 12,
frac.tau = 0.6,
conf_force=NULL,
defaultcut_names=NULL,
cut_type="default",
exclude_cuts = NULL,
replace_med_grf = FALSE,
cont.cutoff = 4,
conf.cont_medians = NULL,
conf.cont_medians_force = NULL,
n.min = 60,
hr.threshold = 1.25,
hr.consistency = 1.0,
sg_focus = "hr",
fs.splits = 1000,
m1.threshold = Inf,
stop.threshold = 0.90,
pconsistency.threshold = 0.90,
showten_subgroups = FALSE,
d0.min = 10, d1.min = 10,
max.minutes = 3,
minp = 0.025,
details = FALSE,
maxk = 2,
by.risk = 12,
plot.sg = FALSE,
plot.grf = FALSE,
max_subgroups_search = 10,
vi.grf.min = -0.2){

  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())
  # Check parallel arguments for subgroup consistency
  if(length(parallel_args) > 0){
    allowed_plans <- c("multisession", "multicore", "callr","sequential")
    plan_type <- parallel_args$plan
    n_workers <- parallel_args$workers
    max_cores <- parallel::detectCores()
    if (is.null(plan_type)) stop("parallel_args$plan must be specified.")
    if (!plan_type %in% allowed_plans) {
      stop("parallel_args$plan must be one of: ", paste(allowed_plans, collapse = ", "))
    }
    if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
      parallel_args$workers <- 1
    } else {
      parallel_args$workers <- min(n_workers, max_cores)
    }
  }


if (!exists("df.analysis") | !is.data.frame(df.analysis)){
    stop("df.analysis does not exists or is not a data.frame")
  } else if (exists("df.analysis") && !is.data.frame(df.analysis)){
    df.analysis <- as.data.frame(df.analysis)
    message("Converting df.analysis to data.frame")
  }

# Revised
  # if (is.null(id.name)) {
  #   df.analysis$id <- seq_len(nrow(df.analysis))
  #   id.name <- "id"
  # } else if (!(id.name %in% names(df.analysis))) {
  #   df.analysis[[id.name]] <- seq_len(nrow(df.analysis))
  # }

 df.analysis <-  add_id_column(df.analysis, id.name)

 var_names <- c(confounders.name,outcome.name,event.name,id.name,treat.name,potentialOutcome.name)
  # Ensure all required variables exist in df.analysis
  missing_vars <- setdiff(var_names, names(df.analysis))
  if(length(missing_vars) > 0) {
    stop("The following variables are missing in df.analysis: ", paste(missing_vars, collapse = ", "))
  }


if(!(sg_focus %in% c("hr","hrMaxSG", "hrMinSG", "maxSG", "minSG"))) stop("sg_focus must be either hr, hrMaxSG (maxSG), or hrMinSG (minSG)")

if(plot.sg && is.null(by.risk)) stop("by.risk must be non-null if plot.sg = TRUE")

# reset stop.threshold (>1) to override stopping
if(showten_subgroups)  stop.threshold <- 1.1

if(!(cut_type %in% c("default","median"))) stop("only default and median cuts")

if(all(confounders.name %in% names(df.analysis)) != TRUE) stop("Not all confounders found in dataset")

if(!is.null(defaultcut_names)){
if(all(defaultcut_names %in% names(df.analysis)) != TRUE) stop("Not all confounders for default cuts found in dataset")
}


if(is.null(confounders.name)) stop("Confounder names (confounders.name) required")
if(is.null(outcome.name) || is.null(event.name) || is.null(treat.name)) stop("outcome.name, event.name, and treat.name required (missing at least 1)")
if(is.null(hr.threshold) || is.null(hr.consistency) || is.null(pconsistency.threshold)) stop("hr.threshold, hr.consistency, pconsistency.threshold required (missing at least 1)")

if(sg_focus %in% c("maxSG","minSG")){
if(stop.threshold < pconsistency.threshold) stop.threshold <- pconsistency.threshold
# Stop searching once pconsistency.threshold is met
# Select max/min subgroup size
 }

# Re-set stop.threshold
# Continue searching
# Select max/min subgroup size (pconsistency threshold is minimum)
if(sg_focus %in% c("hrMaxSG","hrMinSG") && stop.threshold < 1.0) stop.threshold <- 1.0

# Sort data by id
df.analysis <- df.analysis[order(df.analysis[[id.name]]), , drop = FALSE]
# Select relevant columns
temp <- df.analysis[, var_names, drop = FALSE]
# Identify complete cases
complete_idx <- complete.cases(temp)
n_excluded <- sum(!complete_idx)
# Report exclusions
if(n_excluded > 0) {
  message("Total excluded by omitting missing data = ", n_excluded)
}
# Keep only complete cases
df.analysis <- temp[complete_idx, , drop = FALSE]

rm("temp")

t.start_all<-proc.time()[3]

# Initialize outputs
grf_plot <- NULL
grf_cuts <- NULL

# If using grf and not populated then run grf
if(use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts))) {
if(details){
cat("GRF stage for cut selection with dmin,tau=",c(dmin.grf, frac.tau),"\n")
}
grf_res <- NULL
grf_res <- try(
grf.subg.harm.survival(data = df.analysis, confounders.name = confounders.name,
outcome.name = outcome.name, RCT = is.RCT,seedit = seedit,maxdepth = grf_depth,
event.name = event.name, id.name = id.name, treat.name = treat.name, n.min = n.min, dmin.grf = dmin.grf,
frac.tau = frac.tau, details = details)
,TRUE)

# Check if grf_res is valid (not a try-error and not NULL)
if (!inherits(grf_res, "try-error") && !is.null(grf_res)) {
  # If no subgroup found
  if (is.null(grf_res$sg.harm)) {
    use_grf <- FALSE
    if (isTRUE(details)) {
      cat("NO GRF cuts meeting delta(RMST): dmin.grf=", dmin.grf, "\n")
    }
  } else {
    # If subgroup found
    # Check for DiagrammeR availability
    if (requireNamespace("DiagrammeR", quietly = TRUE) && plot.grf) {

      #grf_plot <- plot(grf_res$tree, leaf.labels = c("Control", "Treat"))

      grf_plot <- try(plot(grf_res$tree, leaf.labels = c("Control", "Treat")), silent = TRUE)
      if (inherits(grf_plot, "try-error")) grf_plot <- NULL


    } else {
      #warning("Skipping tree plot --> DiagrammeR packaged required to view tree graph")
      if (isTRUE(details)) {
        cat("DiagrammeR or not creating: skipping tree plot.\n")
      }
      grf_plot <- NULL
    }
    grf_cuts <- grf_res$tree.cuts
  }
} else {
  # If grf_res is invalid, ensure outputs are NULL
  grf_plot <- NULL
  grf_cuts <- NULL
}
}

if(use_grf && !exists("grf_cuts")) warning("GRF cuts not found")

get_argsFS <- formals(get_FSdata)
args_FS <- names(get_argsFS)
# align with args_call_all
args_FS_filtered <- args_call_all[names(args_call_all) %in% args_FS]
# In get_FSdata the data source is "df"
args_FS_filtered$df.analysis <- df.analysis
args_FS_filtered$grf_cuts <- grf_cuts

# Remove
# cat("args_FS_filtered","\n")
# print(names(args_FS_filtered))

# FSdata1 <- do.call(get_FSdata, args_FS_filtered)
# print(names(FSdata1))

#FSdata <- try(do.call(get_FSdata, args_FS_filtered), TRUE)

FSdata <- tryCatch(
  do.call(get_FSdata, args_FS_filtered),
  error = function(e) {
    message("Error in forestsearch: ", e$message)
    return(NULL)
  }
)

if(inherits(FSdata,"try-error")){
warning("FSdata failure")
}


if(inherits(FSdata,"try-error")) stop("FSdata error")

if(!inherits(FSdata,"try-error")){
lassoomit <- FSdata$lassoomit
lassokeep <- FSdata$lassokeep
df <- FSdata$df

Y <- df[,outcome.name]
Event <- df[,event.name]
Treat <- df[,treat.name]

FSconfounders.name <- FSdata$confs_names
confs_labels <- FSdata$confs
if(is.null(df.predict)) df.predict <- df

if(!is.null(vi.grf.min)){
  # Use GRF for screening and ordering
  # Covariates need to be converted to numeric scale
  # original data.frame version
  X <- as.matrix(df[,FSconfounders.name])
  # Convert to numeric
  X <- apply(X,2,as.numeric)
  tau.rmst <- min(c(max(Y[Treat == 1 & Event == 1]),max(Y[Treat == 0 & Event == 1])))
  # For screening we take 0.9*tau.rms
  if(!is.RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X, Y, Treat, Event, horizon = 0.9 * tau.rmst, seed = 8316951)), TRUE)
  if(is.RCT) cs.forest <- try(suppressWarnings(grf::causal_survival_forest(X, Y, Treat, Event, W.hat = 0.5, horizon = 0.9 * tau.rmst, seed = 8316951)), TRUE)

  vi.cs <- round(grf::variable_importance(cs.forest),4)
  vi.cs2 <- data.frame(confs_labels,FSconfounders.name,vi.cs)
  vi.order <- order(vi.cs,decreasing=TRUE)
  vi.cs2 <- vi.cs2[vi.order,]

  conf.screen <- vi.cs2[,2]
  vi_ratio <- vi.cs2[,3] / max(vi.cs2[,3])
  selected.vars <- which(vi_ratio > vi.grf.min)
  conf.screen <- conf.screen[selected.vars]
  # Restrict to max of max_n_confounders
  # Keeping 1st lmax ordered per vi
  lmax <- min(c(length(conf.screen),max_n_confounders))
  conf.screen_limit <- conf.screen[c(1:lmax)]
  conf.screen <- conf.screen_limit
  if(details){
  cat("Number of factors evaluated=",c(lmax),"\n")
  cat("Confounders per grf screening",conf.screen,"\n")
  vi_res <- vi.cs2[selected.vars,]
  # Re-name for printing
  names(vi_res) <- c("Factors","Labels","VI(grf)")
  print(vi_res)
  }
} else {
conf.screen <- FSconfounders.name
}

df.confounders <- df[,conf.screen]
df.confounders <- dummy(df.confounders)

  # name identification as "id" for merging (df.predict) sg membership flags
  id <- df[,c(id.name)]
  df.fs <- data.frame(Y,Event,Treat,id,df.confounders)
  Z <- as.matrix(df.confounders)
  colnames(Z) <- names(df.confounders)

  find.grps <- subgroup.search(Y = Y, Event = Event, Treat = Treat, Z = Z, d0.min = d0.min, d1.min = d1.min, n.min = n.min,
                               hr.threshold = hr.threshold, max.minutes = max.minutes, details = details, maxk = maxk)

  sg.harm <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  grp.consistency <- NULL

  max_sg_est <- find.grps$max_sg_est

  prop_maxk <- find.grps$prop_max_count

  # If no subgroups found then this is end
  t.end_all<-proc.time()[3]
  t.min_all<-(t.end_all-t.start_all)/60

  if(!is.null(find.grps$out.found) && (any(find.grps$out.found$hr.subgroups$HR > hr.consistency)) || any(find.grps$out.found$hr.subgroups$HR > hr.consistency)){# Found something?
  if(plot.sg && is.null(by.risk)) by.risk<-round(max(Y)/12,0)
  if(details) cat("# of candidate subgroups (meeting all criteria) = ",c(nrow(find.grps$out.found$hr.subgroups)),"\n")

  args_sgc <- names(formals(subgroup.consistency))
  # align with args_call_all
  args_sgc_filtered <- args_call_all[names(args_call_all) %in% args_sgc]
  args_sgc_filtered$df <- df.fs
  args_sgc_filtered$hr.subgroups <- find.grps$out.found$hr.subgroups
  args_sgc_filtered$Lsg <- find.grps$L
  args_sgc_filtered$confs_labels <- confs_labels
  args_sgc_filtered$n.splits <- fs.splits
  args_sgc_filtered$stop_Kgroups <- max_subgroups_search

  grp.consistency <- do.call(subgroup.consistency, args_sgc_filtered)

  t.end_all<-proc.time()[3]
  t.min_all<-(t.end_all-t.start_all)/60

if(details){
cat("Minutes forestsearch overall=",c(t.min_all),"\n")
}
# Found something and a subgroup
if(!is.null(grp.consistency$sg.harm)){
sg.harm <- grp.consistency$sg.harm
# data containing id and treatment flag
temp <- grp.consistency$df_flag
# Merge to analysis data and add treatment flag (all.x=TRUE)
df.est_out <- merge(df, temp, by="id", all.x=TRUE)
 # Return df.predict
if(!is.null(df.predict)){
# This does not work if df.predict is test sample in which case
# cannot match to df_flag by id
df.predict_out <- merge(df.predict, temp, by="id", all.x=TRUE)
}
if(!is.null(df.test)){
  df.test_out <- get_dfpred(df.predict = df.test,sg.harm = grp.consistency$sg.harm,version = 2)
   }
}
 } # Found something
out <- list(grp.consistency = grp.consistency,find.grps = find.grps,
confounders.candidate = FSconfounders.name,
confounders.evaluated = confs_labels,
df.est = df.est_out,
df.predict = df.predict_out,
df.test = df.test_out,
minutes_all = t.min_all,
grf_res = grf_res,
sg_focus = sg_focus,
sg.harm = sg.harm,
grf_cuts = grf_cuts,
prop_maxk = prop_maxk,
max_sg_est = max_sg_est,
grf_plot = grf_plot,
args_call_all = args_call_all)
}
# Fsdata NOT successful
if(inherits(FSdata,"try-error")){
  out <- list(sg.harm=NULL)
}
return(out)
}

