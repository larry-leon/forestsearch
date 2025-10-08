
# Check that this can be removed after testing
## Necessary packages
# Add plyr needed for rbind.fill
# requiredPackages <- c("doRNG","doFuture")
# ## Note: cli only needed for code checking and printing when details=TRUE
# ipak <- function(pkg){
#   new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#   if (length(new.pkg)>0){
#     install.packages(new.pkg, dependencies = TRUE, quietly=TRUE, verbose=FALSE, logical.return=FALSE)
#     sapply(pkg, require, character.only = TRUE, quietly=TRUE)
#   }
# }
# ipak(requiredPackages)
# sapply(c(requiredPackages), require, character.only = TRUE, quietly=TRUE)


#' ForestSearch K-Fold Cross-Validation
#'
#' Performs K-fold cross-validation for ForestSearch, evaluating subgroup identification and agreement.
#'
#' @param fs.est ForestSearch results object.
#' @param Kfolds Integer. Number of folds (default: nrow(fs.est$df.est)).
#' @param seedit Integer. Random seed (default: 8316951).
#' @param parallel_args List. Parallelization arguments (plan, workers, show_message).
#' @param sg0.name Character. Name for subgroup 0 (default: "Not recommend").
#' @param sg1.name Character. Name for subgroup 1 (default: "Recommend").
#' @param details Logical. Print details during execution (default: FALSE).
#'
#' @return List with cross-validation results, arguments, timing, and subgroup agreement metrics.
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @importFrom doFuture %dofuture%
#' @export

forestsearch_Kfold <- function(fs.est, Kfolds = nrow(fs.est$df.est), seedit = 8316951, parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
                               sg0.name = "Not recommend", sg1.name = "Recommend", details = FALSE){
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
  # Note: setup_parallel_SGcons is re-purposed from subgroup_consistency
  setup_parallel_SGcons(parallel_args)

t.startk <- proc.time()[3]

fs_args <- fs.est$args_call_all
get_names <- with(fs_args,c(confounders.name, outcome.name, event.name, id.name, treat.name))
## Extract necessary data
dfa <- fs.est$df.est[,c(get_names,"treat.recommend")]
## re-name treat.recommend to original to compare with k-folds
names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"

outcome.name <- fs_args$outcome.name
event.name <- fs_args$event.name
treat.name <- fs_args$treat.name
id.name <- fs_args$id.name

dfnew <- as.data.frame(dfa)
if(Kfolds==nrow(fs.est$df.est)){
df_scrambled <- copy(dfnew)
}
if(Kfolds < nrow(fs.est$df.est)){
df_scrambled <- copy(dfnew)
set.seed(seedit)
id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
df_scrambled <- df_scrambled[id_sample,]
}
folds <- cut(seq(1,nrow(df_scrambled)),breaks = Kfolds,labels = FALSE)
if(details) cat("Range of unique left-out fold sample sizes",c(range(unique(table(folds)))),"\n")
# turn off details in fs_args and reset parallel_args to turn off parallel
fs_args$parallel_args <- list()
fs_args$details <- FALSE
fs_args$plot.sg <- FALSE
est.scale <- fs_args$est.scale

cv_args <- fs_args

# rbind.fill (plyr) replaced with rbind (remove dependence on plyr)
resCV <- foreach::foreach(
  cv_index = seq_len(Kfolds),
  .options.future = list(seed = TRUE, add = c("df_scrambled","folds","fs.est")),
  .combine = "rbind",
  .errorhandling="pass"
) %dofuture% {
  testIndexes <- which(folds == cv_index,arr.ind = TRUE)
  x.test <- df_scrambled[testIndexes, ]
  x.train <- df_scrambled[-testIndexes, ]
  cv_args$df.analysis <- x.train

  fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))

 if(!inherits(fs.train,"try-error") && !is.null(fs.train$sg.harm)){
    sg1 <- fs.train$sg.harm[1]
    sg2 <- fs.train$sg.harm[2]
    df.test <- get_dfpred(df.predict = x.test, sg.harm=fs.train$sg.harm, version = 2)
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(sg1,nrow(df.test))
    df.test$sg2 <- rep(sg2,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
    } else {
    df.test <- x.test
    df.test$cvindex <- rep(cv_index,nrow(df.test))
    df.test$sg1 <- rep(NA,nrow(df.test))
    df.test$sg2 <- rep(NA,nrow(df.test))
    df.test$treat.recommend <- rep(1.0,nrow(df.test))
    dfres <- data.table(df.test[,c(id.name,outcome.name,event.name,treat.name,"treat.recommend","treat.recommend.original","cvindex","sg1","sg2")])
 }
   return(dfres)
}
t.now<-proc.time()[3]
t.min<-(t.now-t.startk)/60
if(length(unique(resCV$id)) != nrow(df_scrambled)) stop("K-fold cross-validation sample size does not equal full analysis dataset")
chk1 <- c(range(unique(table(folds))))
chk2 <- c(range(unique(table(resCV$cvindex))))
if(!all(chk1  == chk2)) stop("Mismatch on range of fold sample sizes")
if(details) cat("Minutes for Cross-validation: (Kfolds,minutes)=",c(Kfolds,t.min),"\n")
resCV <- as.data.frame(resCV)
if(est.scale == "1/hr"){
resCV$treat2 <- 1-resCV[,treat.name]
treat.name <- c("treat2")
sg0.name <- "Recommend"
sg1.name <- "Not recommend"
}
sg_found_count <- sum(ifelse(!is.na(resCV$sg1) | !is.na(resCV$sg2),1,0))
# Proportion of SG's found relative to number of folds (Kfolds)
propn_SG <- 100 * round(c(sg_found_count/Kfolds),3)
if(details && Kfolds == nrow(fs.est$df.est)) cat("% of Kfolds where subgroup was found (valid for n-fold)",c(paste0(propn_SG,"%")),"\n")
return(
list(resCV = resCV, cv_args = cv_args, timing_minutes = t.min,
prop_SG_found=propn_SG, sg_analysis = fs.est$sg.harm, sg0.name = sg0.name,sg1.name = sg1.name,Kfolds = Kfolds)
  )
}


#' ForestSearch Tenfold Cross-Validation Simulation
#'
#' Runs repeated K-fold (default 10-fold) cross-validation simulations for ForestSearch and summarizes performance.
#'
#' @param fs.est ForestSearch results object.
#' @param sims Integer. Number of simulation repetitions.
#' @param Kfolds Integer. Number of folds (default: 10).
#' @param details Logical. Print details during execution (default: TRUE).
#' @param parallel_args List. Parallelization arguments (plan, workers, show_message).
#'
#' @return List with median and full results for sensitivity and subgroup finding metrics.
#'
#' @importFrom future plan
#' @importFrom data.table data.table copy
#' @importFrom foreach foreach
#' @export

forestsearch_tenfold <- function(fs.est, sims, Kfolds = 10, details = TRUE,
                                 parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE)){

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
# Note: setup_parallel_SGcons is re-purposed from subgroup_consistency
setup_parallel_SGcons(parallel_args)

t.start<-proc.time()[3]
## Will summarize cross-validated metrics
## Initialize here (will be concatenated below)
sens_out <- NULL
find_out <- NULL
fs_args <- fs.est$args_call_all
get_names <- with(fs_args,c(confounders.name, outcome.name, event.name, id.name, treat.name))
## Extract necessary data
dfa <- fs.est$df.est[,c(get_names,"treat.recommend")]
## re-name treat.recommend to original to compare with k-folds
names(dfa)[names(dfa) == "treat.recommend"] <- "treat.recommend.original"
dfnew <- as.data.frame(dfa)
# turn off details in fs_args and reset parallel_args to turn off parallel
fs_args$parallel_args <- list()
fs_args$details <- FALSE
fs_args$plot.sg <- FALSE

cv_args <- fs_args

est.scale <- fs_args$est.scale
# Pre-allocate lists to collect results
sens_out_list <- vector("list", sims)
find_out_list <- vector("list", sims)

outcome.name <- fs_args$outcome.name
event.name <- fs_args$event.name
treat.name <- fs_args$treat.name
id.name <- fs_args$id.name
potentialOutcome.name <- fs_args$potentialOutcome.name

simulation_results <- lapply(seq_len(sims), function(ksim) {
  t.startk <- proc.time()[3]
  df_scrambled <- copy(dfnew)
  set.seed(8316951 + 1000 * ksim)
  id_sample <- sample(seq_len(nrow(df_scrambled)), replace = FALSE)
  df_scrambled <- df_scrambled[id_sample, ]
  folds <- cut(seq_len(nrow(df_scrambled)), breaks = Kfolds, labels = FALSE)
  if (details && ksim <= 3) cat("Range of unique left-out fold sample sizes (first 3 simulations)", range(unique(table(folds))), "\n")

  resCV <- foreach::foreach(
    cv_index = seq_len(Kfolds),
    .options.future = list(seed = TRUE, add = c("df_scrambled", "folds", "fs.est")),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    testIndexes <- which(folds == cv_index, arr.ind = TRUE)
    x.test <- df_scrambled[testIndexes, ]
    x.train <- df_scrambled[-testIndexes, ]
    cv_args$df.analysis <- x.train
    fs.train <- suppressWarnings(try(do.call(forestsearch, cv_args), TRUE))
    if (!inherits(fs.train, "try-error")) {
      if (!is.null(fs.train$sg.harm)) {
        if (details && ksim <= 3) {
          cat("Simulation, Fold =", c(ksim, cv_index), "\n")
          print(fs.train$sg.harm)
        }
        sg1 <- fs.train$sg.harm[1]
        sg2 <- fs.train$sg.harm[2]
        df.test <- get_dfpred(df.predict = x.test, sg.harm = fs.train$sg.harm, version = 2)
        df.test$cvindex <- rep(cv_index, nrow(df.test))
        df.test$sg1 <- rep(sg1, nrow(df.test))
        df.test$sg2 <- rep(sg2, nrow(df.test))
        dfres <- data.table(df.test[, c(id.name, outcome.name, event.name, treat.name, "treat.recommend", "treat.recommend.original", "cvindex", "sg1", "sg2")])
      } else {
        if (details && ksim <= 3) cat("Subgroup not found (return ITT), Fold =", c(ksim, cv_index), "\n")
        df.test <- x.test
        df.test$cvindex <- rep(cv_index, nrow(df.test))
        df.test$sg1 <- rep(NA, nrow(df.test))
        df.test$sg2 <- rep(NA, nrow(df.test))
        df.test$treat.recommend <- rep(1.0, nrow(df.test))
        dfres <- data.table(df.test[, c(id.name, outcome.name, event.name, treat.name, "treat.recommend", "treat.recommend.original", "cvindex", "sg1", "sg2")])
      }
    }
    if (inherits(fs.train, "try-error")) {
      df.test <- x.test
      df.test$cvindex <- rep(cv_index, nrow(df.test))
      df.test$sg1 <- rep(NA, nrow(df.test))
      df.test$sg2 <- rep(NA, nrow(df.test))
      df.test$treat.recommend <- rep(1.0, nrow(df.test))
      dfres <- data.table(df.test[, c(id.name, outcome.name, event.name, treat.name, "treat.recommend", "treat.recommend.original", "cvindex", "sg1", "sg2")])
    }
    return(dfres)
  }

  if (length(unique(resCV$id)) != nrow(df_scrambled)) stop("K-fold cross-validation sample size does not equal full analysis dataset")
  chk1 <- range(unique(table(folds)))
  chk2 <- range(unique(table(resCV$cvindex)))
  if (!all(chk1 == chk2)) stop("Mismatch on range of fold sample sizes")

  res <- list(
    resCV = resCV, Kfolds = Kfolds, confounders.name = confounders.name, cv_args = cv_args,
    outcome.name = outcome.name, event.name = event.name, id.name = id.name,
    treat.name = treat.name, sg_analysis = fs.est$sg.harm, sg0.name = "Not recommend", sg1.name = "Recommend"
  )
  out <- forestsearch_KfoldOut(res = res, outall = FALSE, details = FALSE)

  if (details & (ksim %in% c(1:5, 10, 20, 50, 100, 200, 300))) {
    cat("Kfold iteration done=", ksim, "\n")
  }

  list(
    sens_metrics_original = out$sens_metrics_original,
    find_metrics = out$find_metrics
  )
})

# Combine results
sens_out <- do.call(rbind, lapply(simulation_results, `[[`, "sens_metrics_original"))
find_out <- do.call(rbind, lapply(simulation_results, `[[`, "find_metrics"))

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
if(details){
cat("Minutes for Cross-validation: Kfolds, sims, minutes=",c(Kfolds,sims,t.min),"\n")
cat("Projected hours per 100 sims",c((t.min/60)*(100/sims)),"\n")
}
sens_summary <- apply(sens_out,2,median,na.rm=TRUE)
find_summary <- apply(find_out,2,median,na.rm=TRUE)
return(list(sens_summary=sens_summary,find_summary=find_summary,sens_out=sens_out,find_out=find_out))
}


