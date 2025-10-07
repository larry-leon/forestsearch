

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


#library(forestsearch)


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
                   parallel_args = list(plan="multisession", workers = 6, show_message = TRUE)
                   )
})


res_tabs <- sg_tables(fs, ndecimals = 3)

res_tabs$tab_estimates

res_tabs$sg10_out



#plan(sequential)

temp <- get_FSdata(df = dfa, confounders.name = confounders.name, cut_type = "default", use_grf = TRUE,
                   outcome.name = "tte", event.name = "event", use_lasso = FALSE, conf_force=c("age <= 65","meno == 0"))



# registerDoRNG() sets up reproducible random number generation for all subsequent foreach parallel operations in the session.
# You do not need to call it before every foreach or dofuture loop.
# If you restart your R session or change the parallel backend, you should call it again.
# Call once per session/backend setup.
# Not needed before every parallel loop.


NB <- 5
t.start <- proc.time()[3]
# Bootstrap bias-correction
# fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs,nb_boots = NB, show_three = TRUE, details = TRUE, reset_paralle_fs = FALSE, boot_workers =3)
fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = TRUE, details = TRUE, reset_parallel_fs = TRUE,
                                         parallel_args = list(plan = "callr", workers = 6, show_message = TRUE) )
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




