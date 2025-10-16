
rm(list=ls())

suppressMessages(library(weightedSurv))
suppressMessages(library(forestsearch))



library(survival)
df_gbsg <- gbsg
df_gbsg$tte <- df_gbsg$rfstime / 30.4375
df_gbsg$event <- df_gbsg$status
df_gbsg$treat <- df_gbsg$hormon
df_gbsg$grade3 <- ifelse(df_gbsg$grade == "3", 1, 0)
# create id name
df_gbsg$id <- seq_len(nrow(df_gbsg))
# If missing, then created automatically

tte.name <- "tte"
event.name <- "event"
treat.name <- "treat"
id.name <- "id"
arms <- c("treat", "control")

# Baseline factors from which candidate subgroups are formed
confounders.name<-c("age","meno","size","grade3","nodes","pgr","er")

library(doFuture)
library(doRNG)
registerDoFuture()
registerDoRNG()


fs <- forestsearch(df_gbsg,  confounders.name=confounders.name,
                                outcome.name = "tte", treat.name = "treat", event.name = "event", id.name = "id",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.80,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = c("age <= 65"),
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 10, d1.min = 10,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="callr", workers = 100, show_message = TRUE)
)


output_dir <- "Vignettes/results/"
save_results <- dir.exists(output_dir)

NB <- 30

t.start <- proc.time()[3]

fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = TRUE, details = TRUE)

print(names(fs_bc))

print(fs_bc$FSsg_tab)

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes (total) for bootstrap (boots,mins)",c(NB,t.min),"\n")
cat("Projected minutes for 1000",c(t.min*(1000/NB)),"\n")


# Subgroup       n             n1            events        m1     m0     RMST  HR (95% CI)         HR*
#   res_0 "Questionable" "82 (12.0%)"  "26 (31.7%)"  "45 (54.9%)"  "22.9" "43.7" "-14" "1.95 (1.04, 3.67)" "1.5 (0.39,5.81)"
# res_1 "Recommend"    "604 (88.0%)" "220 (36.4%)" "254 (42.1%)" "66.7" "52.6" "9.3" "0.61 (0.47, 0.80)" "0.64 (0.05,7.55)"


if (save_results) {
  filename <- file.path(output_dir,
                        paste0("bootstrap_results_B=",
                               format(NB),
                               ".RData"))
  save(fs_bc, fs, file = filename)
  cat("\nResults saved to:", filename, "\n")
}



