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
                                parallel_args = list(plan="multisession", workers = 12, show_message = TRUE)
)


NB <- 30

fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = TRUE, details = TRUE)

