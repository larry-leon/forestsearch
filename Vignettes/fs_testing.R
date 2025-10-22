
library(weightedsurv)

# library(forestsearch)

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

confounders.name<-c("age","meno","size","grade3","nodes","pgr","er")

library(doFuture)
library(doRNG)
registerDoFuture()
registerDoRNG()


# Example with pre-specified subgroup via conf_force argument
# Note that er and pgr could be important "biomarkers" and we pre-specify "low expression cuts"
# up to the 1st quartiles which are 8 (er) and 7 (pgr);
# Note that by include "low expression" values the complementary levels are automatically considered.
# For example "er <= 0" means "er > 0" evaluated (etc)

system.time({fs <- forestsearch(df_gbsg,  confounders.name=confounders.name,
                   outcome.name = "tte", treat.name = "treat", event.name = "event", id.name = "id",
                   hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                   sg_focus = "hrMaxSG",
                   showten_subgroups = FALSE, details=TRUE,
                   conf_force = c("age <= 65", "er <= 0", "er <= 1", "er <= 2"," er <= 5", "er <= 8",
                                  "pgr <= 0", "pgr <= 1", "pgr <= 2", "pgr <= 5", "pgr <= 7"),
                   max_subgroups_search = 12,
                   cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                   maxk = 2, n.min = 60, d0.min = 12, d1.min = 12,
                   plot.sg = TRUE, by.risk = 12,
                   parallel_args = list(plan="callr", workers = 12, show_message = TRUE)
                   )
})




# Re-produce analysis in paper
system.time({fs_paper <- forestsearch(df_gbsg,  confounders.name=confounders.name,
                                outcome.name = "tte", treat.name = "treat", event.name = "event", id.name = "id",
                                hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.90,
                                sg_focus = "hrMaxSG",
                                showten_subgroups = FALSE, details=TRUE,
                                conf_force = NULL,
                                cut_type = "default", use_grf = TRUE, plot.grf = TRUE, use_lasso = TRUE,
                                maxk = 2, n.min = 60, d0.min = 10, d1.min = 10,
                                plot.sg = TRUE, by.risk = 12,
                                parallel_args = list(plan="callr", workers = 12, show_message = TRUE)
)
})


# reset workeres
plan("sequential")



library(patchwork)

options(warn = -1)

output_dir <- "results/"
save_results <- dir.exists(output_dir)

NB <- 20000

t.start <- proc.time()[3]

fs_bc <- forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = NB, show_three = FALSE, details = TRUE, create_summary = TRUE, create_plots = TRUE)

plan("sequential")

t.now<-proc.time()[3]
t.min<-(t.now-t.start)/60
cat("Minutes (total) for bootstrap (boots,mins)",c(NB,t.min),"\n")
cat("Projected minutes for 1000",c(t.min*(1000/NB)),"\n")

if (save_results) {
  filename <- file.path(output_dir,
                        paste0("gbsg_results_B=",
                               format(NB),
                               ".RData"))
  save(fs_bc, fs, file = filename)
  cat("\nResults saved to:", filename, "\n")
}


# SG estimates
fs_bc$summary$table
# View the diagnostics table
fs_bc$summary$diagnostics_table_gt
#View the formatted timing table
fs_bc$summary$timing$time_table_gt


sg_summary <- fs_bc$summary$subgroup_summary

# Format as gt tables
Bsg_tables <- format_subgroup_summary_tables(sg_summary, nb_boots = NB)

# View tables
Bsg_tables$basic_stats

# too many
Bsg_tables$agreement

# Also a lot
Bsg_tables$factor_freq

Bsg_tables$consistency_dist

Bsg_tables$original_agreement

Bsg_tables$factor_presence


# Print text summary
cat("Subgroups identified in", sg_summary$pct_found, "% of bootstraps\n")
cat("Most common subgroup:", sg_summary$agreement$Subgroup[1], "\n")
cat("  Appeared in", sg_summary$agreement$Percent[1], "% of successful iterations\n")



