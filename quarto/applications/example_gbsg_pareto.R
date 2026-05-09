## Simple example: forestsearch on GBSG with selection_rule = "pareto"
##
## Demonstrates the new selection_rule argument on the GBSG breast-cancer
## dataset.  Mirrors the structure of analysis_gbsg_survival_maxSG.qmd
## but switches sg_focus to "hrMaxSG" so that selection_rule is active.

library(forestsearch)
library(survival)

# ---- Data ----

df.analysis <- gbsg
df.analysis <- within(df.analysis, {
  id          <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1L, 0L)
})

confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")

# ---- Run ----

fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  sg_focus               = "hrMaxSG",
  selection_rule         = "pareto",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

# ---- Inspect ----

cat("\nIdentified subgroup:", paste(fs$sg.harm, collapse = " & "), "\n")
cat("Selection rule used:", fs$args_call_all$selection_rule, "\n\n")

pareto_frontier_table(fs)


fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  sg_focus               = "hrMaxSG",
  selection_rule         = "neighborhood",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)



fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  sg_focus               = "hrMaxSG",
  selection_rule         = "both",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)


fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  use_grf = FALSE,
  sg_focus               = "hrMaxSG",
  selection_rule         = "both",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)


fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  use_grf = FALSE,
  sg_focus               = "hrMaxSG",
  selection_rule         = "neighborhood",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)


# How does it fail on sg_focus = hr with pareto

fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  sg_focus               = "hr",
  selection_rule         = "both",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)



fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  use_grf = FALSE,
  sg_focus               = "hr",
  selection_rule         = "neighborhood",        # <-- new argument
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)


fs <- forestsearch(
  df.analysis            = df.analysis,
  confounders.name       = confounders.name,
  outcome.name           = "time_months",
  treat.name             = "hormon",
  event.name             = "status",
  id.name                = "id",
  hr.threshold           = 1.25,
  hr.consistency         = 1.0,
  pconsistency.threshold = 0.90,
  use_grf = FALSE,
  sg_focus               = "maxSG",
  selection_rule         = "neighborhood",        # <-- new argument
  max_subgroups_search = 20,
  conf.cont_jcuts        = list(er = 10, pgr = 10),
  maxk                   = 2,
  n.min                  = 60,
  d0.min                 = 10,
  d1.min                 = 10,
  fs.splits              = 1000,
  details                = TRUE
)

pareto_frontier_table(fs)


