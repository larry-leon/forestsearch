# ============================================================================
# Demonstration: GRF subgroup LABELS (sg.harm.id) can be grossly wrong at
# grf_depth >= 2, even though GRF's own membership (treat.recommend) is correct.
#
# Construct a DGM whose true harm region is a UNION of two conjunctions reached
# by RIGHT-turns, forcing the policy tree to split twice into two action-1
# leaves -- the case both label builders mishandle:
#   * survival find_leaf_split(): returns ONE split above ONE leaf
#   * GLM .build_sg_harm_id():    emits ALL node splits as "<="
#
# We then show that the SUBGROUP YOU GET BY APPLYING THE PRINTED LABEL has a
# materially different treatment effect than the subgroup GRF actually found --
# i.e. the printed definition would mislead.
# ============================================================================


library(forestsearch)
library(survival)
library(data.table)
library(gt)
library(doFuture)
library(foreach)

# suppressMessages({library(survival); library(data.table); library(grf)
#   library(policytree); library(glmnet)})
# R <- "fs_latest/R"
# invisible(lapply(list.files(R, pattern="[.]R$", full.names=TRUE),
#                  function(f) suppressWarnings(source(f))))

# True harm region: (x1 > 0.0 & x2 > 0.0)  OR  (x1 <= 0.0 & x3 > 0.0)
#   -> both disjuncts require a RIGHT turn (">"), and the region spans 2 leaves.
in_true_harm <- function(X) (X[,1] > 0 & X[,2] > 0) | (X[,1] <= 0 & X[,3] > 0)

hr_in <- function(df, mask, tvar, ...) {
  # observed HR of treatment within mask (survival) -- quick Cox
  dd <- df[mask, , drop=FALSE]
  if (sum(dd[[tvar]]==1) < 5 || sum(dd[[tvar]]==0) < 5) return(NA_real_)
  cph <- tryCatch(coxph(Surv(time, status) ~ trt, data=dd), error=function(e) NULL)
  if (is.null(cph)) return(NA_real_)
  exp(unname(coef(cph)[1]))
}

## ---------------- SURVIVAL ----------------
set.seed(101); n <- 3000
X <- matrix(runif(n*5,-1,1),n,5); colnames(X)<-paste0("x",1:5)
W <- rbinom(n,1,.5); inS <- in_true_harm(X)
tau <- ifelse(inS, 1.3, -0.3)              # strong harm in S, benefit outside
eta <- 0.15*X[,1] + W*tau
te <- rexp(n, exp(eta)); tc <- runif(n,0,1.6)
dS <- data.frame(time=pmin(te,tc), status=as.integer(te<=tc), trt=W,
                 as.data.frame(X), id=seq_len(n)); covs<-paste0("x",1:5)

args <- forestsearch:::.build_grf_survival_args(data=dS, confounders.name=covs, outcome.name="time",
  event.name="status", id.name="id", treat.name="trt", frac.tau=0.80, n.min=60,
  dmin.grf=0.0, is.RCT=TRUE, grf_depth=2, seedit=8316951, return_selected_cuts_only=TRUE)
res <- do.call(grf.subg.harm.survival, args)

# Build the get_dfpred() argument from the structured, path-based definition:
#   * conjunction -> brace-label vector ("{var op val}")
#   * union       -> the disjunctive definition string "(a & b) | (c & d)"
# get_dfpred() handles both (it routes a "|"-containing string to the
# disjunctive evaluator).  Brace-wrapping the union string would break parsing.
def_arg <- function(grf_res) {
  d <- grf_res$sg_def
  if (!is.null(d$labels)) d$labels else d$definition
}

printed_label <- res$sg_def$definition
grf_member  <- res$data$treat.recommend == 0                 # GRF's TRUE subgroup
label_member<- get_dfpred(dS, def_arg(res), version=2)$treat.recommend == 0

cat("================= SURVIVAL =================\n")
cat("PRINTED subgroup definition (sg.harm.id):", printed_label, "\n")
cat(sprintf("GRF's actual subgroup size:        %d\n", sum(grf_member)))
cat(sprintf("Size if you APPLY the printed rule: %d\n", sum(label_member)))
cat(sprintf("Subjects misclassified by label:   %d (%.0f%%)\n",
    sum(grf_member != label_member), 100*mean(grf_member != label_member)))
cat(sprintf("HR in GRF's actual subgroup:        %.2f   <- what GRF found (harm)\n",
    hr_in(dS, grf_member, "trt")))
cat(sprintf("HR in the PRINTED-label subgroup:   %.2f   <- what a reader would compute\n",
    hr_in(dS, label_member, "trt")))
cat(sprintf("HR true-harm region (oracle):       %.2f\n", hr_in(dS, inS, "trt")))

## ---------------- GLM (binary / OR) ----------------
set.seed(101); n <- 3000
X <- matrix(runif(n*5,-1,1),n,5); colnames(X)<-paste0("x",1:5)
W <- rbinom(n,1,.5); inS <- in_true_harm(X)
tau <- ifelse(inS, 1.8, -0.6)
Yb <- rbinom(n,1,plogis(-0.2+0.1*X[,1]+W*tau))
dG <- data.frame(Y=Yb, trt=W, as.data.frame(X), id=seq_len(n))

ga <- forestsearch:::.build_grf_glm_args(data=dG, confounders.name=covs, outcome.name="Y", treat.name="trt",
  id.name="id", outcome_type="binary", n.min=60, dmin.grf=0.0, is.RCT=TRUE, grf_depth=2,
  seedit=8316951, return_selected_cuts_only=TRUE, adverse_outcome=FALSE)
rg <- do.call(grf.subg.harm.glm, ga)
plab <- rg$sg_def$definition
gm  <- rg$data$treat.recommend == 0
lm_ <- get_dfpred(dG, def_arg(rg), version=2)$treat.recommend == 0
or_in <- function(df, mask) { dd<-df[mask,,drop=FALSE]
  if (length(unique(dd$trt))<2) return(NA_real_)
  fit<-tryCatch(glm(Y~trt,binomial,dd),error=function(e)NULL); if(is.null(fit))return(NA_real_)
  exp(unname(coef(fit)[2])) }
cat("\n================= GLM (binary, OR) =================\n")
cat("PRINTED subgroup definition (sg.harm.id):", plab, "\n")
cat(sprintf("GRF's actual subgroup size:        %d\n", sum(gm)))
cat(sprintf("Size if you APPLY the printed rule: %d\n", sum(lm_)))
cat(sprintf("Subjects misclassified by label:   %d (%.0f%%)\n",
    sum(gm != lm_), 100*mean(gm != lm_)))
cat(sprintf("OR in GRF's actual subgroup:        %.2f   <- what GRF found (harm)\n", or_in(dG, gm)))
cat(sprintf("OR in the PRINTED-label subgroup:   %.2f   <- what a reader would compute\n", or_in(dG, lm_)))
cat(sprintf("OR true-harm region (oracle):       %.2f\n", or_in(dG, inS)))
