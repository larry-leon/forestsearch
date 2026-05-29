## README ##

The scripts implement DINA, an estimator of heterogeneous treatment effects defined by the difference in natural parameters. The method is particularly suitable for outcomes following models from the exponential family or the Cox proportional hazards model.

Below is an overview of the scripts. The simulation parameters used in the paper are provided in a separate .xlsx file.

DINA.R
  This script implements the estimator of DINA, which makes use of existing machine learning tools and admits desirable statistical properties.


coef.DINA.R
  This script extracts estimated coefficients from a fitted DINA object.

confint.DINA.R
  This script implements the confidence interval estimation for the DINA coefficients using bootstrap.

predict.DINA.R
  This script generates predictions of treatment effect for new data point based on a fitted DINA model.

summary.DINA.R
  This script summarizes a fitted DINA object, including estimated coefficients, and confidence intervals.

helper.R
  This script contains utility functions used by other scripts.


Example:
Step 1. Source all the R scripts in this directory:

   ```r
   source("helper.R")
   source("DINA.R")
   source("coef.DINA.R")
   source("confint.DINA.R")
   source("predict.DINA.R")
   source("summary.DINA.R")
   ```

Step 2. Fit a DINA model:

   ```r
   fit <- DINA = function(data, family = "poisson")
   ```

Step 3. Extract results:

   ```r
   summary(fit)
   ```

