#' DINA
#'
#' This function implements DINA --- an estimator of heterogeneous treatment effects for general responses. DINA makes use of existing machine learning tools and enjoys desirable statistical properties.
#'
#' @param data list of data. \code{data$X} denotes the covariate matrix: each row stands for an observation, each column stands for a covariate. \code{data$y} denotes the vector of responses: for \code{family = "cox"}, \code{data$y} is a data.frame with two columns: "time" for the vector of survival times, "status" for the vector of censored status. \code{Status = 1} stands for not censored, and \code{status = 0} stands for being censored. \code{data$W} denotes the vector of treatment assignment: \code{W = 1} stands for treatment, and \code{W = 0} stands for control.
#' @param family a character string describing the error distribution and link function to be used in the model. \code{Family} can be one of the following options: "gaussian", "binomial", "poisson", "cox". Default is "gaussian".
#' #' @param prop method for estimating the propensity score. \code{Prop} can be a character string naming the estimator: "logistic regression", "gradient boosting". \code{Prop} can also be a function and \code{prop(X)} should output the estimated propensity scores for each row of the covariate matrix \code{X}. Default is "logistic regression".
#' @param paramsW list of parameters for propensity score estimation.
#' @param baseline method for estimating the natural parameter functions in the treatment and control groups. \code{Baseline} can be a character string naming the estimator: "regression", "gradient boosting". \code{Baseline} can also be a list with two members: \code{eta1} denotes the estimator of the natural parameter function in the treatment group,   \code{eta0} denotes the estimator of the natural parameter function in the control group. The functions should take in a covariate matrix \code{X} and output the estimated quantities for each row of the covariate matrix \code{X}. Default is "regression".
#' @param params list of parameters for the natural parameter function estimation.
#' @param cross.fitting logical value indicating whether to use cross-fitting. If \code{TRUE}, the estimation of the nuisance functions and causal parameters are based on independent observations. Default is \code{TRUE}.
#' @param k number of folds for cross-fitting if \code{cross.fitting = TRUE}. Default is 2.
#' @param variance.method method for estimating the variance matrix. \code{Method} can be a character string naming the estimator: "CLT", "bootstrap". Default is "CLT".
#' @param B the number of bootstrap resamples. This should be a single positive integer. Default is 100.
#'
#' @return an object of class \code{DINA}.
#'  \itemize{
#'  \item{estimator}{vector of estimated coefficients.}
#'  \item{Sigma}{variance matrix of the coefficients.}
#'  \item{params}{list of hyper-parameters used.}
#' }
#'
#' @export
DINA = function(data, family = "gaussian", prop = "logistic regression", paramsW = NULL, baseline = "gradient boosting", params = list(), cross.fitting = TRUE, k = 2, variance.method = "CLT", B = 100){
  # preprocess
  if(variance.method == "CLT"){B = NULL}
  if(variance.method == "bootstrap" && B < 2){
    stop("To estimate the variance matrix by bootstrap, we need more than one bootstrap resample.")
  }
  n = dim(data$X)[1]; d = dim(data$X)[2]
  if(family != "cox"){data$Y = as.matrix(data$Y)}

  # hyper-parameters for nuisance function estimation
  if(family == "gaussian"){
    params$objective = "reg:squarederror"
    params$eval_metric = "rmse"
    link.function = function(x){x}; var.function = function(x){x * 0 + 1}
  }
  if(family == "binomial"){
    params$objective = "binary:logistic"
    params$eval_metric = "logloss"
    link.function = function(x){log(x/(1-x))}; var.function = function(x){x*(1-x)}
  }
  if(family == "poisson"){
    params$objective = "count:poisson"
    params$eval_metric = "poisson-nloglik"
    link.function = function(x){log(x)}; var.function = function(x){x}
  }
  if(family == "cox"){
    params$objective = "survival:cox"
    params$eval_metric = "cox-nloglik"
    link.function = NULL; var.function = NULL
  }
  DINA.terms = list(); DINA.terms$family = family; DINA.terms$prop = prop; DINA.terms$baseline = baseline; DINA.terms$cross.fitting.fold = k; DINA.terms$variance.method = variance.method; DINA.terms$variance.method.B = B

  # no cross-fitting
  if(!cross.fitting){
    result = DINA.fit(X = data$X, Y = data$Y, W = data$W, family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, k = k, link.function = link.function, var.function = var.function) # estimation
    if(variance.method == "CLT"){result$terms = DINA.terms; class(result) = "DINA"; return(result)}
    if(variance.method == "bootstrap"){
      estimator = matrix(0, nrow = B, ncol = d+1)
      for(b in 1:B){
        index.resample = sample(n, n, replace = T)
        estimator[b,] = DINA.fit(X = data$X[index.resample,,drop = F], Y = data$Y[index.resample,,drop = F], W = data$W[index.resample], family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, k = k, link.function = link.function, var.function = var.function)$estimator
      }
      result$Sigma = cov(estimator)
      result$terms = DINA.terms
      class(result) = "DINA"; return(result)
    }
  }

  # cross-fitting
  split.fold = split.uniform(n = n, k = k, replace = F) # sample splitting
  result = DINA.fit(X = data$X, Y = data$Y, W = data$W, family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, k = k, index.resample = split.fold$index.resample, foldIndex = split.fold$foldIndex, link.function = link.function, var.function = var.function)
  if(variance.method == "CLT"){
    result$terms = DINA.terms; class(result) = "DINA"
    return(result)
  }
  if(variance.method == "bootstrap"){
    estimator = matrix(0, nrow = B, ncol = d+1)
    for(b in 1:B){
      split.fold = split.uniform(n = n, k = k, replace = T)
      estimator[b,] = DINA.fit(X = data$X, Y = data$Y, W = data$W, family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, k = k, index.resample = split.fold$index.resample, foldIndex = split.fold$foldIndex, link.function = link.function, var.function = var.function)$estimator # estimation
    }
    result$Sigma = cov(estimator) # use the covariance matrix estimated by bootstrap
    result$terms = DINA.terms
    class(result) = "DINA"; return(result)
  }
}


#' DINA.bagging
#'
#' This function implements DINA averaged over bootstrap resamples. The estimator DINA is executed over multiple bootstrap resamples, and the final estimator will be the average of the multiple realizations.
#'
#' @param data list of data. \code{data$X} denotes the covariate matrix: each row stands for an observation, each column stands for a covariate. \code{data$y} denotes the vector of responses: for \code{family = "cox"}, \code{data$y} is a data.frame with two columns: "time" for the vector of survival times, "status" for the vector of censored status. \code{Status = 1} stands for not censored, and \code{status = 0} stands for being censored. \code{data$W} denotes the vector of treatment assignment: \code{W = 1} stands for treatment, and \code{W = 0} stands for control.
#' @param family a character string describing the error distribution and link function to be used in the model. \code{Family} can be one of the following options: "gaussian", "binomial", "poisson", "cox". Default is "gaussian".
#' #' @param prop method for estimating the propensity score. \code{Prop} can be a character string naming the estimator: "logistic regression", "gradient boosting". \code{Prop} can also be a function and \code{prop(X)} should output the estimated propensity scores for each row of the covariate matrix \code{X}. Default is "logistic regression".
#' @param paramsW list of parameters for propensity score estimation.
#' @param baseline method for estimating the natural parameter functions in the treatment and control groups. \code{Baseline} can be a character string naming the estimator: "regression", "gradient boosting". \code{Baseline} can also be a list with two members: \code{eta1} denotes the estimator of the natural parameter function in the treatment group,   \code{eta0} denotes the estimator of the natural parameter function in the control group. The functions should take in a covariate matrix \code{X} and output the estimated quantities for each row of the covariate matrix \code{X}. Default is "regression".
#' @param params list of parameters for the natural parameter function estimation.
#' @param k number of folds for cross-fitting. Default is 2.
#' @param B the number of bootstrap resamples. This should be a single positive integer no smaller than 100. Default is 100.
#' @param variance.method method for estimating the variance matrix. \code{Method} can be a character string naming the estimator: "infinitesimal jackknife", "jackknife". Default is "infinitesimal jackknife".
#' @param variance.bias.correction logical value indicating whether to execute bias correction. If \code{TRUE}, an estimator of the bias will be subtracted from the original estimator. Default is \code{TRUE}.
#'
#' @return an object of class \code{DINA}.
#'  \itemize{
#'  \item{estimator}{vector of estimated coefficients.}
#'  \item{Sigma}{variance matrix of the coefficients.}
#' }
#'
#' @export
DINA.bagging = function(data, family = "gaussian", prop = "logistic regression", paramsW = NULL, baseline = "gradient boosting", params = list(), k = 2,  B = 100, variance.method = "infinitesimal jackknife", variance.bias.correction = T){
  # preprocess
  if(B < 100){stop("We recommend using at least 100 base DINA estimators to compute the bagging estimator.")}
  n = dim(data$X)[1]; d = dim(data$X)[2]
  if(family != "cox"){data$Y = as.matrix(data$Y)}

  # hyper-parameters for nuisance function estimation
  if(family == "gaussian"){
    params$objective = "reg:squarederror"
    params$eval_metric = "rmse"
    link.function = function(x){x}; var.function = function(x){x * 0 + 1}
  }
  if(family == "binomial"){
    params$objective = "binary:logistic"
    params$eval_metric = "logloss"
    link.function = function(x){log(x/(1-x))}; var.function = function(x){x*(1-x)}
  }
  if(family == "poisson"){
    params$objective = "count:poisson"
    params$eval_metric = "poisson-nloglik"
    link.function = function(x){log(x)}; var.function = function(x){x}
  }
  if(family == "cox"){
    params$objective = "survival:cox"
    params$eval_metric = "cox-nloglik"
  }

  bootstrap.membership = matrix(0, nrow = B, ncol = n)
  estimator = matrix(0, nrow = B, ncol = d+1)
  result = list()
  for(b in 1:B){
    # standard cross-fitting estimator based on one bootstrap resample
    split.fold = split.weight(n = n, k = k) # sample splitting of a bootstrap resample with duplicates
    result = DINA.fit(X = data$X, Y = data$Y, W = data$W, family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, k = k, index.resample = split.fold$index.resample, foldIndex = split.fold$foldIndex, link.function = link.function, var.function = var.function) # estimation; here!!!
    estimator[b,] = result$estimator
    bootstrap.membership[b,] = split.fold$weight
  }
  result$estimator = apply(estimator, 2, mean)

  # variance matrix
  result$Sigma = DINA.bagging.inference(bagging.estimator = estimator, bootstrap.membership = bootstrap.membership, method = variance.method, bias.correction = variance.bias.correction)
  result$Sigma = psd.projection(result$Sigma)
  DINA.terms = list(); DINA.terms$family = family; DINA.terms$prop = prop; DINA.terms$baseline = baseline; DINA.terms$cross.fitting.fold = k; DINA.terms$bootstrap.sample.size = B; DINA.terms$variance.method = variance.method; DINA.terms$variance.bias.correction = variance.bias.correction
  result$terms = DINA.terms
  class(result) = "DINA"; return(result)
}


#' DINA.fit
#'
#' This function is a helper for the function \code{DINA}.
#'
#' @param X covariate matrix. Each row stands for an observation, each column stands for a covariate.
#' @param Y vector of response. For \code{family = "cox"}, \code{Y} is a data.frame with two columns: "time" for the vector of survival times, "status" for the vector of censored status. \code{Status = 1} stands for not censored, and \code{status = 0} stands for being censored.
#' @param W vector of treatment assignment. \code{W = 1} stands for treatment, and \code{W = 0} stands for control.
#' @param family a character string describing the error distribution and link function to be used in the model. \code{Family} can be one of the following options: "gaussian", "binomial", "poisson", "cox". Default is "gaussian".
#' @param prop method for estimating the propensity score. \code{Prop} can be a character string naming the estimator: "logistic regression", "gradient boosting". \code{Prop} can also be a function and \code{prop(X)} should output the estimated propensity scores for each row of the covariate matrix \code{X}. Default is "logistic regression".
#' @param paramsW list of parameters for propensity score estimation.
#' @param baseline method for estimating the natural parameter functions in the treatment and control groups. \code{Baseline} can be a character string naming the estimator: "regression", "gradient boosting". \code{Baseline} can also be a list with two members: \code{eta1} denotes the estimator of the natural parameter function in the treatment group,   \code{eta0} denotes the estimator of the natural parameter function in the control group. The functions should take in a covariate matrix \code{X} and output the estimated quantities for each row of the covariate matrix \code{X}. Default is "regression".
#' @param params list of parameters for the natural parameter function estimation.
#' @param link.function link function from the conditional mean to the natural parameter.
#' @param var.function variance function from the conditional mean to the variance function.
#' @param k number of folds for cross-fitting. Default is 2.
#' @param index.resample vector of subject indices in the permuted sample. We permute the original dataset and each subject will appear exactly once in the permuted sample. The permutation is done to break inherent ordering in the original dataset (if there is any) and makes sure the resulted folds are exchangeable.
#' @param foldIndex list of fold memberships. The j-th member is a vector containing all subjects of the bootstrap resample assigned to the j-th fold. If \code{foldIndex = NULL}, no sample splitting is done and the estimation of nuisance functions and causal parts are both based on the full dataset.
#' @return an object of class \code{DINA}.
#'  \itemize{
#'  \item{estimator}{vector of estimated coefficients.}
#'  \item{Sigma}{variance matrix of the coefficients.}
#' }
#'
#' @export
DINA.fit = function(X, Y, W, family, prop = "logistic regression", paramsW = NULL, baseline = "gradient boosting", params = list(), k = 2, index.resample = NULL, foldIndex = NULL, link.function = NA, var.function = NA){
  # preprocess
  if(family != "cox"){Y = as.matrix(Y)}
  n = dim(X)[1]; d = dim(X)[2]
  if(is.character(prop)){
    if(prop == "gradient boosting"){
      if(is.null(paramsW)){
        paramsW = list()
        paramsW$nrounds = 100; paramsW$params = list(eta = 0.1)
      }
    }
  }
  if(is.character(baseline)){
    if(baseline == "gradient boosting"){
      if(is.null(params$params)){params$params = list(eta = 0.1)}
      if(is.null(params$nrounds)){params$nrounds = 100}
    }
  }
  # record
  result = list()
  estimator = Sigma1 = Sigma2 = list()

  # no sample splitting
  if(is.null(foldIndex)){
    nuisance.function = DINAEstimator.step1(X, Y, W, family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, link.function = link.function, var.function = var.function) # estimate the nuisance functions
    causal.estimator = DINAEstimator.step2(X = X, Y = Y, W = W, family = family, nuisance.function = nuisance.function) # estimate the causal part
    result$estimator = causal.estimator$estimator
    result$Sigma = solve(causal.estimator$Sigma1) %*% causal.estimator$Sigma2 %*% solve(causal.estimator$Sigma1)/n # compute the variance matrix of the coefficients
    class(result) = "DINA"; return(result)
  }

  # cross-fitting
  for(i in 1:k){
    nuisance.function = DINAEstimator.step1(X[index.resample[-foldIndex[[i]]],], Y[index.resample[-foldIndex[[i]]],], W[index.resample[-foldIndex[[i]]]], family = family, prop = prop, paramsW = paramsW, baseline = baseline, params = params, link.function = link.function, var.function = var.function) # estimate the nuisance functions
    causal.estimator = DINAEstimator.step2(X = X[index.resample[foldIndex[[i]]],], Y = Y[index.resample[foldIndex[[i]]],], W = W[index.resample[foldIndex[[i]]]], family = family, nuisance.function = nuisance.function) # estimate the causal part
    estimator[[i]] = causal.estimator$estimator
    Sigma1[[i]] = causal.estimator$Sigma1
    Sigma2[[i]] = causal.estimator$Sigma2
  }
  estimator.ave = apply(matrix(unlist(estimator), byrow = T, nrow = k), 2, mean, na.rm = T) # average the estimated coefficients over folds
  names(estimator.ave) = names(estimator[[1]])
  result$estimator = estimator.ave

  Sigma1.ave = apply(array(unlist(Sigma1), dim = c(d+1,d+1,k)), c(1,2), mean)
  Sigma2.ave = apply(array(unlist(Sigma2), dim = c(d+1,d+1,k)), c(1,2), mean)
  result$Sigma = solve(Sigma1.ave) %*% Sigma2.ave %*% solve(Sigma1.ave)/n # compute the variance matrix of the coefficients
  class(result) = "DINA"
  return(result)
}


#' DINAEstimator.step1
#'
#' This function implements the first step of DINA: estimation of the nuisance functions.
#'
#' @param X matrix of covariates.
#' @param Y vector of responses.
#' @param W vector of treatment assignment.
#' @param family a character string describing the error distribution and link function to be used in the model. \code{Family} can be one of the following options: "gaussian", "binomial", "poisson", "cox". Default is "gaussian".
#' @param prop method for estimating the propensity score. \code{Prop} can be a character string naming the estimator: "logistic regression", "gradient boosting". \code{Prop} can also be a function and \code{prop(X)} should output the estimated propensity scores for each row of the covariate matrix \code{X}. Default is "logistic regression".
#' @param paramsW list of parameters for propensity score estimation.
#' @param baseline method for estimating the natural parameter functions in the treatment and control groups. \code{Baseline} can be a character string naming the estimator: "regression", "gradient boosting". \code{Baseline} can also be a list with two members: \code{eta1} denotes the estimator of the natural parameter function in the treatment group,   \code{eta0} denotes the estimator of the natural parameter function in the control group. The functions should take in a covariate matrix \code{X} and output the estimated quantities for each row of the covariate matrix \code{X}. Default is "regression".
#' @param params list of parameters for the natural parameter function estimation.
#' @param link.function link function from the conditional mean to the natural parameter.
#' @param var.function variance function from the conditional mean to the variance function.
#'
#' @return a list of data.
#'  \itemize{
#'  \item{a}{the estimated nuisance function a(x).}
#'  \item{nu}{the estimated nuisance function nu(x).}
#' }
#'
#' @export
DINAEstimator.step1 = function(X, Y, W, family, prop, paramsW = NULL, baseline = "gradient boosting", params, link.function = NA, var.function = NA){
  # preprocess
  n = dim(X)[1]; d = dim(X)[2]
  # estimate the fundamental quantities
  # estimate the propensity score
  if(is.character(prop)){
    if(prop == "logistic regression"){
      data.prop = data.frame(W, X); colnames(data.prop) = c("W", paste("X", seq(1,d), sep = ""))
      prop.model = glm(W~., data = data.prop, family = "binomial")
    } else if(prop == "gradient boosting"){
      prop.model = xgboost::xgboost(data = X, label = W, nrounds = paramsW$nrounds, params = paramsW$params, objective = "binary:logistic", eval_metric = "logloss", verbose = 0)
    }
  }else if(is.function(prop)){
    prop.model = prop # input the propensity score model
  }else{
    stop("please specify a method for estimating the propensity score")
  }

  # estimate the natural parameter functions
  if(family == "gaussian"){
    if(is.character(baseline)){
      if(baseline == "regression"){
        data.baseline = data.frame(Y, X); colnames(data.baseline) = c("Y", paste("X", seq(1,d), sep = ""))
        m.model = glm(Y~., data = data.baseline, family = family)
      }else if(baseline == "gradient boosting"){
        m.model = xgboost::xgboost(data = X, label = Y, nrounds = params$nrounds, params = params$params, objective = params$objective, eval_metric = params$eval_metric, verbose = 0)
      }
    }else if(is.function(baseline)){m.model = baseline
    }else{
      stop("please specify a method for estimating the marginal mean function")
    }

    # construct a(x) and nu(x)
    nuisance.function = function(X.test, prop.model.test = prop.model, m.model.test = m.model, prop.test = prop, baseline.test = baseline){
      data.test = data.frame(X.test); colnames(data.test) = c(paste("X", seq(1,d), sep = ""))
      if(is.character(prop.test)){
        if(prop.test == "logistic regression"){
          e.test = predict(prop.model.test, newdata = data.test, type = "response")
        } else if(prop.test == "gradient boosting"){
          e.test = predict(prop.model.test, newdata = X.test)}
      }else if(is.function(prop.test)){e.test = prop.model.test(X.test)}
      if(is.character(baseline.test)){
        if(baseline.test == "gradient boosting"){m.test = predict(m.model.test, newdata = X.test)}
        if(baseline.test == "regression"){m.test = predict(m.model.test, newdata = data.test, type = "response")}
      }else if(is.function(baseline.test)){m.test = m.model.test(X.test)}
      result = list(); result$a = e.test; result$nu = m.test
      return(result)
    }
  }else if(family %in% c("binomial", "poisson")){
    if(is.character(baseline)){
      if(baseline == "regression"){
        data.baseline = data.frame(Y, X); colnames(data.baseline) = c("Y", paste("X", seq(1,d), sep = ""))
        eta1.model = glm(Y~., data = data.baseline[W==1,], family = family)
        eta0.model = glm(Y~., data = data.baseline[W==0,], family = family)
      }else if(baseline == "gradient boosting"){
        eta1.model = xgboost::xgboost(data = X[W==1, , drop = F], label = Y[W==1], nrounds = params$nrounds, params = params$params, objective = params$objective, eval_metric = params$eval_metric, verbose = 0)
        eta0.model = xgboost::xgboost(data = X[W==0, , drop = F], label = Y[W==0], nrounds = params$nrounds, params = params$params, objective = params$objective, eval_metric = params$eval_metric, verbose = 0)
      }
    }else if(is.function(baseline)){eta1.model = baseline$eta1; eta0.model = baseline$eta0}

    nuisance.function = function(X.test, prop.model.test = prop.model, eta1.model.test = eta1.model, eta0.model.test = eta0.model, prop.test = prop, baseline.test = baseline){
      data.test = data.frame(X.test); colnames(data.test) = c(paste("X", seq(1,d), sep = ""))
      if(is.character(prop.test)){
        if(prop.test == "logistic regression"){e.test = predict(prop.model.test, newdata = data.test, type = "response")
        } else if(prop.test == "gradient boosting"){
          e.test = predict(prop.model.test, newdata = X.test)
        }
      }else if(is.function(prop.test)){e.test = prop.model.test(X.test)}
      if(is.character(baseline.test)){
        if(baseline.test == "regression"){
          mu1.test = predict(eta1.model.test, newdata = data.test, type = "response")
          mu0.test = predict(eta0.model.test, newdata = data.test, type = "response")
        }else if(baseline.test == "gradient boosting"){
          mu1.test = predict(eta1.model.test, newdata = X.test)
          mu0.test = predict(eta0.model.test, newdata = X.test)
        }
      }else if(is.function(baseline.test)){
        mu1.test = baseline.model.test$eta1(X.test)
        mu0.test = baseline.model.test$eta0(X.test)
      }
      a.test = (e.test * var.function(mu1.test))/(1e-6 + (1-e.test) * var.function(mu0.test)); a.test = a.test/(1 + a.test)
      nu.test = a.test * link.function(mu1.test) + (1-a.test) * link.function(mu0.test)

      result = list(); result$a = a.test; result$nu = nu.test
      return(result)
    }
  }else if (family == "cox"){
    if(is.character(baseline)){
      if(baseline == "regression"){
        data.baseline = data.frame(time = Y$time, status = Y$status, X = X); colnames(data.baseline) = c("time", "status", paste("X", seq(1,d), sep = ""))
        eta1.model = survival::coxph(survival::Surv(time, status) ~ ., data = data.baseline[W == 1,]) # no intercept: intercept is absorbed into the baseline hazard function
        eta0.model = survival::coxph(survival::Surv(time, status) ~ ., data = data.baseline[W == 0,])
        eta1.hat = predict(eta1.model, type = "lp")
        eta0.hat = predict(eta0.model, type = "lp")
      }else if(baseline == "gradient boosting"){
        eta1.model = xgboost::xgboost(data = X[W==1, , drop = F], label = Y$time[W==1] * (2 * Y$status[W==1] - 1), nrounds = params$nrounds, params = params$params, objective = "survival:cox", eval_metric = params$eval_metric, verbose = 0)
        eta0.model = xgboost::xgboost(data = X[W==0, , drop = F], label = Y$time[W==0] * (2 * Y$status[W==0] - 1), nrounds = params$nrounds, params = params$params, objective = "survival:cox", eval_metric = params$eval_metric, verbose = 0)
        eta1.hat = log(predict(eta1.model, newdata = X[W==1, , drop = F]))
        eta0.hat = log(predict(eta0.model, newdata = X[W==0, , drop = F]))
      }
    }else if(is.function(baseline)){
      eta1.model = baseline$eta1; eta0.model = baseline$eta0
      eta1.hat = eta1.model(X[W==1, , drop = F])
      eta0.hat = eta0.model(X[W==0, , drop = F])
    }
    # censoring probabilities
    not.censor.prob1.model = not.censoring.probability(X = X[W == 1, , drop = F], time = Y$time[W == 1], status = Y$status[W == 1], eta = eta1.hat, censorType = params$censorType, params = params)
    not.censor.prob0.model = not.censoring.probability(X = X[W == 0, , drop = F], time = Y$time[W == 0], status = Y$status[W == 0], eta = eta0.hat, censorType = params$censorType, params = params)

    nuisance.function = function(X.test, prop.model.test = prop.model, eta1.model.test = eta1.model, eta0.model.test = eta0.model, not.censor.prob1.model.test = not.censor.prob1.model, not.censor.prob0.model.test = not.censor.prob0.model, prop.test = prop, baseline.test = baseline){
      data.test = data.frame(X.test); colnames(data.test) = c(paste("X", seq(1,d), sep = ""))
      if(is.character(prop.test)){
        if(prop.test == "logistic regression"){e.test = predict(prop.model.test, newdata = data.test, type = "response")
        } else if(prop.test == "gradient boosting"){
          e.test = predict(prop.model.test, newdata = X.test)
        }
      }else if(is.function(prop.test)){e.test = prop.model.test(X.test)}
      if(is.character(baseline.test)){
        if(baseline.test == "regression"){
          eta1.test = predict(eta1.model.test, newdata = data.test, type = "lp")
          eta0.test = predict(eta0.model.test, newdata = data.test, type = "lp")
        }else if(baseline.test == "gradient boosting"){
          eta1.test = log(predict(eta1.model.test, newdata = X.test))
          eta0.test = log(predict(eta0.model.test, newdata = X.test))
        }
      }else if(is.function(baseline.test)){
        eta1.test = baseline.model.test$eta1(X.test)
        eta0.test = baseline.model.test$eta0(X.test)
      }
      not.censoring.probability1 = not.censor.prob1.model.test(eta.test = eta1.test, X.test = X.test)
      not.censoring.probability0 = not.censor.prob0.model.test(eta.test = eta0.test, X.test = X.test)

      a.test = (e.test * not.censoring.probability1)/(1e-6 + (1-e.test) * not.censoring.probability0); a.test = a.test/(1 + a.test)
      nu.test = a.test * eta1.test + (1-a.test) * eta0.test

      result = list(); result$a = a.test; result$nu = nu.test
      return(result)
    }
  }
  return (nuisance.function)
}


#' DINAEstimator.step2
#'
#' This function implements the second step of DINA: estimation of the causal parameter given nuisance function estimators.
#'
#' @param X matrix of covariates.
#' @param Y vector of responses.
#' @param W vector of treatment assignment.
#' @param nuisance.function nuisance function estimators. If we call nuisance.function.estimator = nuisance.function(X), then nuisance.function.estimator is a list with two members: \code{a} denotes the estimated nuisance function a(x), \code{nu} denotes the estimated nuisance function nu(x).
#'
#' @export
DINAEstimator.step2 = function(X, Y, W, family, nuisance.function){
  nuisance.function.estimator = nuisance.function(X)

  if(family != "cox"){
    data.HTE = data.frame(Y, nuisance.function.estimator$nu, W - nuisance.function.estimator$a, I(X)); colnames(data.HTE) = c("Y", "nu", "resW", "X")
    glm.model = tryCatch(glm(Y ~ resW + X:resW + offset(nu) - 1, data = data.HTE, family = family), error=function(e) e, warning=function(w) w)
    if(is(glm.model,"warning")){
      stop("Warning from glm.fit.")
    } else {estimator = coef(glm.model)}
    result = list(); result$estimator = estimator
    result$Sigma1 = t(cbind(1, X)) %*% diag((data.HTE$resW * residuals(glm.model, type = "response") / residuals(glm.model, type = "pearson"))^2) %*% cbind(1, X)/length(Y)
    result$Sigma2 = t(cbind(1, X)) %*% diag((data.HTE$resW * residuals(glm.model, type = "response"))^2) %*% cbind(1, X)/length(Y)
  }else{
    # Cox model
    data.HTE = data.frame(Y$time, Y$status, nuisance.function.estimator$nu, W - nuisance.function.estimator$a, I(X)); colnames(data.HTE) = c("time", "status", "nu", "resW", "X")
    glm.model = tryCatch(survival::coxph(survival::Surv(time, status) ~ resW + X:resW + offset(nu) - 1, data = data.HTE), error=function(e) e, warning=function(w) w)
    if(is(glm.model,"warning")){estimator = rep(NA, d+1)
    } else {estimator = coef(glm.model)}
    result = list(); result$estimator = estimator
    result$Sigma1 = result$Sigma2 = solve(glm.model$var) / length(Y$time)
  }
  return(result)
}


#' DINA.bagging.inference
#'
#' This function computes the variance matrix of bagged DINA estimators.
#'
#' @param bagging.estimator matrix of DINA estimators. Each row stands for a bootstrap resample, each column stands for a feature.
#' @param bootstrap.membership matrix of how many times a unit appears in a bootstrap resample. Each row stands for a bootstrap resample, each column stands for a unit.
#' @param method method for estimating the variance matrix. \code{Method} can be a character string naming the estimator: "infinitesimal jackknife", "jackknife". Default is "infinitesimal jackknife".
#' @param bias.correction logical value indicating whether to execute bias correction. If \code{TRUE}, an estimator of the bias will be subtracted from the original estimator. Default is \code{TRUE}.
#'
#' @return variance matrix of bagged DINA estimators.
#'
#' @export
DINA.bagging.inference = function(bagging.estimator, bootstrap.membership, method = "infinitesimal jackknife", bias.correction = T){
  # preprocess
  indexNA = which(apply(is.na(bagging.estimator), 1, sum) > 0)
  if(length(indexNA) > 0){
    bagging.estimator = bagging.estimator[-indexNA,]
    bootstrap.membership = bootstrap.membership[-indexNA, ]
  }
  n = dim(bootstrap.membership)[2]
  B = dim(bootstrap.membership)[1]

  if(method == "infinitesimal jackknife"){
    bagging.estimator.center = scale(bagging.estimator, center = T, scale = F)
    Cov = t(bootstrap.membership - 1) %*% bagging.estimator.center/B # compute the covariance of the bootstrap.membership and the bagging.estimator, Cov is a n by (d+1) matrix,
    Sigma = t(Cov) %*% Cov
    if(bias.correction){
      Sigma = Sigma - n / B^2 * t(bagging.estimator.center) %*% bagging.estimator.center
    }
  }else if(method == "jackknife"){
    bagging.estimator.center = scale(bagging.estimator, center = T, scale = F)
    bootstrap.membership = 1 - pmin(bootstrap.membership, 1)
    Delta = diag(1/pmax(apply(bootstrap.membership, 2, sum), 1)) %*% t(bootstrap.membership) %*% bagging.estimator.center # compute the difference of the average of the estimators without the i-th sample and the bagging average, Delta is a n by (d+1) matrix
    Sigma = t(Delta) %*% Delta * (n-1)/n
    if(bias.correction){
      Sigma = Sigma - (exp(1) - 1) * n / B^2 * t(bagging.estimator.center) %*% bagging.estimator.center
    }
  }
  return(Sigma)
}


