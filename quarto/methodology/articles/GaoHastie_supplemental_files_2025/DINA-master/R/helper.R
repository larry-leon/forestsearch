# helper functions for DINA

#' dataGnr
#'
#' This function generates the simulated data.
#'
#' @param n sample size.
#' @param d number of covariates.
#' @param alpha coefficient vector for the control natural parameter function. Without further specification, the control natural parameter function is the sum of two parts: a linear function of the covariates \code{X} (we include the intercept as the last column of \code{X}) with the coefficient vector \code{alpha}, and a non-linear part determined by the type of misspecification \code{typeMisMu0} and the magnitude \code{delta}.
#' @param beta coefficient vector for the heterogeneous treatment effect. Without further specification, the heterogeneous treatment effect is a linear function of the covariates \code{X} (we include the intercept as the last column of \code{X}) with the coefficient vector \code{beta}.
#' @param theta coefficient vector for the propensity score. Without further specification, the propensity score follows a logistic regression model of the covariates \code{X} (we include the intercept as the last column of \code{X}) with the coefficient vector \code{theta}.
#' @param e constant propensity score value. If \code{e} is not zero, then the propensity score is set at \code{e} for all units. Default is zero.
#' @param typeX distribution of the covariates. Default is the uniform distribution U(0,1).
#' @param typeProp model of the propensity score. \code{TypeProp} is used if \code{e = 0}, i.e., the propensity score is not constant. The propensity score can depend on the covariates in the following ways: "logistic", "linear". Default is "logistic".
#' @param family family of responses. The response can be one of the following families: "gaussian", "binomial", "poisson", "exponential", "gamma", "cox".
#' @param typeMisMu0 type of the non-linear part of the control natural parameter function. The non-linear part can be one of the following types: "quadratic", "interaction", "three-interaction", "sin", "cos", "abs", "unobserved", "unobserved linear", "unobserved interaction".
#' @param delta the magnitude of the non-linear part of the control natural parameter function.
#' @param params list of parameters. For censoring, \code{params$censorType} specifies the distribution of the censoring times. The censoring times can be one of the following types: NULL (no censoring), "exponential", "cox", "uniform", "fixed".
#'
#' @return a list of data.
#'  \itemize{
#'  \item{X}{covariate matrix.}
#'  \item{Y}{vector of response.}
#'  \item{W}{vector of treatment assignment.}
#'  \item{prop}{vector of propensity score.}
#'  \item{mu0}{vector of control natural parameter function.}
#'  \item{tau}{heterogeneous treatment effect.}
#' }
#'
#' @examples
#'data.gaussian = dataGnr(n = 100, d = 1, alpha = rep(1,2), beta = rep(0.5,2), e = 0.5, delta = 1, family = "gaussian", typeMisMu0 = "quadratic")
#'data.binomial = dataGnr(n = 100, d = 5, alpha = rep(1,6), beta = rep(0.5,6), e = 0.5, delta = 1, family = "binomial", typeMisMu0 = "quadratic")
#'params = list(); params$censorType = "uniform"; params$T0 = 5
#'data.cox = dataGnr(n = 100, d = 5, alpha = rep(1,6), beta = rep(0.5,6), e = 0.5, delta = 1, family = "cox", typeMisMu0 = "quadratic", params = params)
#'
#' @export
dataGnr = function(n, d, alpha, beta, e = 0, theta = NULL, delta = 0, typeX = "uniform", typeProp = "logistic", family, typeMisMu0 = "interaction", params = list()){
  # X
  if(typeX == "uniform") {
    X = matrix(runif(d * n, -1, 1), ncol = d)
  } else if(typeX == "gaussian"){
    X = matrix(rnorm(d * n, 0, 1), ncol = d)
  } else if(typeX == "mixUniform") {
    X = matrix(runif(d * n, -1, 1), ncol = d)
    if(is.null(params$mixtureGap)){params$mixtureGap = (max(X[,1]) - min(X[,1]))/2}
    if(is.null(params$mixtureProp)){params$mixtureProp = 0.5}
    X[,1] = X[,1] + c(rep(params$mixtureGap/2, floor(n * params$mixtureProp)), rep(-params$mixtureGap/2, n-floor(n * params$mixtureProp)))
    X = X[sample(n),,drop = F] # randomly permute the rows of X
  }
  # prop and W
  if(e != 0) {prop = rep(e, n)
  } else if (typeProp == "logistic") {
    prop = cbind(X, 1) %*% theta
    prop = c(exp(prop)/(1 + exp(prop)))
  } else if(typeProp == "linear"){prop = c(cbind(X, 1) %*% theta)}
  W = rbinom(n, 1, prop)
  # tau
  tau = c(cbind(X, 1) %*% beta)
  # mu0
  mu0 = c(cbind(X, 1) %*% alpha)
  if (typeMisMu0 == "quadratic") {mu0 = mu0 + delta * (X[,1])^2
  } else if (typeMisMu0 == "interaction") {
    mu0 = mu0 + delta * X[,1] * X[,2]
  } else if (typeMisMu0 == "three-interaction") {
    mu0 = mu0 + delta * X[,1] * X[,2] * X[,3]
  } else if (typeMisMu0 == "sin") {mu0 = mu0 + delta * sin(X[,1])
  } else if (typeMisMu0 == "cos") {mu0 = mu0 + delta * cos(X[,1])
  } else if (typeMisMu0 == "abs") {mu0 = mu0 + delta * abs(X[,1])
  } else if (typeMisMu0 == "unobserved") {
    mu0 = mu0 + delta * runif(n,-1,1)
  } else if (typeMisMu0 == "unobserved linear") {
    mu0 = mu0 + delta * X[,1] * runif(n,-1,1)
  } else if (typeMisMu0 == "unobserved interaction") {
    mu0 = mu0 + delta * X[,1] * X[,2] * runif(n,-1,1)
  }
  # Y
  if(family == "gaussian"){
    if(is.null(params$sigma)){params$sigma = 1}
    Y = W * tau + mu0 + rnorm(n, 0, params$sigma)
  } else if (family == "binomial"){
    Y = rbinom(n, 1, prob = exp(W*tau+mu0)/(1+exp(W*tau+mu0)))
  } else if (family == "poisson"){
    Y = rpois(n = n, lambda = exp(W * tau + mu0))
  } else if (family == "exponential"){
    Y = rexp(n = n, rate = exp(W * tau + mu0))
  } else if (family == "gamma"){
    if(is.null(params$shape)){params$shape = 1}
    Y = rgamma(n, shape = params$shape, rate = exp(W*tau+mu0))
  } else if (family == "cox"){
    # Weibull distribution
    if(is.null(params$lambda)){params$lambda = 1}
    if(is.null(params$rho)){params$rho = 1}
    time = (-log(runif(n = n))/(params$lambda * exp(mu0 + W*tau)))^(1/params$rho)
    # censoring times
    if(!is.null(params$censorType)){
      if(params$censorType == "exponential"){
        if(is.null(params$rateC)){params$rateC = 1}
        C = rexp(n = n, rate = params$rateC)
      } else if (params$censorType == "cox"){
        if(is.null(params$censorLambda)){params$censorLambda = 1}
        if(is.null(params$censorRho)){params$censorRho = 1}
        C = (-log(runif(n = n))/params$censorLambda)^(1/params$censorRho)
      } else if (params$censorType == "uniform"){
        if(is.null(params$T0)){params$T0 = max(Y)+1}
        C = runif(n = n, 0, params$T0)
      } else if (params$censorType == "fixed"){
        if(is.null(params$T0)){params$T0 = max(Y)+1}
        C = rep(params$T0, n)}
      status = as.numeric(time <= C)
      time = pmin(time, C)
    } else {
      status = rep(1, length(time)) # no censoring
    }
    # response
    Y = data.frame(time = time, status = status)
  }

  record = list(); record$X = X; record$Y = Y; record$prop = prop; record$W = W; record$mu0 = mu0; record$tau = tau
  return (record)
}


#' split.uniform
#'
#' This function splits samples into folds with equal fold size.
#'
#' @param n sample size.
#' @param k number of folds. Default is 2.
#' @param replace logical value indicating the whether the sampling should be with replacement. If \code{TRUE}, sampling will be with replacement. Default is \code{FALSE}.
#'
#' @return a list indicating the resulting folds.
#'  \itemize{
#'  \item{index.resample}{vector of subject indices in the permuted sample. We permute the original dataset and each subject will appear exactly once in the permuted sample. The permutation is done to break inherent ordering in the original dataset (if there is any) and makes sure the resulted folds are exchangeable.}
#'  \item{foldIndex}{list of fold memberships. The j-th member is a vector containing all subjects of the bootstrap resample assigned to the j-th fold.}
#' }
#'
#' @examples
#' set.seed(318)
#' test = split.uniform(n, k)
#' test$index.resample[test$foldIndex[[1]]]
#'
#' @export
split.uniform = function(n, k = 2, replace = F){
  if(n < 2 * k){stop("please use a smaller k.")}
  result = list()
  result$index.resample = sample(n, n, replace = replace) # randomly permute the indices
  result$foldIndex = list()
  for(i in 1:k){
    result$foldIndex[[i]] = which((1 + seq(0, n - 1) %% k) == i)
  }
  return(result)
}


#' split.weight
#'
#' This function splits a sample with duplicates into non-overlap folds with equal size.
#'
#' @param n sample size.
#' @param k number of folds. Default is 2.
#'
#' @return a list indicating the resulting folds.
#'  \itemize{
#'  \item{index.resample}{vector of subject indices in the bootstrap resample. For example, assume there are 3 subjects in the dataset and consider a bootstrap resample consisting of two copies of the first subject and one copy of the second subject. In this way, \code{index.resample} is the vector [1,1,2].}
#'  \item{foldIndex}{list of fold memberships. The j-th member is a vector containing all subjects of the bootstrap resample assigned to the j-th fold.}
#'  \item{weight}{vector of weights for the original dataset. The j-th coordinate stands for the times that the j-th subject of the original dataset appears in the bootstrap resample.}
#' }
#'
#' @examples:
#' set.seed(318)
#' test = split.weight(n, k)
#' test$index.resample[test$foldIndex[[1]]]
#' test$index.resample[test$foldIndex[[2]]]
#' test$weight
split.weight = function(n, k = 2){
  if(n < 2 * k){stop("please use a smaller k.")}
  result = list()
  while(TRUE){
    foldIndex = list()
    weight = rmultinom(n = 1, size = n, prob = rep(1/n,n))[,1]

    index = sample(n, n, replace = F)
    index.resample = rep(index, weight[index]) # randomly permute the indices

    # split subjects appearing at least twice in the bootstrap resample
    index2 = which(weight[index]>1)
    if(length(index2) > 0){
      foldSize2 = rep(floor(length(index2)/k), k) # what happen if it is zero
      if(length(index2) %% k != 0){
        foldSize2[seq(1,length(index2) %% k)] = foldSize2[seq(1,length(index2) %% k)] + 1
      }
      cumFoldSize2 = c(0,cumsum(foldSize2))
      for(i in 1:(k-1)){
        foldIndex[[i]] = which(index.resample %in% index[index2[seq(1+cumFoldSize2[i], cumFoldSize2[i+1])]])
      }
      foldIndex[[k]] = which(index.resample %in% index[index2[seq(1+cumFoldSize2[k], length(index2))]])
    }else{
      for(i in 1:k){foldIndex[[i]] = numeric(0)}
    }

    # split subjects appearing exactly once in the bootstrap resample
    index1 = which(weight[index]==1)
    foldSize1 = floor(n/k) - as.numeric(lapply(foldIndex, length)) # while true
    if(min(foldSize1) < 1){next} # ensure equal weight sums across folds
    cumFoldSize1 = c(0, cumsum(foldSize1))
    for(i in 1:(k-1)){
      foldIndex[[i]] = c(foldIndex[[i]], which(index.resample %in% index[index1[seq(1+cumFoldSize1[i], cumFoldSize1[i+1])]]))
    }
    foldIndex[[k]] = c(foldIndex[[k]], which(index.resample %in% index[index1[seq(1+cumFoldSize1[k], length(index1))]]))
    break
  }

  result$index.resample = index.resample; result$foldIndex = foldIndex; result$weight = weight
  return(result)
}


#' not.censoring.probability
#'
#' This function estimates the probability of not getting censored.
#'
#' @param censorType type of censoring. Types of censoring can be one of the following options: NULL (no censoring), "exponential", "cox", "uniform", "fixed", "unknown". For censoring types except NULL and "unknown", we compute the probability of getting censored in closed form or numerically. If the type of censoring is "unknown" or "exponential", a classifier is trained with the response whether a unit is censored and all available covariates as predictors.
#'
#' @return a function takes in covariate values and outputs the probabilities of not getting censored.
#' @export
not.censoring.probability = function(X = NULL, time = NULL, status = NULL, eta = NULL, censorType = NULL, params = list()){
  if(is.null(censorType)){
    not.censoring.probability.test = function(eta.test = NULL, X.test){
      prob = rep(1, dim(X.test)[1])
      return(prob)
    }
  }else if(censorType %in% c("unknown")){
    prob.model = xgboost::xgboost(data = X, label = status, nrounds = params$nrounds, params = params$params, objective = "reg:logistic", eval_metric = "logloss", verbose = 0)
    not.censoring.probability.test = function(eta.test = NULL, X.test, prob.model.test = prob.model){
      prob = predict(prob.model.test, newdata = X.test)
      return(prob)
    }
  }else if(censorType == "fixed"){
    cumulative.hazard = base.cumulative.hazard(eta = eta, time = time, status = status, time.test = params$T0)
    not.censoring.probability.test = function(eta.test, X.test = NULL, cumulative.hazard.test = cumulative.hazard){
      prob = 1 - exp(-exp(eta.test) * cumulative.hazard.test)
      return(prob)
    }
  }else if(censorType == "uniform"){
    time.test = seq(0, params$T0, length.out = 500)[-1]
    cumulative.hazard = base.cumulative.hazard(eta = eta, time = time, status = status, time.test = time.test)
    censor.density = rep(1/params$T0, length(time.test))
    not.censoring.probability.test = function(eta.test, X.test = NULL,  cumulative.hazard.test = cumulative.hazard, censor.density.test = censor.density, time.test2 = time.test){
      prob = not.censoring.probability.helper(eta = eta.test, cumulative.hazard = cumulative.hazard.test, censor.density = censor.density.test, time = time.test2)
      return(prob)
    }
  }else if(censorType == "exponential"){
    time.test = seq(0, max(time) * 2, length.out = 500)[-1]
    cumulative.hazard = base.cumulative.hazard(eta = eta, time = time, status = status, time.test = time.test) # eta, time, status are computed for the training data; time.test are the time spots to be evaluated at
    censor.density = dexp(time.test, rate = params$rateC)
    not.censoring.probability.test = function(eta.test, X.test = NULL, cumulative.hazard.test = cumulative.hazard, censor.density.test = censor.density, time.test2 = time.test){
      prob = not.censoring.probability.helper(eta = eta.test, cumulative.hazard = cumulative.hazard.test, censor.density = censor.density.test, time = time.test2)
      return(prob)
    }
  }
  return(not.censoring.probability.test)
}


#' not.censoring.probability.helper
#'
#' This function is a helper for the function \code{not.censoring.probability}.
#'
#' @param eta vector of linear predictors in a fitted Cox model.
#' @param cumulative.hazard vector of baseline cumulative hazard at \code{time}.
#' @param censor.density vector of the density of censoring times evaluated at \code{time}.
#' @param time vector of survival times.
#'
#' @return vector of the probabilities of not getting censored.
#'
#' @export
not.censoring.probability.helper = function(eta, cumulative.hazard, censor.density, time){
  time.order = order(time, decreasing = F)
  time = sort(time, decreasing = F)
  time.cdf = 1 - exp(-matrix(exp(eta), ncol = 1) %*% matrix(cumulative.hazard[time.order], nrow = 1)) # each column stands for a time, and each row stands for a unit
  prob = c(time.cdf %*% (censor.density[time.order] * diff(c(0, time), 1)))
  return(prob)
}


#' base.cumulative.hazard
#'
#' This function computes the baseline cumulative hazard for a fitted Cox model.
#'
#' @param eta vector of linear predictors in a fitted Cox model.
#' @param time vector of survival times.
#' @param status vector of censored status. \code{Status = 1} stands for not censored, and \code{status = 0} stands for being censored.
#' @param time.test vector of time to evaluate the baseline cumulative hazard at.
#'
#' @return vector of baseline cumulative hazard.
#'
#' @export
base.cumulative.hazard = function(eta, time, status, time.test = NULL){
  # preprocess
  if(is.null(time.test)){time.test = time} # if not time.test is input, then evaluate the baseline cumulative hazard at time
  time.order = order(time, decreasing = T)
  lambda = rev(1/cumsum(exp(eta[time.order]))) # compute lambda_i = 1/(sum_{t_k >= t_i} exp(eta_k)); lambda is ordered increasingly in time
  Lambda = cumsum(lambda * rev(status[time.order])) # compute Lambda_i = (sum_{t_k <= t_i, k not censored} exp(eta_k)); Lambda is ordered increasingly in time should be there
  Lambda.test = approx(x = rev(time[time.order]), y = Lambda, xout = time.test, method = "linear", rule = 2, ties = "mean")$y # predicted cumulative baseline hazards at time.test based on the linear interpolation of Lambda. Outside the interval of time, the value at the closest time extreme is used
  return(Lambda.test)
}


#' cox.inference
#'
#' This function computes the first and second order derivatives of the partial likelihood of a Cox model with regard to the regression coefficients.
#'
#' @param eta vector of linear predictors in a fitted Cox model.
#' @param time vector of survival times.
#' @param status vector of censored status. \code{Status = 1} stands for not censored, and \code{status = 0} stands for being censored.
#' @param X covariate matrix. Each row stands for an observation, each column stands for a covariate.
#'
#' @return  a list of data.
#'  \itemize{
#'  \item{U}{the first order derivative of the partial likelihood.}
#'  \item{I}{the second order derivative of the partial likelihood.}
#' }
#'
#' @export
cox.inference = function(eta, time, status, X){
  time.order = order(time, decreasing = F)
  eta.sort = eta[time.order] # eta ordered increasingly w.r.t. the survival response
  X.order = X[time.order,,drop = F]
  weight.matrix = t(replicate(length(time), exp(eta.sort))) # each row stands for a unit, each column stands for a candidate of the risk set
  weight.matrix[lower.tri(weight.matrix, diag = F)] = 0 # set weights not in the risk set as zero
  weight.matrix = diag(1/apply(weight.matrix, 1, sum)) %*% weight.matrix # normalize the row weight sum to be one
  weight.matrix = diag(status[time.order]) %*% weight.matrix # set the weights of censored units to be zero
  X.order.average = weight.matrix %*% X.order # compute the average of covariates in the risk set with weights proportional to the hazards

  U = apply(diag(status[time.order]) %*% X.order - X.order.average, 2, sum) # compute the first order derivative of the partial likelihood

  I = t(X.order) %*% diag(apply(weight.matrix, 2, sum)) %*% X.order - t(X.order.average) %*% X.order.average # compute the second order derivative of the partial likelihood

  result = list(); result$U = U; result$I = I
  return(result)
}


#' psd.projection
#'
#' This function computes the matrix projection onto the positive-semi-definite cone with respect to the Frobenius or operator norm.
#'
#' @param Sigma symmeteric matrix.
#'
#' @return the projected matrix in the positive-semi-definite cone.
#'
#' @export
psd.projection = function(Sigma){
  Sigma.eigen = eigen(Sigma)
  Sigma.psd = Sigma.eigen$vectors %*% diag(pmax(0, Sigma.eigen$values)) %*% t(Sigma.eigen$vectors)
  return(Sigma.psd)
}

