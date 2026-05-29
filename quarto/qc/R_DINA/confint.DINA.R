#' confint.DINA
#'
#' This function computes confidence intervals for one or more parameters in a DINA object. The method assumes normality.
#'
#' @method confint DINA
#' @importFrom stats confint
#'
#' @param object object of class "DINA".
#' @param level the confidence level required. Default is 0.95.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a matrix with columns giving lower and upper confidence limits for each parameter. These will be labeled as (1-\code{level})/2 and 1 - (1-\code{level})/2 in % (by default 2.5% and 97.5%).
#' @export
confint.DINA = function(object, level = 0.95, ...){
  # preprocess
  if(class(object) != "DINA"){stop("Please input an object of class DINA.")}

  # compute confidence intervals
  c = qnorm(1 - (1-level)/2)
  result = cbind(object$estimator - c * sqrt(diag(object$Sigma)), object$estimator + c * sqrt(diag(object$Sigma)))

  # postprocess
  rownames(result) = names(coef(object))
  colnames(result) = paste(100 * c((1-level)/2, 1 - (1-level)/2), "%", sep = "")
  return(result)
}
