#' summary.DINA
#'
#' \code{summary} method for class "DINA".
#'
#' @param object a DINA object
#'
#' @return  a list of summary statistics of the given DINA \code{object}.
#'  \itemize{
#'  \item{coefficients}{a matrix with columns for the estimated coefficient, its standard error, and the corresponding (two-sided) p-value. Each row stands for an predictor (intercept included).}
#'  \item{params}{hyper-parameters used to obtain the given DINA \code{object}.}
#' }
#'
#' @export
summary.DINA = function(object){
  result = list()
  result$coefficients =  cbind(Estimate = object$estimator,
                               `Std. Error` = sqrt(diag(object$Sigma)),
                               `P-value` = 2 * pnorm(-abs(object$estimator)/sqrt(diag(object$Sigma))))
  result$params = params

  return(result)
}


