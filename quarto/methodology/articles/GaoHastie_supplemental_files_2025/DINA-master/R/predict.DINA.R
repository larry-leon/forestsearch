#' predict.DINA
#'
#' This function computes the predicted heterogeneous treatment effect based on a DINA object.
#' @method predict DINA
#' @importFrom stats predict
#'
#' @param object object of class "DINA".
#' @param newdata data frame or matrix of covariates with which to predict.
#' @param ... further arguments passed to or from other methods.
#'
#' @return a vector of predictions.
#'
#' @export
predict.DINA = function(object, newdata, ...){
  # preprocess
  if(class(object) != "DINA"){stop("Please input an object of class DINA.")}
  if(!class(newdata)[1] %in% c("array", "matrix", "data.frame")){stop("Please input a matrix or data frame as the newdata.")}

  # predict
  result = c(cbind(1, newdata) %*%  object$estimator)

  return(result)
}
