#' coef.DINA
#'
#' This function extracts the coefficients of the heterogeneous treatment effect from a DINA object.
#' @method coef DINA
#' @importFrom stats coef
#'
#' @param object object of class "DINA".
#' @param ... further arguments passed to or from other methods.
#'
#' @return a named vector of coefficients extracted from the DINA object.
#'
#' @export
coef.DINA = function(object, ...){
  # preprocess
  if(class(object) != "DINA"){stop("Please input an object of class DINA.")}

  # predict
  result = object$estimator
  names(result) = c("(Intercept)", paste("X", seq(1,length(result)-1), sep = ""))

  return(result)
}


