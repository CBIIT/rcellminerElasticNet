#' Summary of Elastic Net Results 
#'
#' @param object an object of type enResults
#' @param ... additional parameters 
#' @return Nothing
#' 
#' @concept rcellminerElasticNet
#' @export
summary.enResults <- function(object, ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("\nAlpha:\n", object$alpha, "\n")
  cat("\nLambda:\n", object$lambda, "\n")
  
  cat("\nFolds:\n", max(object$foldIds), "\n")
  cat("\nTraining Runs:\n", dim(object$featureWtMat)[2], "\n")
}
