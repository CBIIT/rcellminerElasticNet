#' Apply a fit linear model to predict the response based on new data.
#' 
#' @param model A fit linear model (currently restricted to models of class enResults).
#'  If this parameter is provided, coeffVec and yIntercept parameters are not needed,
#'  and will be ignored if provided.
#' @param useModelYIntercept A boolean value indicating whether to use the intercept
#'  term provided by the model (default = TRUE)
#' @param coeffVec A named vector of coefficient weights (which need not be specified
#'  if model parameter is provided instead).
#' @param yIntercept An intercept term (default = 0)
#' @param newData A matrix with predictor variable data over observations (specified
#'  along rows or columns).
#' 
#' @return A vector of predicted response values.
#'  
#' @concept rcellminerElasticNet
#' @export
predictWithLinRegModel <- function(model=NULL, useModelYIntercept, coeffVec=NULL, 
                                   yIntercept=0, newData){
  if (is.null(model) && is.null(coeffVec)){
    stop("Either model or coeffVec parameter must be provided to predict response.")
  }
  if (!is.null(model)){
    if ("enResults" %in% class(model)){
      coeffVec <- model$predictorWts
      if (useModelYIntercept){
        yIntercept <- model$yIntercept
      }
    } else{
      stop("model parameter must be of class enResults.")
    }
  }
  
  if (all(names(coeffVec) %in% rownames(newData))){
    newData <- t(newData)
  }
  if (!all(names(coeffVec) %in% colnames(newData))){
    stop("All variable names must appear in rownames or column names of newData.")
  }
  
  N <- nrow(newData)
  X <- cbind(rep(1, N), newData)
  betaVec <- setNames(rep(0, ncol(X)), colnames(X))
  betaVec[names(coeffVec)] <- coeffVec
  betaVec[1] <- yIntercept
  y <- setNames(as.numeric(X %*% betaVec), rownames(X))
  
  return(y)
}