#' Estimate the prediction error of a linear model using repeated k-fold
#' cross validation.
#' 
#' @param X A data matrix with observation data along the rows and predictor data
#'  along the columns.
#' @param y A response vector.
#' @param nFolds The number of folds used in cross-validation.
#' @param nRepeats The number of cross-validation repeats (with different fold
#'  assignments).
#' @param cvFoldIds An nObs x nReplicates matrix with integer fold identifiers 
#'  along the columns.
#' 
#' @return A list with the following elements
#' \itemize{
#' \item{cvPred} {A vector of cross-validation predicted values (averaged over cross validation
#'  repeats).}
#' \item{cvPredR} {The Pearson's correlation between the above vector of predicted response 
#'  values and the actual response values.}
#' \item{cvPredRsqared} {The square of the Pearson's correlation between the above vector of
#'  predicted response values and the actual response values.}
#' \item{cvMeanSqErr} {The mean squared error (averaged over the results from nFolds x nRepeats
#'  sets of cross-validation predictions).}
#' \item{cvSdMeanSqErr} {The standard deviation of the above (nFolds x nRepeats) cross validation
#'  mean squared error values.}
#' \item{cvPredMat} {Matrix with predicted cross-validation predictions}
#' \item{cvMeanSqErrMat} {Matrix with the cross-validation mean squared errors}
#' }
#' 
#' @concept rcellminerElasticNet
#' @export
getLmCvFit <- function(X, y, nFolds=10, nRepeats=10, cvFoldIds=NULL){
  if (is.null(cvFoldIds)){
    cvFoldIds <- getCvFoldIds(nObs=nrow(X), nFolds=nFolds, nRepeats=nRepeats)
  } else{
    nFolds <- max(as.numeric(cvFoldIds))
    nRepeats <- ncol(cvFoldIds)
  }
  
  lmData <- data.frame(cbind(y, X))
  colnames(lmData) <- c("y", paste0("x", 1:ncol(X)))
  
  cvPredMat <- matrix(NA, nrow=nrow(X), ncol=nRepeats)
  rownames(cvPredMat) <- rownames(X)
  cvMeanSqErrMat <- matrix(NA, nrow=nFolds, ncol=nRepeats)
  
  for (j in seq_len(nRepeats)){
    for (i in seq_len(nFolds)){
      iTest <- which(cvFoldIds[, j] == i)
      iTrain <- setdiff((1:nrow(lmData)), iTest)
      lmFit <- lm(formula = y ~ ., data = lmData[iTrain, ])
      yPred <- predict(lmFit, newdata = lmData[iTest, -1, drop=FALSE])
      
      cvPredMat[iTest, j] <- yPred
      cvMeanSqErrMat[i, j] <- mean((lmData[iTest, "y"] - yPred)^2)
    }
  }
  
  lmCvFitResults <- list()
  lmCvFitResults$cvPred <- rowMeans(cvPredMat)
  lmCvFitResults$cvPredR <- cor.test(lmCvFitResults$cvPred, lmData$y)$estimate
  lmCvFitResults$cvPredRsqared <- (lmCvFitResults$cvPredR)^2
  lmCvFitResults$cvMeanSqErr <- mean(as.numeric(cvMeanSqErrMat))
  lmCvFitResults$cvSdMeanSqErr <- sd(as.numeric(cvMeanSqErrMat))
  lmCvFitResults$cvPredMat <- cvPredMat
  lmCvFitResults$cvMeanSqErrMat <- cvMeanSqErrMat
  
  return(lmCvFitResults)
}