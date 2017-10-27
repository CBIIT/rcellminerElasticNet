#' Add features that are correlated with a set of elastic net predictors.
#' 
#' @param enResults an elastic net results object.
#' @param dataMatList a list of matrices with feature data organized along the rows,
#' and feature names accessible via rownames(dataMatList).  These matrices must contain
#' data for all features used in the original elastic net call.
#' @param corThreshold a correlation threshold.
#' @return an elastic net results object containing all data in the enResults parameter,
#' plus two added elements:
#' \itemize{
#'   \item{predictorCorFeatures} {is a list of character vectors indexed by predictor name.  
#' Each character vector associated with a predictor contains the names of all features correlated 
#' above corThreshold with that predictor. (Correlations are with respect to features used as 
#' input to the original elastic net call, as recorded in rownames(enResults$featureWtMat).)}
#'   \item{predictorsAndCorFeatures} {is a character vector containing all predictors together 
#' with features correlated with them above corThreshold.}
#' }
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
#' @importFrom rcellminer selectCorrelatedRows getFeatureDataFromMatList
addPredictorCorrelatedFeatures <- function(enResults, dataMatList, corThreshold=0.8){
  features <- setdiff(rownames(enResults$featureWtMat), "(Intercept)")
  featureData <- getFeatureDataFromMatList(features, dataMatList, excludeMissingFeatures=TRUE)
  
  if (!identical(features, rownames(featureData))){
    missingFeatures <- setdiff(features, rownames(featureData))
    stop(paste("There is no data for the following features in dataMatList: ", 
               paste(missingFeatures, collapse=" "), ".", sep=""))
  }
  
  predNames <- names(enResults$predictorWts)
  names(predNames) <- predNames
  
  enResults$predictorCorFeatures <- lapply(predNames, function(pred){
    predData <- getFeatureDataFromMatList(pred, dataMatList)
    corFeatures <- rownames(selectCorrelatedRows(Y=predData, featureData, 
                                                 corThreshold=corThreshold, useAbsCor=FALSE))
    return(setdiff(corFeatures, pred))
  })
  
  enResults$predictorsAndCorFeatures <- unique(c(names(enResults$predictorWts),
                                                 unname(c(enResults$predictorCorFeatures, recursive = TRUE))))
  
  return(enResults)
}
