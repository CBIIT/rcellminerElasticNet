#' Cross validate the entire elastic net procedure to obtain a set of predicted
#' response values for strictly held out data.
#'
#' @param nFolds The number of folds used in cross-validation.
#' @param nRepeats The number of cross-validation repeats (with different fold
#'  assignments).
#' @param cvFoldIds An nObs x nReplicates matrix with integer fold identifiers 
#'  along the columns.
#' @param minFeatureResponseAbsCor Minimum absolute correlation between feature
#'  and response (used to filter rows of featureMat).
#' @param featureMat p x n matrix with input feature vectors along rows.
#' @param responseVec n-dimensional response vector to be predicted using a sparse
#'   linear combination of input feature vectors specified in featureMat.
#' @param standardize Logical flag for feature variable standardization (across
#'   observations), passed in to glmnet functions (glmnet documentation: if
#'   variables are already in the same units, standardization may not be necessary).
#' @param standardizeY Logical flag for response variable standardization across
#' observations.
#' @param fitIntercept Logical flag indicating whether intercept term should be fit.
#' @param alphaVals a vector of alpha values to be optimized over  
#' @param lambdaVals a vector of lambda values to be optimized over
#' @param nFoldsForParamSelection the number of cross-validation folds to perform
#' @param nTrainingRuns number of training runs to perform
#' @param minFeatureFrequencyPctl a fractional value (0-1). Features in the x 
#'   percentile defined by this parameter after the removal of features with 
#'   zero weights are retained.        
#' @param cumCorNumPredictors the maximum number of predictors to be returned 
#' @param verbose a boolean, whether debugging information should be displayed 
#' @param useLambdaMin a boolean, whether to use the minimum lambda (lambda.min) 
#'   from cross-validation as the optimum lambda, if FALSE then the largest value 
#'   of lambda such that error is within 1 standard error of the minimum (lambda.1se)
#'   is used.
#' @param useOneStdErrRule Use one standard error rule for model selection, i.e.,
#'   select smallest model (in terms of number of predictors) for which the estimated
#'   error (by cross-validation) is within one standard error of the minimum estimated
#'   error.
#' @param id a optional string identifier for the EN run
#' @param useModelYIntercept A logical value indicating if the model intercept term
#'  should be used for prediction of response values held out in each CV fold.
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
#' \item{cvPredMat} {An nObservations x nRepeats matrix with entry i, j indicating 
#'  cross-validation predicted response for the i-th observation in the j-th
#'  cross-validation repeat.}
#' \item{cvMeanSqErrMat} {An nFolds x nRepeats matrix with entry i, j indicating the 
#'  the mean cross validation error in the i-th fold of the j-th cross-validation
#'  repeat.}
#' \item{cvFoldEnPredictorWts} {A feature x en-run matrix, with cells indicating the weight 
#'  associated with a particular feature in a particular elastic net run.  (A zero 
#'  entry indicates that a given feature was not selected in the run associated with 
#'  the matrix column.)}
#' \item{cvFoldEnPredictorTab} {A data frame with summary information for predictors selected
#'  over cv-fold elastic net runs (selection frequency, average weight, sign of 
#'  predictor weights).}
#' }
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
cvElasticNet <- function(nFolds=10, nRepeats=10, cvFoldIds=NULL,
                        minFeatureResponseAbsCor=0,
                        featureMat, responseVec,
                        standardize=TRUE, standardizeY=FALSE, fitIntercept=TRUE,
                        alphaVals=seq(0.2, 1, length=9), lambdaVals=NULL, 
                        nFoldsForParamSelection=10, nCvRepeats=10,
                        nTrainingRuns=200, minFeatureFrequencyPctl=0.95,
                        cumCorNumPredictors=10, useLambdaMin=TRUE, 
                        useOneStdErrRule=FALSE, id="", 
                        keepEnResults=FALSE, useModelYIntercept=FALSE,
                        verbose=TRUE) {
  if (minFeatureResponseAbsCor > 0){
    featureMat <- selectCorrelatedRows(Y=responseVec, X=featureMat, 
                                       corThreshold = minFeatureResponseAbsCor)
  }
  
  if (is.null(cvFoldIds)){
    cvFoldIds <- getCvFoldIds(nObs=ncol(featureMat), nFolds=nFolds, nRepeats=nRepeats)
  } else{
    nFolds <- max(as.numeric(cvFoldIds))
    nRepeats <- ncol(cvFoldIds)
  }
  
  cvPredMat <- matrix(NA, nrow=ncol(featureMat), ncol=nRepeats)
  rownames(cvPredMat) <- colnames(featureMat)
  cvMeanSqErrMat <- matrix(NA, nrow=nFolds, ncol=nRepeats)
  enRunTimes <- matrix(NA, nrow=nFolds, ncol=nRepeats)
  
  allEnResults <- vector(mode = "list", length = nFolds*nRepeats)
  enResultIds <- vector(mode = "character", length = nFolds*nRepeats)
  
  k <- 0
  for (j in seq_len(nRepeats)){
    for (i in seq_len(nFolds)){
      k <- k + 1
      iTest <- which(cvFoldIds[, j] == i)
      iTrain <- setdiff((1:ncol(featureMat)), iTest)
      
      if (verbose){
        writeLines(paste0("XVal Repeat: ", j, "."))
        writeLines(paste0("XVal Fold: ", i, "."))
        writeLines("Test Set Indices:")
        writeLines(as.character(iTest))
        writeLines("Random Seed:")
        writeLines(as.character(.Random.seed))
        if (require(pryr)){
          writeLines(paste0("cvElasticNet() allEnResults: ", 
                     pryr::object_size(allEnResults) / (10^6), " MB."))
        }
      }
      
      runTimeInfo <- system.time({
      enResults <- elasticNet(featureMat=featureMat[, iTrain, drop=FALSE], 
                              responseVec=responseVec[iTrain],
                              standardize=standardize, 
                              standardizeY=standardizeY, 
                              fitIntercept=fitIntercept,
                              alphaVals=alphaVals, 
                              lambdaVals=lambdaVals, 
                              nFoldsForParamSelection=nFoldsForParamSelection, 
                              nCvRepeats=nCvRepeats,
                              nTrainingRuns=nTrainingRuns, 
                              minFeatureFrequencyPctl=minFeatureFrequencyPctl,
                              cumCorNumPredictors=cumCorNumPredictors, 
                              useLambdaMin=useLambdaMin, 
                              useOneStdErrRule=useOneStdErrRule, 
                              id=paste0(id, "_cvfold_", i, "_repeat_", j), 
                              verbose=verbose)
      })
      enRunTimes[i, j] <- runTimeInfo["elapsed"]
      
      if (!keepEnResults){
        enResults$featureWtMat <- NULL
      }
      
      allEnResults[[k]] <- enResults
      enResultIds[k] <- enResults$id
      
      yPred <- predictWithLinRegModel(model=enResults, 
                                      useModelYIntercept = useModelYIntercept,
                                      newData=featureMat[, iTest, drop=FALSE])
      
      cvPredMat[iTest, j] <- yPred
      cvMeanSqErrMat[i, j] <- mean((responseVec[iTest] - yPred)^2)
    }
  }
  
  cvEnResults <- list()
  cvEnResults$cvPred <- rowMeans(cvPredMat)
  cvEnResults$cvPredR <- cor.test(cvEnResults$cvPred, responseVec)$estimate
  cvEnResults$cvPredRsqared <- (cvEnResults$cvPredR)^2
  cvEnResults$cvMeanSqErr <- mean(as.numeric(cvMeanSqErrMat))
  cvEnResults$cvSdMeanSqErr <- sd(as.numeric(cvMeanSqErrMat))
  cvEnResults$cvPredMat <- cvPredMat
  cvEnResults$cvMeanSqErrMat <- cvMeanSqErrMat
  
  # summarize feature selection over cv folds ---------------------------------
  names(allEnResults) <- enResultIds
  
  # Get the total set of features selected in one or more fold-specific 
  # elastic net runs.
  allEnFeatures <- c(lapply(allEnResults, function(x) {names(x$predictorWts)}),
                           recursive=TRUE)
  allEnFeatures <- unique(unname(allEnFeatures))
  
  # Construct a feature x en-run matrix, with cells indicating the weight associated
  # with a particular feature in a particular elastic net run.  (A zero entry indicates
  # that a given feature was not selected in the run associated with the matrix column).
  cvFoldEnPredictorWts <- matrix(0, nrow=length(allEnFeatures), ncol=length(allEnResults))
  rownames(cvFoldEnPredictorWts) <- allEnFeatures
  colnames(cvFoldEnPredictorWts) <- names(allEnResults)
  for (cvFoldName in names(allEnResults)){
    enPredictorWts <- allEnResults[[cvFoldName]][["predictorWts"]]
    
    stopifnot(all(names(enPredictorWts) %in% rownames(cvFoldEnPredictorWts)))
    
    cvFoldEnPredictorWts[names(enPredictorWts), cvFoldName] <- enPredictorWts
  }
  
  # For a high-level summary, average feature weights along rows (i.e., across
  # all fold-specific elastic net runs).
  cvFoldEnPredictorAvgWt <- rowMeans(cvFoldEnPredictorWts)
  cvFoldEnPredictorAvgWt <- cvFoldEnPredictorAvgWt[order(abs(cvFoldEnPredictorAvgWt),
                                                         decreasing = TRUE)]
  
  # reorder rows of fold-specific predictor matrix by magnitude of avg. weight
  # (over all folds).
  cvFoldEnPredictorWts <- cvFoldEnPredictorWts[names(cvFoldEnPredictorAvgWt), ]
  
  nEnResults <- ncol(cvFoldEnPredictorWts)
  cvFoldEnPredictorFreq <- apply(cvFoldEnPredictorWts, MARGIN = 1,
                                 FUN = function(x) { sum(x != 0) / nEnResults })
  
  cvFoldEnPredictorTab <- data.frame(PREDICTOR=names(cvFoldEnPredictorAvgWt),
                                     AVG_WEIGHT=cvFoldEnPredictorAvgWt,
                                     FREQUENCY=cvFoldEnPredictorFreq,
                                     SIGN=NA,
                                     stringsAsFactors = FALSE)
  rownames(cvFoldEnPredictorTab) <- cvFoldEnPredictorTab$PREDICTOR
  
  for (predName in rownames(cvFoldEnPredictorTab)){
    predWts <- cvFoldEnPredictorWts[predName, ]
    nonZeroPredWts <- predWts[predWts != 0]
    if (all(nonZeroPredWts > 0)){
      cvFoldEnPredictorTab[predName, "SIGN"] <- "+"
    } else if (all(nonZeroPredWts < 0)){
      cvFoldEnPredictorTab[predName, "SIGN"] <- "-"
    } else{
      cvFoldEnPredictorTab[predName, "SIGN"] <- "+/-"
    }
  }
  
  cvEnResults$cvFoldEnPredictorWts <- cvFoldEnPredictorWts
  cvEnResults$cvFoldEnPredictorTab <- cvFoldEnPredictorTab
  #----------------------------------------------------------------------------
  
  if (keepEnResults){
    cvEnResults$allEnResults <- allEnResults
  }
  
  cvEnResults$enRunTimes <- enRunTimes
  
  return(cvEnResults)
}