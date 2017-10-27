#' The elastic net (EN) function with multiple training run averaging   
#' 
#' Applies the elastic net regression algorithm to learn a sparse linear model
#' for predicting a response vector from a set of input feature vectors.
#'
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
#' @param obsFractionForModelSelection The fraction of the number of observations.
#'   This is to limit the maximum possible number of predictors considered during 
#'   model selection (default = 0.75).
#' @param id a optional string identifier for the EN run
#' 
#' @return a list with members: 
#' \itemize{
#'   \item{predictorWts} {Coefficient weights for the selected predictors.}
#'   \item{predictorSelectionFreq} {The selection frequency for the final set of predictors, i.e.,
#'    the fraction of training set iterations in which the predictor was selected.}
#'   \item{cvPredRsqVals} {A vector with the k-th entry indicating the square of the Pearson's
#'    correlation coefficient between the actual response and the cross-validation predicted 
#'    response with linear models based on the top k elastic net selected predictors 
#'    (by coeff. weight magnitude).}
#'   \item{cvPredRVals} {A vector with the k-th entry indicating the Pearson's
#'    correlation coefficient between the actual response and the cross-validation
#'    predicted response with linear models based on the top k elastic net selected predictors 
#'    (by coeff. weight magnitude).}
#'   \item{cvMeanSqErrVals} {A vector with the k-th entry indicating the mean cross-validation error
#'    estimate for a model based on the top k predictors (by coeff. weight magnitude).}
#'   \item{cvSdMeanSqErrVals} {A vector with the k-th entry indicating the standard deviation of
#'    the cross-validation error estimate for a model based on the top k predictors 
#'    (by coeff. weight magnitude).}
#'   \item{predictedResponse} {A matrix with the k-th column giving the predicted response with
#'    respect to the top k predictors (by coefficient weight magnitude).}
#'   \item{predictedResponseCor} {a vector whose nth element gives the correlation of 
#'     the response with the linear combination of the first n predictors, with 
#'     coefficient weights given in predictorWts.}
#'   \item{call} {the command used to call the function with parameters described}
#'   \item{alpha} {the optimized alpha value used}
#'   \item{lambda} {the optimized lambda value used}
#'   \item{foldIdsParamSelection} {The cross-validation fold IDs used during elastic net parameter 
#'    selection.}
#'   \item{foldIdsModelSelection} {The cross-validation fold IDs used during final model selection
#'    (selecting among models with different numbers of predictors).}
#'   \item{cvm} {mean cross-validated errors for alpha values tried}
#'   \item{featureWtMat} {A matrix with the feature weights for all features across all training 
#'   runs (with each run based on a different random data subset).}
#' }
#' 
#' @author Vinodh Rajapakse 
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
#' @importFrom glmnet cv.glmnet glmnet
elasticNet <- function(featureMat, responseVec,
                       standardize=TRUE, standardizeY=FALSE, fitIntercept=TRUE,
                       alphaVals=seq(0.2, 1, length=9), lambdaVals=NULL, 
                       nFoldsForParamSelection=10, nCvRepeats=10,
                       nTrainingRuns=200, minFeatureFrequencyPctl=0.95,
                       cumCorNumPredictors=10, useLambdaMin=TRUE, 
                       useOneStdErrRule=FALSE, id="", 
                       obsFractionForModelSelection=0.75,
                       verbose=TRUE) {
  # ------[ (alpha, lambda) parameter selection ]---------------------------------------------------
  if (verbose){
    cat("Selecting parameters (alpha, lambda).\n")
  }
  
  optEnParams <- selectElasticNetParams(featureMat, responseVec, standardize=standardize,
                                        standardizeY=standardizeY, fitIntercept=fitIntercept,
                                        alphaVals=alphaVals, lambdaVals=lambdaVals, 
                                        nFolds=nFoldsForParamSelection, nRepeats=nCvRepeats, 
                                        useLambdaMin=useLambdaMin, verbose=verbose)
  
  optAlpha <- optEnParams$alpha
  optLambda <- optEnParams$lambda
  #optCvOutput <- cvObjs[[as.character(optAlpha)]]
  # -----------------------------------------------------------------------------------------------
  
  # ------[ fit models and compute feature effect sizes ]------------------------------------------
  if (verbose){
    cat("Fitting elastic net models and computing feature effect sizes.\n")
  }
  
  trainingSetFraction <- 1 - (1/nFoldsForParamSelection)
  # The additional row is used to track the beta0 coefficient    
  featureWtMat <- matrix(0, nrow=(nrow(featureMat)+1), ncol=nTrainingRuns)
  trainingSetSize <- floor(trainingSetFraction*length(responseVec))
  
  if(verbose) {
    pb <- txtProgressBar(min=1, max=nTrainingRuns, style=3)    
  }
  
  for (j in (1:nTrainingRuns)) {
    #if (verbose){
    #  cat("training iteration: ")
    #  cat(j)
    #  cat(".\n")
    #}
    
    # NOTE: This is not bootstrapping; selection is not done with replacement 
    selector <- sample(x=(1:length(responseVec)), size=trainingSetSize, replace=FALSE)
    
    if (standardizeY){
      glmnetY <- scale(responseVec[selector])
    } else{
      glmnetY <- responseVec[selector]
    }
    glmnetOutput <- glmnet(x=t(featureMat)[selector, ], y=glmnetY,
      standardize=standardize, intercept=fitIntercept, alpha=optAlpha, lambda=lambdaVals)
    betaAsMat <- predict(glmnetOutput, type="coef", s=optLambda)
    featureWtMat[, j] <- as.numeric(betaAsMat)
    
    if(verbose) {
      setTxtProgressBar(pb, j)
    }
  }
  
  rownames(featureWtMat) <- rownames(betaAsMat)
  colnames(featureWtMat) <- 1:nTrainingRuns
  
  featureNames <- rownames(featureMat)
  stopifnot(identical(featureNames, rownames(featureWtMat)[2:nrow(featureWtMat)]))

  featureNonzeroWtPct <- vector(mode="numeric", length=length(featureNames))
  names(featureNonzeroWtPct) <- featureNames
  
  for (name in featureNames){
    featureNonzeroWtPct[name] <- sum(featureWtMat[name, ] != 0)/nTrainingRuns
  }

  # ------[ filter features ]------------------------------------------

  
  if (all(featureNonzeroWtPct == 0)){
    stop("Check input data (no features could be selected during training runs).")
  }
  # (1) Exclude features that were not selected in any training run.
  reducedFeatureSetNonzeroWtPct <- featureNonzeroWtPct[featureNonzeroWtPct != 0]
  
  # (2) Exclude features whose weights are not consistent in sign.
  omega <- vector(mode="numeric", length=length(reducedFeatureSetNonzeroWtPct))
  names(omega) <- names(reducedFeatureSetNonzeroWtPct)
  for (name in names(omega)){
    numPos <- sum(featureWtMat[name, ] > 0)
    numNeg <- sum(featureWtMat[name, ] < 0)
    numNonzero <- sum(featureWtMat[name, ] != 0)
    omega[name] <- (numPos - numNeg)/numNonzero
  }
  
  reducedFeatureSetNonzeroWtPct <- reducedFeatureSetNonzeroWtPct[names(omega[abs(omega) == 1])]
  
  # (3) Retain features that were most frequently selected during the training runs.
  minFreq <- quantile(reducedFeatureSetNonzeroWtPct, probs=minFeatureFrequencyPctl)
  reducedFeatureSetNonzeroWtPct <- reducedFeatureSetNonzeroWtPct[reducedFeatureSetNonzeroWtPct >= minFreq]
  
  # (4) Order retained features by magnitude of average weight over training runs.
  # NOTE: drop=FALSE below, so that if there is really just one predictor, the
  # subsetting of featureWtMat still returns a matrix, even if it has just one row.
  reducedFeatureSetAvgWts <- rowSums(featureWtMat[names(reducedFeatureSetNonzeroWtPct), , drop=FALSE])/nTrainingRuns
  reducedFeatureSetAvgWts <- reducedFeatureSetAvgWts[order(abs(reducedFeatureSetAvgWts), decreasing=TRUE)]

  reducedFeatureSetNonzeroWtPct <- reducedFeatureSetNonzeroWtPct[names(reducedFeatureSetAvgWts)]
  # -----------------------------------------------------------------------------------------------
  
  # ------[ select best model ]--------------------------------------------------------------------
  
  predictorSet <- names(reducedFeatureSetAvgWts)
  nPredictors <- length(predictorSet)
  cvPredRVals <- setNames(rep(NA, nPredictors), predictorSet)
  cvPredRsqVals <- setNames(rep(NA, nPredictors), predictorSet)
  cvMeanSqErrVals <- setNames(rep(NA, nPredictors), predictorSet)
  cvSdMeanSqErrVals <- setNames(rep(NA, nPredictors), predictorSet)
  
  cvFoldIds <- getCvFoldIds(nObs=ncol(featureMat), nFolds=nFoldsForParamSelection, 
                            nRepeats=nCvRepeats)
  
  maxModelSelectionPredictors <- min(nPredictors, 
                                     floor(length(responseVec) * obsFractionForModelSelection))

  for (i in seq_len(maxModelSelectionPredictors)){
    
    lmCvFit <- tryCatch(getLmCvFit(X=t(featureMat[predictorSet[1:i], , drop=FALSE]), y=responseVec, 
                          nFolds=nFoldsForParamSelection, nRepeats=nCvRepeats, cvFoldIds=cvFoldIds),
                        error = function(x) { return(NULL) },
                        warning = function(x) { return(NULL) })

    if (!is.null(lmCvFit)){
      cvPredRVals[i] <- lmCvFit$cvPredR
      cvPredRsqVals[i] <- lmCvFit$cvPredRsqared
      cvMeanSqErrVals[i] <- lmCvFit$cvMeanSqErr
      cvSdMeanSqErrVals[i] <- lmCvFit$cvSdMeanSqErr
    }
  }
  
  iMinEstErr <- which.min(cvMeanSqErrVals)
  if (length(iMinEstErr) == 1){
    #----[cross-validation results are available]------------------------------
    if (useOneStdErrRule){
      maxErr <- cvMeanSqErrVals[iMinEstErr] + cvSdMeanSqErrVals[iMinEstErr]
      k <- 1
      while (is.na(cvMeanSqErrVals[k]) || (cvMeanSqErrVals[k] >= maxErr)){
        # Because of right condition above, k must be less than iMinEstErr.
        k <- k + 1
      }
      nSelectedPredictors <- k
    } else{
      nSelectedPredictors <- iMinEstErr
    }
    #--------------------------------------------------------------------------
  } else{
    # We get here only if getLmCvFit() fails for every candidate
    # k predictor set.
    nSelectedPredictors <- length(predictorSet)
  }

  # -----------------------------------------------------------------------------------------------
  
  # ------[ organize outputs ]---------------------------------------------------------------------
  output <- list()
  stopifnot(identical(names(reducedFeatureSetAvgWts), names(reducedFeatureSetNonzeroWtPct)))
  output$predictorWts_preModelSelection <- reducedFeatureSetAvgWts
  output$predictorWts <- reducedFeatureSetAvgWts[1:nSelectedPredictors]
  output$yIntercept <- mean(featureWtMat["(Intercept)", ])
  output$predictorSelectionFreq_preModelSelection <- reducedFeatureSetNonzeroWtPct
  output$predictorSelectionFreq <- reducedFeatureSetNonzeroWtPct[1:nSelectedPredictors]
  
  #----[cross-validation results for nested predictor sets]-------------------------
  # The k-th entry in the vectors below gives the result for the model based
  # on the first k predictors (indicated in the names attribute of each vector).
  output$cvPredRVals <- cvPredRVals
  output$cvPredRsqVals <- cvPredRsqVals
  output$cvMeanSqErrVals <- cvMeanSqErrVals
  output$cvSdMeanSqErrVals <- cvSdMeanSqErrVals
  #---------------------------------------------------------------------------------
  
  #----[full training data model predicted vs. actual response correlations]--------
  updatedCumCorNumPredictors <- min(length(output$predictorWts), cumCorNumPredictors)
  
  output$predictedResponse <- matrix(0, nrow=length(responseVec), ncol=updatedCumCorNumPredictors)
  rownames(output$predictedResponse) <- names(responseVec)
  colnames(output$predictedResponse) <- as.character(1:updatedCumCorNumPredictors)
  
  for (n in (1:updatedCumCorNumPredictors)){
    output$predictedResponse[, as.character(n)] <- predictWithLinRegModel(
      coeffVec = output$predictorWts[1:n], 
      yIntercept = output$yIntercept,
      newData = featureMat)
  }
  
  output$predictedResponseCor <- vector(mode="numeric", length=ncol(output$predictedResponse))
  names(output$predictedResponseCor) <- colnames(output$predictedResponse)
  for (j in (1:ncol(output$predictedResponse))){
    output$predictedResponseCor[colnames(output$predictedResponse)[j]] <- cor.test(
      output$predictedResponse[, j], as.numeric(responseVec))$estimate
  }
  #---------------------------------------------------------------------------------
  
  # Get function call
  output$call <- match.call(expand.dots=TRUE)

  # Optimal alpha/lambda values
  output$alpha <- optAlpha 
  output$lambda <- optLambda 

  # Cross-validation fold IDs used 
  output$foldIdsParamSelection <- optEnParams$cvFoldIds
  output$foldIdsModelSelection <- cvFoldIds
  
  # Mean cross-validated errors for alpha values tried 
  output$cvm <- optEnParams$cvInfoTab$MIN_CV_ERROR
  names(output$cvm) <- optEnParams$cvInfoTab$ALPHA
  
  # Cross-validation data
  #output$cvOutput <- optCvOutput
  
  # The complete feature weight matrix across all training runs
  output$featureWtMat <- featureWtMat
  
  output$id <- id
  
  # Append extra class for print/summary functions
  class(output) <- c("enResults", class(output))
  # -----------------------------------------------------------------------------------------------

  return(output)
}




