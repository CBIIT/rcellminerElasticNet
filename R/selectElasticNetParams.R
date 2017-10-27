#' Select the optimal elastic net parameters based on cross-validation
#' error estimates.
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
#' @param nFolds The number of folds used in cross-validation.
#' @param nRepeats The number of cross-validation repeats (with different fold
#'  assignments).
#' @param useLambdaMin a boolean, whether to use the minimum lambda (lambda.min) 
#'   from cross-validation as the optimum lambda, if FALSE then the largest value 
#'   of lambda such that error is within 1 standard error of the minimum (lambda.1se)
#'   is used.
#' 
#' @return A list with the following elements:
#' \itemize{
#'  \item{alpha} {The optimal alpha parameter.}
#'  \item{lambda} {The optimal lambda parameter.}
#'  \item{minCvError} {The minimum average cross validation error.}
#'  \item{cvInfoTab} {A data frame summarizing minimum cross validation errors and lambda
#'   parameter selections for each possible choice of the alpha parameter.}
#'  \item{cvFoldIds} {A matrix indicating the fold ids for each cross validation repeat
#'   along the columns.}
#' }
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
selectElasticNetParams <- function(featureMat, responseVec,
                                   standardize=TRUE, standardizeY=FALSE, fitIntercept=TRUE,
                                   alphaVals=seq(0.2, 1, length=9), lambdaVals=NULL, 
                                   nFolds=10, nRepeats=10, useLambdaMin=TRUE, verbose=TRUE){
  if (is.null(lambdaVals)){
    # We run cv.glmnet() once to obtain a sequence of lambda values based on the data.
    if (standardizeY){
      glmnetY <- scale(responseVec)
    } else{
      glmnetY <- responseVec
    }
    lambdaVals <- cv.glmnet(x=t(featureMat), y=glmnetY,
      nfolds=nFolds, standardize=standardize, intercept=fitIntercept)$lambda
  } 
  cvFoldIds <- getCvFoldIds(nObs=ncol(featureMat), nFolds=nFolds, nRepeats=nRepeats)
  
  # We construct a 3-dimensional array, with each (i, j, k) index specifying the mean
  # cross-validation error for the i-th alpha value, the j-th lambda value,
  # and the k-th repeat of the n-fold cross validation procedure.
  cvErrors <- array(0, dim=c(length(alphaVals), length(lambdaVals), nRepeats))
  
  for (k in (1:nRepeats)){
    # Note that with each cross-validation repeat, we *must* use the same
    # fold split when evaluating the various (alpha, lambda) combinations.
    foldIds <- cvFoldIds[, k]
    for (i in seq_along(alphaVals)){
      alpha <- alphaVals[i]
      cvResult <- cv.glmnet(x=t(featureMat), y=glmnetY,
        alpha=alpha, lambda=lambdaVals, foldid=foldIds, 
        standardize=standardize, intercept=fitIntercept)
      cvmVec <- cvResult$cvm
      if (length(cvmVec) < length(lambdaVals)){
        # Handle instance in which glmnet does not consider full set of lambda
        # values (typically near end of path) because model results are not
        # sufficiently different from one value of lambda to the next.
        tmp <- rep(NA, length(lambdaVals))
        tmp[1:length(cvmVec)] <- cvmVec
        cvmVec <- tmp
      }
      cvErrors[i, , k] <- cvmVec
    }
  }
  
  # Now we collapse the 3-array to a 2-D (alpha value x lambda value) array by
  # averaging the values in the z-dimension (cross-validation repeat).
  meanCvErrorsOverReps <- apply(cvErrors, MARGIN=c(1, 2), FUN=mean)
  
  # For each row (alpha parameter setting), we record the minimum (average)
  # cross validation error, and the lambda value which produced it.
  cvInfoTab <- data.frame(ALPHA=alphaVals, MIN_CV_ERROR=NA, LAMBDA_MIN=NA)
  cvInfoTab$MIN_CV_ERROR <- apply(meanCvErrorsOverReps, MARGIN=1, FUN = min, na.rm=TRUE)
  cvInfoTab$LAMBDA_MIN <- apply(meanCvErrorsOverReps, MARGIN=1, FUN=function(x){
    lambdaVals[which.min(x)]
  })
  
  # Finally, we extract and return the alpha and lambda value that yielded the
  # overall minimum average cross validation error (with the later derived from
  # multiple cross-validation repeats with distinct fold splits).
  iOptRow <- which.min(cvInfoTab$MIN_CV_ERROR)
  optEnParams <- list()
  optEnParams$alpha <- cvInfoTab$ALPHA[iOptRow]
  optEnParams$lambda <- cvInfoTab$LAMBDA_MIN[iOptRow]
  optEnParams$minCvError <- cvInfoTab$MIN_CV_ERROR[iOptRow]
  optEnParams$cvInfoTab <- cvInfoTab
  optEnParams$cvFoldIds <- cvFoldIds
  
  return(optEnParams)
}