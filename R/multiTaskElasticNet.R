#' The multi-response elastic net (EN) function with multiple training run averaging.
#' 
#' Applies the elastic net regression algorithm to learn a sparse linear model
#' for predicting a response matrix from a set of input feature vectors.
#'
#' @param X n x p matrix with input feature vectors along columns.
#' @param Y n x k matrix of responses to be predicted using sparse
#'   linear combinations of input feature vectors specified in X
#' @param alphaVals a vector of alpha values to be optimized over  
#' @param lambdaVals a vector of lambda values to be optimized over
#' @param nFoldsForParamSelection the number of cross-validation folds to perform
#' @param nTrainingRuns number of training runs to perform
#' @param minFeatureFrequencyPctl a fractional value (0-1). Features in the x 
#'   percentile defined by this parameter after the removal of features with 
#'   zero weights are retained.        
#' @param verbose a boolean, whether debugging information should be displayed 
#' @param useLambda1se a boolean, whether to use the largest value of lambda such that the error 
#'   is within 1 standard error of the minimum (lambda.1se) as the optimum lambda;
#'   if FALSE then lambda.min is used.
#' @param id a optional string identifier for the EN run
#' 
#' @return a list of enResults objects, each with members: 
#' \itemize{
#'   \item{predictorWts} {weights for selected predictors}
#'   \item{predictorSelectionFreq} {the frequency of the selected predictors across EN iterations}
#'   \item{predictedResponse} {cumulative predicted response across the predictors}
#'   \item{predictedResponseCor} {a vector whose nth element gives the correlation of 
#'     the response with the linear combination of the first n predictors, with 
#'     coefficient weights given in predictorWts.}
#'   \item{call} {the command used to call the function with parameters described}
#'   \item{alpha} {the optimized alpha value used}
#'   \item{lambda} {the optimized lambda value used}
#'   \item{foldIds} {the cross-validation IDs used}
#'   \item{cvm} {mean cross-validated errors for alpha values tried}
#'   \item{cvOutput} {an object containing all cross-validation information for the optimized alpha}
#'   \item{featureWtMat} {a matrix with the feature weights for all features across all training runs}
#' }
#' 
#' @author Vinodh Rajapakse 
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
#' @importFrom glmnet cv.glmnet glmnet
multiTaskElasticNet <- function(X, Y, 
                                alphaVals=seq(0.2, 1, length=9), lambdaVals=NULL, nFoldsForParamSelection=4,
                                nTrainingRuns=200, minFeatureFrequencyPctl=0.95, verbose=TRUE, 
                                useLambda1se=TRUE, id="") {
  # ------[ (alpha, lambda) parameter selection ]---------------------------------------------------
  if (verbose){
    cat("Selecting parameters (alpha, lambda).\n")
  }
  
  names(alphaVals) <- as.character(alphaVals)
  xvalResults <- data.frame(bestLambda=numeric(length(alphaVals)),
                            bestLambdaMeanXvalError=numeric(length(alphaVals)),
                            stringsAsFactors=FALSE)
  rownames(xvalResults) <- names(alphaVals)
  
  # We run cv.glmnet() once to obtain set of fixed set of fold assignments to be used over 
  # cross-validation runs with different settings of alpha parameter.
  foldIds <- cv.glmnet(x=X, y=Y, family="mgaussian", nfolds=nFoldsForParamSelection, keep=TRUE)$foldid
  
  cvObjs <- lapply(alphaVals, FUN=function(alpha){ 
    cv.glmnet(x=X, y=Y, family="mgaussian", alpha=alpha, lambda=lambdaVals, foldid=foldIds) })
  
  xvalResults$bestLambda <- vapply(cvObjs, 
                                   FUN=function(x) x[[ifelse(useLambda1se, yes="lambda.1se", no="lambda.min")]],
                                   FUN.VALUE=numeric(1)
  )
  
  xvalResults$bestLambdaMeanXvalError <- vapply(cvObjs, 
                                                FUN=function(x) x[["cvm"]][ which(x[["lambda"]] == x[[ifelse(useLambda1se, yes="lambda.1se", no="lambda.min")]])[1] ],
                                                FUN.VALUE=numeric(1)
  )
  
  iMin <- which.min(xvalResults$bestLambdaMeanXvalError)
  optAlpha <- as.numeric(rownames(xvalResults)[iMin])
  optLambda <- xvalResults$bestLambda[iMin]
  optCvOutput <- cvObjs[[as.character(optAlpha)]]
  # -----------------------------------------------------------------------------------------------
  
  # ------[ fit models and compute feature effect sizes ]------------------------------------------
  if (verbose){
    cat("Fitting elastic net models and computing feature effect sizes.\n")
  }
  
  # The additional row is used to track the beta0 coefficient    
  coefArray <- array(0, dim=c((ncol(X)+1), ncol(Y), nTrainingRuns))
  trainingSetFraction <- 1 - (1/nFoldsForParamSelection)
  trainingSetSize <- floor(trainingSetFraction*nrow(Y))
  
  for (k in seq_len(nTrainingRuns)){
    if (verbose){
      cat(paste("training iteration:", k), sep="\n")
    }
    
    # NOTE: This is not bootstrapping; selection is not done with replacement. 
    selector <- sample(x=seq_len(nrow(Y)), size=trainingSetSize, replace=FALSE)
    
    glmnetOutput <- glmnet(x=X[selector, ], y=Y[selector, ], family="mgaussian", alpha=optAlpha, lambda=lambdaVals)
    
    # YColCoefs is a list of coef. vectors (1 column matrices) for each response along columns of Y.
    YColCoefs <- coef(glmnetOutput, s=optLambda)
    
    coefArray[, , k] <- as.matrix(as.data.frame(lapply(YColCoefs, as.numeric)))
  }
  
  coefNames <- rownames(YColCoefs[[1]])
  YColNames <- names(YColCoefs)
  dimnames(coefArray) <- list(coefNames, YColNames, as.character(seq_len(nTrainingRuns)))
  
  # Extract and average y-intercept ('beta_0') values over training runs.
  avgBeta0 <- as.numeric(apply(coefArray["(Intercept)", , , drop=FALSE], MARGIN=c(1,2), FUN=mean))
  names(avgBeta0) <- dimnames(coefArray)[[2]]
  
  coefArray <- coefArray[setdiff(dimnames(coefArray)[[1]], "(Intercept)"), , ]
  
  # ------[ filter features ]------------------------------------------
  
  # featureSelFreq is a num(coefs) x num(response) array of selection frequencies.
  featureSelFreq <- apply(coefArray, MARGIN=c(1, 2), FUN=function(x) sum(x != 0)/nTrainingRuns)
  # minFeatureSelFreq is a vector of minimum selection frequencies (over responses).
  minFeatureSelFreq <- apply(featureSelFreq, MARGIN=1, FUN=min)
  
  # (1) Exclude features that were not selected for all responses in at least one training run.
  coefArray <- coefArray[names(minFeatureSelFreq[minFeatureSelFreq != 0]), , , drop=FALSE]
  
  # (2) Exclude features whose weights are inconsistent in sign - for any response.
  nonZeroWtsHaveSameSign <- function(x){
    numPos <- sum(x > 0)
    numNeg <- sum(x < 0)
    numNonZero <- sum(x != 0)
    return(abs((numPos - numNeg)/numNonZero) == 1)
  }
  
  featuresHaveConsistentSigns <- apply(apply(coefArray, MARGIN=c(1,2), FUN=nonZeroWtsHaveSameSign), 
                                       MARGIN=1, FUN=all)
  coefArray <- coefArray[which(featuresHaveConsistentSigns), , , drop=FALSE]
  
  # (3) Retain features that were most frequently selected during the training runs.
  minFeatureSelFreq <- minFeatureSelFreq[dimnames(coefArray)[[1]]]
  minFreq <- quantile(minFeatureSelFreq, probs=minFeatureFrequencyPctl)
  
  minFeatureSelFreq <- minFeatureSelFreq[minFeatureSelFreq >= minFreq]
  coefArray <- coefArray[names(minFeatureSelFreq), , , drop=FALSE]
  
  # -----------------------------------------------------------------------------------------------
  
  # ------[ organize outputs ]---------------------------------------------------------------------
  
  # We return a list of enResults objects - one for each reponse vector along a column of Y.
  # Most elements of these enResults objects are common to all responses, and we compute some of
  # these shared components first.  Correlations of the predicted response with each response,
  # are response-specific, and are computed within the loop below.
  
  #----[ construct shared elements of enResults objects ]--------------
  # Note that we form a single vector of average feature weights by (1) averaging the feature
  # weights for each feature-response pair over all training runs and (2) averaging the latter
  # feature-response averages (over responses) to obtain a single feature-specific weight.
  avgFeatureWts <- apply(apply(coefArray, MARGIN=c(1, 2), FUN=mean), MARGIN=1, FUN=mean)
  avgFeatureWts <- avgFeatureWts[order(abs(avgFeatureWts), decreasing=TRUE)]
  avgFeatureSelFreq <- apply(featureSelFreq[names(avgFeatureWts), , drop=FALSE], MARGIN=1, FUN=mean)
  
  avgFeatureWtVsRun <- apply(coefArray, MARGIN=c(1,3), FUN=mean)
  
  augmentedPredMat <- cbind(as.vector(array(1, dim=nrow(Y))), X[, names(avgFeatureWts)])
  rownames(augmentedPredMat) <- rownames(Y)
  colnames(augmentedPredMat) <- c("(Intercept)", names(avgFeatureWts))
  
  numPredictors <- length(avgFeatureWts)
  predictedResponse <- matrix(0, nrow=nrow(Y), ncol=numPredictors)
  rownames(predictedResponse) <- rownames(Y)
  colnames(predictedResponse) <- as.character(seq_len(numPredictors))
  
  for (n in seq_len(numPredictors)){
    restrictedBeta <- vector(mode="numeric", length=ncol(augmentedPredMat))
    names(restrictedBeta) <- colnames(augmentedPredMat)
    restrictedBeta["(Intercept)"] <- mean(avgBeta0)
    restrictedBeta[names(avgFeatureWts[1:n])] <- avgFeatureWts[1:n]
    
    predictedResponse[, as.character(n)] <- as.numeric(augmentedPredMat %*% restrictedBeta)
  }
  #--------------------------------------------------------------------
  
  #----[ construct list of enResults objects ]-------------------------
  enResultsList <- vector(mode="list", length=dim(coefArray)[2])
  names(enResultsList) <- dimnames(coefArray)[[2]]
  
  for (responseName in names(enResultsList)){
    output <- list()
    
    output$predictorWts <- avgFeatureWts
    output$predictorSelectionFreq <- avgFeatureSelFreq
    output$predictedResponse <- predictedResponse
    
    output$predictedResponseCor <- apply(X=predictedResponse, MARGIN=2, 
                                         FUN=function(x) cor.test(x, Y[, responseName])$estimate)
    
    # Get function call
    output$call <- match.call(expand.dots=TRUE)
    
    # Optimal alpha/lambda values
    output$alpha <- optAlpha 
    output$lambda <- optLambda 
    
    # Cross-validation IDs used 
    output$foldIds <- foldIds
    
    # Mean cross-validated errors for alpha values tried 
    output$cvm <- xvalResults$bestLambdaMeanXvalError 
    names(output$cvm) <- alphaVals
    
    # Cross-validation data
    output$cvOutput <- optCvOutput
    
    # The complete feature weight matrix across all training runs
    output$featureWtMat <- avgFeatureWtVsRun
    
    output$id <- id
    
    # Append extra class for print/summary functions
    class(output) <- c("enResults", class(output))
    
    enResultsList[[responseName]] <- output
  }
  #--------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------------
  
  return(enResultsList)
}
