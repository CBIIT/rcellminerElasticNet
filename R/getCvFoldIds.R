#' Return a matrix of fold assignments for use in cross-validation. 
#' 
#' @param nObs The number of observations in the data set.
#' @param nFolds The number of folds used in cross-validation.
#' @param nRepeats The number of cross-validation repeats (with different fold
#'  assignments).
#' @param id a optional string identifier for the EN run
#' 
#' @return An nObs x nReplicates matrix with integer fold identifiers along
#'  the columns.
#'  
#' @concept rcellminerElasticNet
#' @export
#' 
getCvFoldIds <- function(nObs, nFolds=10, nRepeats=10){
  if (nObs == nFolds){
    # Leave one out.
    return(matrix(1:nObs, nrow=nObs))
  }
  
#   if (!require(cvTools)){
#     stop("Please install cvTools package.")
#   }
  tmp <- rcellminerElasticNet::xValFolds(nObs, K=nFolds, R=nRepeats, type = "random")
  
  cvFoldIds <- matrix(NA, nrow=nObs, ncol=nRepeats)
  for (j in (1:nRepeats)){
    obsIndexPermutation <- tmp$subsets[, j]
    stopifnot(identical(sort(obsIndexPermutation), 1:nObs))
    cvFoldIds[, j] <- tmp$which[obsIndexPermutation]
  }
  
  return(cvFoldIds)
}
