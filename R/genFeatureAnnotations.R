#' Generate feature annotations 
#' 
#' @param molDb molecular data used for elastic net results
#' @param elNetResults elastic net results object 
#' @param responseVec response vector used in the elastic net results
#' @param verbose a boolean, whether to display debugging information
#' @return a vector of strings with feature annotations in the form: cr=CORRELATION 
#'   cm=CUMULATIVE_CORRELATION fq=FEATURE_FREQUENCY
#'   \itemize{
#'     \item{cr} {individual correlation of feature}
#'     \item{cm} {cumulative correlation of features up to this feature}
#'     \item{fq} {frequency of feature in multiple training runs}
#'   }
#'   
#' @concept rcellminerElasticNet
#' @export
#' 
#' @importFrom rcellminer getMolDataType
genFeatureAnnotations <- function(molDb, elNetResults, responseVec, verbose=FALSE) {
  numPlotPredictors <- length(elNetResults$predictorWts)
  featureAnnot <- vector(mode="character", length=numPlotPredictors)
  names(featureAnnot) <- names(elNetResults$predictorWts[1:numPlotPredictors])
  
  i <- 0
  
  for (featureName in names(featureAnnot)){
    i <- i + 1
    if(verbose) {
      cat("featureName: ", featureName, "\n") 
    }  
    
    cumCor <- round(elNetResults$predictedResponseCor[i], digits=2)
    predProfile <- molDb[[getMolDataType(featureName)]][featureName, ]
    
    if(verbose) {
      cat("predProfile Length: ", length(predProfile), "\n")
      cat("responseVec Length: ", length(responseVec), "\n")      
    }
    
    drugPredCor <- round(cor.test(x=predProfile, y=responseVec)$estimate, digits=2)
    
    selectionFreq <- round(elNetResults$predictorSelectionFreq[featureName], digits=2)
    
    featureCaption <- paste("cr=", drugPredCor, " cm=", cumCor, " fq=", selectionFreq, sep="")
    featureAnnot[featureName] <- featureCaption
  }  
  
  return(featureAnnot)
}
