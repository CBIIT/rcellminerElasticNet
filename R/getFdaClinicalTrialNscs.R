#' Get a list of NSCs that are FDA approved on in clinical trials
#' 
#' @param file a file name, if provided a file will be generated at this location with the results
#' 
#' @return a vector of NSCs annotated with an FDA status of FDA approved or 
#'   in clinical trials for rcellminer data
#' 
#' @examples 
#' getFdaClinicalTrialNscs()
#' 
#' @concept rcellminerElasticNet
#' @export
getFdaClinicalTrialNscs <- function(file=NULL) {
  require(rcellminer)
  
  nci60Annot <- rcellminer::getFeatureAnnot(rcellminerData::drugData)[["drug"]]
  fdaCtAnnot <- nci60Annot[which(nci60Annot$FDA_STATUS %in% c("Clinical trial", "FDA approved")), ]
  fdaCtNscs <- fdaCtAnnot$NSC
  
  if(!is.null(file)) {
    write.table(fdaCtNscs, file=file, row.names=FALSE, quote=FALSE, col.names="nsc")
  }
  
  return(fdaCtNscs)
}
