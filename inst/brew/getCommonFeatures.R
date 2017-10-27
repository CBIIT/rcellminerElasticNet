library(rcellminer)

molDB <- getMolDataMatrices()

dataSets <- c("swa", "exp", "mut")

molDB <- molDB[dataSets] 

molDBFeatures <- lapply(molDB, function(x) {unname(removeMolDataType(rownames(x)))})
commonFeatures <- Reduce(intersect, molDBFeatures)

write.table(commonFeatures, file=file.path(.lmp, "dockerfiles", "swathAnalysis", "biowulf", "commonFeatures.txt"), 
            quote=FALSE, row.names=FALSE, col.names="commonFeatures")
