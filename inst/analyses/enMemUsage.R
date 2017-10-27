#--------------------------------------------------------------------------------------------------
# LOAD ADDITIONAL PACKAGES NEEDED FOR PERFORMANCE ANALYSIS.
#--------------------------------------------------------------------------------------------------
library(pryr)

#--------------------------------------------------------------------------------------------------
# SET UP ELASTIC NET TEST INPUTS.
#--------------------------------------------------------------------------------------------------
# Set random seed
set.seed(1)

library(rcellminer)
library(rcellminerData)
library(rcellminerElasticNet)
library(impute)

molDB <- getMolDataMatrices()
drugActData <- exprs(getAct(rcellminerData::drugData))

# Drug NSC
drug <- "102816"

# Included datasets
includeDatasets <- c("exp", "mut")

# Excluded features and cell lines
excludeCellLines <- c("ME:MDA_N") 
#excludeFeatures <- c("expSLFN11")
excludeFeatures <- NULL

## Restrict gene set 
# load(file.path(.lmp, "gene_set_pathway_analysis", "data", "ddr_genelist.Rdata"))
# restrictedGeneSet <- rownames(ddr.list)
restrictedGeneSet <- NULL

# Show debugging data?
verbose <- TRUE

# THRESHOLDS 
# Features must have absolute correlation greater than corThreshold to be
# included in elastic net.
corThreshold <- 0
# Elastic net results will be augmented to include information about input features 
# that are correlated with selected predictors above predictorCorFeatureThreshold.
predCorFeatureThreshold <- 0.8
minNumResponsiveLines <- 3
responseThreshold <- 0.5

minNumMutations <- 4

# Results parameters
maxPlotPredictors <- 60

filePrefix <- "DATA_exp_mut_NSC_102816"
outputDir <- "~/Downloads/"
outputFilePath <- "enResults.csv"

tic()

#----[prepare drug data]-----------------------------------------------------------------
# Make sure data is available for drug.
stopifnot(all(is.element(drug, rownames(drugActData))))

# Exclude selected cell lines.
drugActData <- drugActData[, -which(colnames(drugActData) %in% excludeCellLines)]

# Exclude cell lines w/o activity data for drug of interest
drugActProfile <- setNames(as.numeric(drugActData[drug, ]), colnames(drugActData))
drugActProfile <- drugActProfile[!is.na(drugActProfile)]
drugActData <- drugActData[, names(drugActProfile)]

# NOTE DRUG ACTIVITY PROFILES WITH AN INSUFFICIENT NUMBER OF RESPONSIVE LINES.
numResponsive <- sum(drugActProfile > responseThreshold, na.rm=TRUE)
if (numResponsive < minNumResponsiveLines){
  filePrefix <- paste0("invalid_", filePrefix)
}

#----[prepare molecular data]------------------------------------------------------------
# Restrict molecular data to drug data cell line set prepared above.
molDB <- lapply(molDB, function(x) {
  x[, names(drugActProfile)]
})

# Restrict molecular data to specified subtypes.
datasetIdx <- (names(molDB) %in% includeDatasets)
molDB <- molDB[datasetIdx]

# Impute missing gene expression data. 
if ("exp" %in% names(molDB)){
  # NOTE: invisible() is used to suppress function print statements
  molDB[["exp"]] <- invisible(impute.knn(molDB[["exp"]])$data)
}

# Transform to z-score data.
if ("swa" %in% names(molDB)){
  molDB[["swa"]] <- t(scale(t(molDB[["swa"]])))
}

# Filter mutation profiles with fewer than minNumMutations.
if ("mut" %in% names(molDB)){
  countMutations <- function(zScoreMutPattern){
    minVal <- min(zScoreMutPattern, na.rm = TRUE)
    maxVal <- max(zScoreMutPattern, na.rm = TRUE)
    midVal <- (minVal + maxVal)/2
    mutCount <- 0
    if (minVal < maxVal){
      mutCount <- sum(zScoreMutPattern > midVal)
    }
    return(mutCount)
  }
  
  mutCounts <- apply(molDB$mut, MARGIN = 1, FUN = countMutations)
  molDB$mut <- molDB$mut[(mutCounts >= minNumMutations), ]
  
  # Transform to z-score data.
  molDB[["mut"]] <- t(scale(t(molDB[["mut"]])))
}

#----[prepare elastic net feature data]--------------------------------------------------
# Filter datasets
featureSet <- selectCorrelatedRowsFromMatrices(Y=drugActProfile, 
                                               XList=molDB, 
                                               corThreshold=corThreshold)

# Filter features 
if(!is.null(restrictedGeneSet)) {
  featureSet <- restrictFeatureMat(geneSet=restrictedGeneSet, featureMat=featureSet)  
}

featureSet <- featureSet[setdiff(rownames(featureSet), excludeFeatures), ]

# debug -----------
# featureSet <- featureSet[1:100, ]
# debug -----------

# Included feature sets
names(molDB)


#--------------------------------------------------------------------------------------------------
# RUN BASIC ELASTIC NET.
#--------------------------------------------------------------------------------------------------
elNetResults <- elasticNet(featureMat=featureSet, responseVec=drugActProfile, id=drug,
                           alphaVals=seq(0.2, 1, length = 10), useLambdaMin=TRUE,
                           cumCorNumPredictors=maxPlotPredictors, verbose=TRUE,
                           standardize = FALSE, useOneStdErrRule = TRUE)

#--------------------------------------------------------------------------------------------------
# RUN CV ELASTIC NET.
#--------------------------------------------------------------------------------------------------
set.seed(2)
cvFoldIds <- getCvFoldIds(nObs = length(drugActProfile), nFolds = length(drugActProfile))

elNetResults$fullEnCvResults <- cvElasticNet( 
  cvFoldIds=cvFoldIds, keepEnResults=FALSE,
  featureMat=featureSet, responseVec=drugActProfile, id=drug,
  alphaVals=seq(0.2, 1, length = 10), useLambdaMin=TRUE,
  cumCorNumPredictors=maxPlotPredictors, verbose=TRUE,
  standardize = FALSE, useOneStdErrRule = TRUE)

#--------------------------------------------------------------------------------------------------
# DETERMINE WHICH DRUGS DO NOT HAVE ELASTIC NET RESULTS.
#--------------------------------------------------------------------------------------------------
require(rcellminer)
require(rcellminerElasticNet)
require(stringr)

# Generated elastic net results ---------------------------------------------------------
enDataDir <- "~/Dropbox/TMP/ENDEV/BIOWULF/"
enResultsFiles <- dataFiles <- list.files(path=enDataDir, pattern = "*.Rdata")

# Expected elastic net results ----------------------------------------------------------
drugAnnot <- getFeatureAnnot(rcellminerData::drugData)[["drug"]]
fdaDrugAnnot <- drugAnnot[which(drugAnnot$FDA_STATUS == "FDA approved"), ]
ctDrugAnnot  <- drugAnnot[which(drugAnnot$FDA_STATUS == "Clinical trial"), ]

fdaCtNscs <- c(rownames(fdaDrugAnnot), rownames(ctDrugAnnot))

dataTypes <- c("exp", "mut", "swa", "exp_mut", "exp_swa", "mut_swa", "exp_mut_swa")

expectedEnResultsFiles <- lapply(setNames(fdaCtNscs, fdaCtNscs), function(nsc) {
  paste0("DATA_", dataTypes, "_NSC_", nsc, ".Rdata")
})
expectedEnResultsFiles <- unname(c(expectedEnResultsFiles, recursive = TRUE))

# Comparison ----------------------------------------------------------------------------
missingEnResults <- NULL
for (enFile in expectedEnResultsFiles){
  invalidEnFile <- paste0("invalid_", enFile)
  if (!(enFile %in% enResultsFiles)){
    if (!(invalidEnFile %in% enResultsFiles)){
      missingEnResults <- c(missingEnResults, enFile)
    }
  }
}

#--------------------------------------------------------------------------------------------------