library(rcellminer)
library(rcellminerElasticNet)

# SET PARAMETERS ----
setwd("~/Downloads")

drug <- "756717"
maxPlotPredictors <- 60
excludeCellLines <- c("ME:MDA-N")
pdfFilenameWide <- "~/Downloads/wide.pdf"
pdfFilenameTall <- "~/Downloads/tall.pdf"

# LOAD RCELLMINER DATA ----
molDB <- getMolDataMatrices()
drugActData <- exprs(getAct(rcellminerData::drugData))

# LOAD EN RESULTS ----
swathRData <- file.path("~/Downloads", "swathOutputData", "swathOutput_041816", "DATA_exp_mut_swa_NSC_756717.Rdata")
load(swathRData)

# SETUP DRUG DATA ----
drugActData <- drugActData[, -which(colnames(drugActData) %in% excludeCellLines)]

drugActProfile <- setNames(as.numeric(drugActData[drug, ]), colnames(drugActData))
drugActProfile <- drugActProfile[!is.na(drugActProfile)]
drugActData <- drugActData[, names(drugActProfile)]

# SETUP MOL DATA ----
molDB <- lapply(molDB, function(x) {
  x[, names(drugActProfile)]
})

# MAKE PLOT ----
numPlotPredictors <- min(length(elNetResults$predictorWts), maxPlotPredictors)

featureAnnot <- genFeatureAnnotations(molDB, elNetResults, as.numeric(drugActData[drug, ]))
featureAnnot <- featureAnnot[1:numPlotPredictors]

plotElasticNet(drugName = drug, weights = elNetResults$predictorWts[1:numPlotPredictors], 
               drugAct = drugActData, molDb = molDB, numCellLines = length(drugActProfile), 
               numFeatures = numPlotPredictors, featureAnnotations = NULL, showCellLineLabels = TRUE, 
               pdfFilename = pdfFilenameWide, thresholdValues = TRUE, topMarWt = 2.4, bottomMarWt = 1.5, 
               leftMarHeatmap = 6, verbose = FALSE, pdfWidth = 30)

plotElasticNet(drugName = drug, weights = elNetResults$predictorWts[1:numPlotPredictors], 
               drugAct = drugActData, molDb = molDB, numCellLines = length(drugActProfile), 
               numFeatures = numPlotPredictors, featureAnnotations = NULL, showCellLineLabels = TRUE, 
               pdfFilename = pdfFilenameTall, thresholdValues = TRUE, topMarWt = 2.4, bottomMarWt = 1.5, 
               leftMarHeatmap = 6, verbose = FALSE, pdfWidth = 10)


