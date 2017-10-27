#--------------------------------------------------------------------------------------------------
# SET UP
#--------------------------------------------------------------------------------------------------
library(rcellminerElasticNet)
library(stringr)
library(gplots)

enDataDir <- "~/Dropbox/TMP/ENDEV/BIOWULF/"

dataFiles <- list.files(path=enDataDir, pattern = "*.Rdata")

nscSet <- NULL
for (fName in dataFiles){
  tmp <- str_split(string = fName, pattern = "NSC_")[[1]][2]
  nsc <- str_split(string = tmp, pattern = "\\.")[[1]][1]
  if (!(nsc %in% nscSet)){
    nscSet <- c(nscSet, nsc)
  }
}

dataTypes <- c("exp", "mut", "swa", "exp_mut", "exp_swa", "mut_swa", "exp_mut_swa")

cvCorMat <- matrix(NA, nrow=length(nscSet), ncol=length(dataTypes))
rownames(cvCorMat) <- nscSet
colnames(cvCorMat) <- dataTypes

#--------------------------------------------------------------------------------------------------
# FILL CORRELATION MATRIX (FULL ELASTIC CV-PREDICTED VS. ACTUAL CORRELATIONS).
#--------------------------------------------------------------------------------------------------
badFiles <- NULL
for (nsc in nscSet){
  for (dType in dataTypes){
    rdataFile <- paste0("DATA_", dType, "_NSC_", nsc, ".Rdata")
    if (rdataFile %in% dataFiles){
      load(paste0(enDataDir, .Platform$file.sep, rdataFile)) # object name: elNetResults.
      cvPredRsquared <- elNetResults$fullEnCvResults$cvPredRsqared
      if (!is.null(cvPredRsquared)){
        cvCorMat[nsc, dType] <- elNetResults$fullEnCvResults$cvPredR
      } else{
        badFiles <- c(badFiles, rdataFile)
      }
      #rm(elNetResults)
    }
  }
}

#--------------------------------------------------------------------------------------------------
# MAKE HEAT MAP.
#--------------------------------------------------------------------------------------------------



# Note: Possible alternatives to heatmap.2 (gplots):
# --- https://plot.ly/ggplot2/ggdendro-dendrograms/

# Note: Row/Column clustering is currently not done because of NA values in data matrix.

plotTitle <- "[Fill in]"
outputFilePath <- "inst/extdata/en_cvresults_heatmap.pdf"
pdf(file = outputFilePath, width = 8, height = 16)
par(cex.main=0.95)
heatmap.2(na.exclude(cvCorMat),
          #margins = c(40,30),
          #dendrogram="none", 
          #Colv = FALSE, 
          #Rowv = FALSE,
          main=plotTitle,
          density.info="histogram",
          key=TRUE,
          symkey=FALSE,
          trace="none",
          na.color=par("bg"),
          cexRow=0.5,
          cexCol=0.7
          )
dev.off()

#--------------------------------------------------------------------------------------------------