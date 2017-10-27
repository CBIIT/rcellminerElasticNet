#--------------------------------------------------------------------------------------------------
# rm(list=ls())
# ADJUST TO ELASTIC NET EXAMPLES DIRECTORY PATH
#setwd("/Users/rajapaksevn/Dropbox/lmp_nci/ElasticNetExps/examples/YPREQS/")

library("impute")

source("/Users/lunaa/Dropbox/workspace/lmp_nci/ElasticNetExps/src/ElasticNet.R")
source("/Users/lunaa/Dropbox/workspace/lmp_nci/Utils/src/GenUtils.R")

# LOAD PROTEIN DATA AND FILTER--------------------------------------------------

source("/Users/lunaa/Dropbox/workspace/lmp_nci/NCI60GlobalProteomeComparison/R/amin.reformat.R")

files.dir <- "/Users/lunaa/Dropbox/drug_target_tmp/nci60_global_proteomic_data/amin.satoshi.proteomics/" 
raw.amin <- paste(files.dir, "amin.minus.big.xls/Table_S3.txt", sep="")
dir.amin <- dirname(raw.amin)
amin.file <- paste(dir.amin,"amin.reformat",sep="/")
mat <- amin.reformat.genenames(raw.amin, FALSE)
amin <- amin.reformat.celllinenames(mat,tissues)

# Fix special cases 1:(http://dtp.nci.nih.gov/docs/misc/common_files/NCI-ADRres.html)
colnames(amin) <- sub("BR:MCF7ADR", "OV:NCIADRRES", colnames(amin))
colnames(amin) <- sub("BR:MDAMB435", "ME:MDAMB435", colnames(amin))

amin.select <- amin.select.IBAQ.LFQ(amin,IBAQ=FALSE)

# LOAD NCI MOLDB ---------------------------------------------------------------

load("/Users/lunaa/Dropbox/workspace/lmp_nci/RData/moaToCompoundListMap2.Rdata")
load("/Users/lunaa/Dropbox/NCILMP_DATA/lmpdbPub.Rdata")
drugDB <- pubDrugDB

molDB <- elNetMolDataNCI60
molDB$msp <- amin.select

drugActData <- drugDB$act

# FIX PROTEIN DATA CELL LINE NAMES ---------------------------------------------

abbrCellLinesSimp1 <- gsub("[^A-Z0-9]", "", colnames(amin.select))
abbrCellLinesSimp2 <- gsub("[^A-Z0-9]", "", colnames(molDB[[1]]))
cellLineMatch <- match(abbrCellLinesSimp1, abbrCellLinesSimp2)

colnames(amin.select) <- colnames(molDB[[1]])[cellLineMatch]
rownames(amin.select) <- paste0("msp", rownames(amin.select))

molDB$msp <- as.matrix(amin.select)

#-----[SPECIFY DRUG SUBSET OVER WHICH BATCH ELASTIC NET WILL BE RUN]-------------------------------
# (Raltitrexed, Cytarabine, Nutlin, Bleomycin, Topotecan, Vemurafenib)
#drugSet <- c("639186", "63878", "756876", "125066", "609699", "757438")
#drugSet <- c("609699")

drugSet <- sample(unlist(moaToCompoundListMap2), 10)
stopifnot(all(is.element(drugSet, rownames(drugActData))))
verbose <- TRUE

#-----[SET ELASTIC NET PARAMETERS]-----------------------------------------------------------------
corThreshold <- 0.3
minNumResponsiveLines <- 3
responseThreshold <- 0.5
# Excluded features and cell lines
excludedCellLines <- "ME.MDA_N"
excludeFeatures <- c("mir", "pro", "mda")
outputFilePath <- "TOPOTECAN_Alpha_6_8_10_copexpmut.csv"
#--------------------------------------------------------------------------------------------------
# BATCH ELASTIC NET
#--------------------------------------------------------------------------------------------------
# SET UP RESULT STORAGE
numDrugs <- length(drugSet)
numFeaInResults <- 5
resultTypes <- c("Name", "Wt", "CumCor")

# Make empty data.frame of the right size
batchElNetResults <- as.data.frame(matrix(NA, nrow=numDrugs, ncol=numFeaInResults*length(resultTypes)), stringsAsFactors=FALSE)

tmpNames <- NULL
colCnt <- 1

for(i in 1:numFeaInResults) {
  for(j in 1:length(resultTypes)) {
    tmpName <- paste("Feature", i, resultTypes[j], sep="")
    tmpNames <- c(tmpNames, tmpName)
    
    # Only the name should be character
    if(j %% length(resultTypes) == 1) {
      batchElNetResults[,colCnt] <- vector(mode="character", length=numDrugs)
    } else {
      batchElNetResults[,colCnt] <- vector(mode="numeric", length=numDrugs)      
    }  
  }  
  
  colCnt <- colCnt + 1
}

names(batchElNetResults) <- tmpNames
rownames(batchElNetResults) <- drugSet

# RUN ELASTIC NET---------------------------------------------------------------
tic()
set.seed(1)
curDrug <- 0
for (drug in rownames(batchElNetResults)){
  if (verbose){
    curDrug <- curDrug + 1
    cat("processing drug: ", curDrug, "\n")
  }
  drugActProfile <- as.numeric(drugActData[drug, ])
  names(drugActProfile) <- colnames(drugActData)
  
  # SKIP DRUG ACTIVITY PROFILES WITH AN INSUFFICIENT NUMBER OF RESPONSIVE LINES.
  numResponsive <- sum(drugActProfile > responseThreshold, na.rm=TRUE)
  if (numResponsive < minNumResponsiveLines){
    next
  }
  
  # EXCLUDE CELL LINES WITH MISSING ACTIVITY DATA OR EXCLUDED BY USER 
  # AND RECORD VALID CELL LINES TO ACCESS CORRESPONDING MOLECULAR DATA.
  drugActProfile <- drugActProfile[!is.na(drugActProfile) & !(names(drugActProfile) %in% excludedCellLines)]
  validCellLines <- names(drugActProfile)
  
  featureSets <- list()
  featureNames <- NULL

  # Find features within each feature set that correlate above the threshold 
  # and the feature names to the featureSet
  for(fea in names(molDB)) {
  	if(!(fea %in% excludeFeatures)) {
  		featureSets[[fea]] <- selectCorrelatedRows(y=drugActProfile, X=molDB[[fea]][, validCellLines], corThreshold=corThreshold)
  		
  		if (nrow(featureSets[[fea]]) > 0) {
  			featureNames <- c(featureNames, rownames(featureSets[[fea]]))
  		}
  	}
  }
  
  # Construct the final featureSet 
  featureSet <- matrix(0, nrow=length(featureNames), ncol=ncol(featureSets[[1]]))
  rownames(featureSet) <- featureNames
  colnames(featureSet) <- colnames(featureSets[[1]])

  for(fea in names(molDB)) {
  	if(!(fea %in% excludeFeatures)) {
  		if(nrow(featureSets[[fea]]) > 0) {
  			featureSet[rownames(featureSets[[fea]]), ] <- featureSets[[fea]]
  		}
  	}
  }
  
  try({
    #elNetResults <- elasticNet(featureMat=featureSet, responseVec=drugActProfile, alphaVals=seq(0.6, 0.8, length=10))
    # ADJUSTED CALL TO RECORD CUMULATIVE CORRELATIONS FOR EXPANDED SET OF PREDICTORS.  
    elNetResults <- elasticNet(featureMat=featureSet, responseVec=drugActProfile, alphaVals=seq(0.6, 0.8, length=10), cumCorNumPredictors=34)
    
    save(elNetResults, file=paste("./", drug, ".Rdata", sep=""))
  
    for(i in 1:numFeaInResults) {
      for(j in 1:length(resultTypes)) {
        if(resultTypes[j] == "Name") {
          batchElNetResults[drug, paste("Feature", i, resultTypes[j], sep="")] <- names(elNetResults$predictorWts[i])        
        }
        
        if(resultTypes[j] == "Wt") {
          batchElNetResults[drug, paste("Feature", i, resultTypes[j], sep="")] <- elNetResults$predictorWts[i]       
        }
        
        if(resultTypes[j] == "CumCor") {
          batchElNetResults[drug, paste("Feature", i, resultTypes[j], sep="")] <- elNetResults$predictedResponseCor[i]       
        }
      }
    }
  })
}

toc()

write.csv(batchElNetResults, file=outputFilePath)

# END OF BATCH ELASTIC NET.
#--------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
# PLOT RESULTS
#--------------------------------------------------------------------------------------------------
setwd("/Users/rajapaksevn/Dropbox/lmp_nci/ElasticNetExps/examples/FGS_MANUSCRIPT/")

source("/Users/rajapaksevn/Dropbox/lmp_nci/ElasticNetExps/src/plotElasticNet.R")
source("/Users/rajapaksevn/Dropbox/lmp_nci/ElasticNetExps/src/ElasticNet.R")

load("/Users/rajapaksevn/VNR/DATA/LMPDATA/lmpdb.Rdata")

molDB <- elNetMolDataNCI60
drugActData <- drugDB$act

# (Raltitrexed, Cytarabine, Nutlin, Bleomycin, Topotecan)
#drugSet <- c("639186", "63878", "756876", "125066", "609699")

#drugSet <- sapply( list.files("./RESULTS/"), FUN=function(s) {strsplit(x=s, split="\\.")[[1]][1]} )

drugSet <- c("609699")
maxPlotPredictors <- 50

# outputFileLines <- vector(mode="character", length=(length(drugSet) + 1))
# names(outputFileLines) <- c("NSC", drugSet)
# outputFileLines[1] <- "NSC"
# for (i in (1:maxPlotPredictors)){
#   predHeader <- paste("Pred", i, sep="")
#   weightHeader <- paste("Weight", i, sep="")
#   corHeader <- paste("Cor", i, sep="")
#   cumCorHeader <- paste("CumCor", i, sep="")
#   selectFreqHeader <- paste("SelectFreq", i, sep="")
#   outputFileLines[1] <- paste(outputFileLines[1], predHeader, weightHeader, corHeader, cumCorHeader, selectFreqHeader, sep="\t")
# }

for (nsc in drugSet){
  try(
{
  drugProfile <- as.numeric(drugActData[nsc, ])
  load(paste("./RESULTS/", nsc, ".Rdata", sep=""))
  
  elNetPlotFile <- paste("./PLOTS/", nsc, ".elNetPlot.pdf", sep="")
  
  numPlotPredictors <- min(length(elNetResults$predictorWts), maxPlotPredictors)
  featureAnnot <- vector(mode="character", length=numPlotPredictors)
  names(featureAnnot) <- names(elNetResults$predictorWts[1:numPlotPredictors])
  
  i <- 0
  iPred <- 0
  ln <- nsc
  
  for (featureName in names(featureAnnot)){
    i <- i + 1
    iPred <- iPred + 1
    cumCor <- round(elNetResults$predictedResponseCor[i], digits=2)
    predProfile <- molDB[[getMolDataType(featureName)]][featureName, ]
    drugPredCor <- round(cor.test(x=predProfile, y=drugProfile)$estimate, digits=2)
    selectionFreq <- round(elNetResults$predictorSelectionFreq[featureName], digits=2)
    
    featureCaption <- paste("cr=", drugPredCor, " cm=", cumCor, " fq=", selectionFreq, sep="")
    featureAnnot[featureName] <- featureCaption
    
    ln <- paste(ln, featureName, elNetResults$predictorWts[featureName], drugPredCor, cumCor, selectionFreq, sep="\t")
  }
  
  if (iPred < maxPlotPredictors){
    for (k in (iPred:maxPlotPredictors)){
      ln <- paste(ln, "NA", "NA", "NA", "NA", "NA", sep="\t")
    }
  }
  
  #outputFileLines[nsc] <- ln
  
  plotElasticNet(drugName=nsc, weights=elNetResults$predictorWts[1:numPlotPredictors], 
                 drugAct=drugActData, 
                 molDb=molDB, 
                 numCellLines=ncol(molDB[[1]]), 
                 numFeatures=numPlotPredictors,
                 featureAnnotations=featureAnnot,
                 showCellLineLabels=TRUE, 
                 pdfFilename=elNetPlotFile, thresholdValues=TRUE)
}
  )
}


#writeLines(text=outputFileLines, con="RaltitrexedEtcElNetResults.txt")

#--------------------------------------------------------------------------------------------------
# WRITE OUT RAW DATA
#--------------------------------------------------------------------------------------------------
setwd("/Users/rajapaksevn/Dropbox/lmp_nci/ElasticNetExps/examples/FGS_MANUSCRIPT/")

load("/Users/rajapaksevn/VNR/DATA/LMPDATA/lmpdb.Rdata")

molDB <- elNetMolDataNCI60
drugActData <- drugDB$act

# (Raltitrexed, Cytarabine, Nutlin, Bleomycin, Topotecan)
drugSet <- c("639186", "63878", "756876", "125066", "609699")

for (nsc in drugSet){
  fileName <- paste(nsc, ".elnetdata.csv", sep="")
  load(paste("./RESULTS/", nsc, ".Rdata", sep=""))
  
  rowNames <- c(nsc, names(elNetResults$predictorWts))
  elnetData <- matrix(0, nrow=length(rowNames), ncol=ncol(drugActData))
  rownames(elnetData) <- rowNames
  colnames(elnetData) <- colnames(drugActData)
  
  elnetData[nsc, ] <- as.numeric(drugActData[nsc, ])
  
  for (predName in names(elNetResults$predictorWts)){
    elnetData[predName, ] <- molDB[[getMolDataType(predName)]][predName, ]
  }
  
  elnetData <- as.data.frame(elnetData)
  
  write.csv(x=elnetData, file=fileName, quote=FALSE)
}

#--------------------------------------------------------------------------------------------------