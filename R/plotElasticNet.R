#' Plot elastic net results
#' 
#' Creates a plot for elastic net results similar to the GDSC plots (http://www.cancerrxgene.org/help/)
#' 
#' @param drugName a string with the name of drug
#' @param weights a named vector of values with the elastic net weights for a set of features
#' @param drugAct a matrix of drug activity values, a value for each drug tested for each cell line
#' @param molDb a list of matrix of supported data types. The list must use the following prefixes: 
#'   exp: gene expression, mir: miRNA, mut: mutations, pro: protein, mda: metadata, 
#'   cop: copy number variations
#' @param numCellLines an integer of the number of cell lines to be plotted. This function will take the top (by drug activity)
#'     numCellLines/2 and the bottom numCellLines/2 for plotting
#' @param numFeatures an integer of the number of features to be plotted; default 10
#' @param featureAnnotations a vector of optional annotations that will be appended to feature names in parentheses
#' @param showCellLineLabels a boolean whether to show cell line labels; default TRUE
#' @param pdfFilename a string filename for saving the image to a PDF
#' @param thresholdValues a boolean to threshold values to an upper limit of 3 and a lower bound of -3 
#' @param isMutBin a boolean whether the mutation values are binarized (default: TRUE)
#' @param axisScaling a scaling factor for the the feature labels, if not use a 
#'   default will be used depending on if the there annotations
#' @param bottomMarWt a user selected top margin for the the feature weights (default: 2.75)
#' @param topMarWt a user selected top margin for the the feature weights (default: 4)
#' @param leftMarHeatmap the left margin for the the heatmap and drug activity; 
#'   can be increased to handle longer labels (default: 5)
#' @param pdfWidth the width of the PDF (default: 10)
#' @param responseLabel a string describing the response vector (default: "Drug Activity")
#' @param verbose a boolean, whether to display debugging information
#' 
#' @return if verbose is TRUE then a list with the plotted drug activity values, heatmap values, and 
#'   cell line order will be returned 
#' 
#' @details This function uses layout() and may conflict if embedded in other plots using layout().
#'   If the pdfFilename argument is used then the image will not be displayed on the screen.
#' 
#' @return plotElasticNet does not return values; a plot is produced 
#'
#' @note Make sure the numFeatures is less than the length of the weights
#' 
#' @author Augustin Luna
#' 
#' @references GDSC Website: http://www.cancerrxgene.org/help/
#' @examples 
#' #plotElasticNet("PD-0325901", weights, binMutCCLE, drugActCCLE, genExpCCLE, 400, 10)
#' 
#' @concept rcellminerElasticNet
#' @export
#' 
#' @importFrom gplots barplot2
#' @importFrom gdata startsWith
plotElasticNet <- function(drugName, weights, drugAct, molDb, numCellLines, numFeatures=10, 
                           featureAnnotations=NULL, showCellLineLabels=TRUE, 
                           pdfFilename=NULL, thresholdValues=FALSE, isMutBin=TRUE, 
                           axisScaling=NULL, topMarWt=4, bottomMarWt=2.75, 
                           leftMarHeatmap=5, pdfWidth=10, responseLabel="Drug Activity",
                           verbose=TRUE) {
  # FIXME: Force R to not show scientific notation in axis
  options("scipen"=0)
  
  if (is.data.frame(drugAct)){
    drugAct <- as.matrix(drugAct)
  }
  
  ### START TESTING PARAMETERS
  #numCellLines <- 400
  #numFeatures <- 10 
  # Reverse so that the highest weighted feature appears at top of diagram
  #weights <- rev(testFea[1:numFeatures])
  #drugName <- "PD-0325901"
  #mutFea <- binMutCCLE, 
  #drugAct <- drugActCCLE, 
  #genExp <- genExpCCLE
  ### END TESTING PARAMETERS
  
  ### START TESTING PARAMETERS
  # 	numCellLines <- 60
  # 	numFeatures <- 5 
  # 	weights <- elNetResults$nonzeroFeatureWts[1:5]
  # 	drugName <- "609699"
  # 	mutFea <- molDB$mut
  # 	newNames <- as.vector(sapply(rownames(mutFea), function(x) {paste("mut", x, sep="")}))
  #     rownames(mutFea) <- newNames    
  # 	drugAct <- drugActNCI60 
  # 	
  # 	genExp <- molDB$exp
  # 	newNames <- as.vector(sapply(rownames(genExp), function(x) {paste("exp", x, sep="")}))
  #     rownames(genExp) <- newNames	
  #   
  #     copNum <- molDB$copy
  # 	newNames <- as.vector(sapply(rownames(copNum), function(x) {paste("cop", x, sep="")}))
  #     rownames(copNum) <- newNames	
  ### END TESTING PARAMETERS
  
  # Extract predictor types if a full predictor database given ---------------------
  if(!is.null(molDb)) {
    if("mir" %in% names(molDb)) {
      mirExp <- molDb[["mir"]]
    } else {
      mirExp <- NULL
    }
    
    if("cop" %in% names(molDb)) {
      copNum <- molDb[["cop"]]
    } else {
      copNum <- NULL
    }
    
    if("exp" %in% names(molDb)) {
      genExp <- molDb[["exp"]]
    } else {
      genExp <- NULL
    }

	  if("pro" %in% names(molDb)) {
		  proVal <- molDb[["pro"]]
	  } else {
      proVal <- NULL
	  }
    
    if("swa" %in% names(molDb)) {
      swaVal <- molDb[["swa"]]
    } else {
      swaVal <- NULL
    }
	
	  if("mda" %in% names(molDb)) {
		  mdaVal <- molDb[["mda"]]
	  } else {
      mdaVal <- NULL
	  }
	
    if("mut" %in% names(molDb)) {
      mutFea <- molDb[["mut"]]
    } else {
      mutFea <- NULL
    }
  }
  
  # Reverse weights so the highest weight is on top
  weights <- rev(weights)
  
  # Reverse feature annotations so they correspond to the weights correctly
  featureAnnotations <- rev(featureAnnotations)
  
  # Get rows with mutations
  mutIdx <- grep("^mut", names(weights))
  if(length(mutIdx) != 0) {
    mutRows <- mutFea[names(weights)[mutIdx],]    
  } else {
    mutRows <- NULL
  }
  
  # Get values from the drug activity data -------------------------------------
  # Account for odd numbers of cell lines
  if(numCellLines %% 2 == 0) {
    topCellLinesValues <- sort(drugAct[drugName, ], decreasing=TRUE)[1:(numCellLines/2)]    
    topCellLinesNames <- names(sort(drugAct[drugName, ], decreasing=TRUE)[1:(numCellLines/2)])
  } else {
    topCellLinesValues <- sort(drugAct[drugName, ], decreasing=TRUE)[1:((numCellLines/2)+1)]  
    topCellLinesNames <- names(sort(drugAct[drugName, ], decreasing=TRUE)[1:((numCellLines/2)+1)])
  }

  bottomCellLinesValues <- rev(sort(drugAct[drugName, ], decreasing=FALSE)[1:(numCellLines/2)])
  bottomCellLinesNames <- rev(names(sort(drugAct[drugName, ], decreasing=FALSE)[1:(numCellLines/2)]))
  
  maxCellLineValue <- round(max(topCellLinesValues), 2) 
  minCellLineValue <- round(min(bottomCellLinesValues), 2)
  midCellLineValue <- round((max(topCellLinesValues) + min(bottomCellLinesValues)) / 2, 2)
  
  cellLineOrder <- c(topCellLinesNames, bottomCellLinesNames)
  
  # Initialize drug activity and expression data
  heatMapValues <- matrix(NA, ncol=numCellLines, nrow=numFeatures)
  # The unione of the top and bottom cell lines may not always equal all the cell lines togethers
  drugActValues <- drugAct[drugName, cellLineOrder]	
  
  if(verbose) {
    str(topCellLinesNames)
    cat("heatMapValues: ncol: ", ncol(heatMapValues), " nrow: ", nrow(heatMapValues), "\n")
    cat("numCellLines: ", numCellLines, " cellLineOrder: ", length(cellLineOrder), "\n")
  }
  
  # Put data into heatmap matrix
  cnt <- 1
  
  # Keep track of not mutation features 
  nonMutFea <- NULL
  
  for(feature in names(weights)) {
    if(startsWith(feature, "mut")) {
      if(verbose) cat("Mut0: ", feature, "\n")
      heatMapValues[cnt,] <- mutFea[feature, cellLineOrder]
    }
    
    if(startsWith(feature, "exp")) {      
      tmp <- genExp[feature, cellLineOrder]
      
      if(verbose) cat("GeneExp: NCol: ", length(tmp), " NRow: ", nrow(tmp), "\n")
      
      heatMapValues[cnt,] <- genExp[feature, cellLineOrder]
      nonMutFea <- c(nonMutFea, cnt)
    }
    
    if(startsWith(feature, "cop")) {
      heatMapValues[cnt,] <- copNum[feature, cellLineOrder]
      nonMutFea <- c(nonMutFea, cnt)
    }
    
    if(startsWith(feature, "mir")) {
      heatMapValues[cnt,] <- mirExp[feature, cellLineOrder]
      nonMutFea <- c(nonMutFea, cnt)
    }

  	if(startsWith(feature, "pro")) {
  		heatMapValues[cnt,] <- proVal[feature, cellLineOrder]
  		nonMutFea <- c(nonMutFea, cnt)
  	}
  	
    if(startsWith(feature, "swa")) {
      # Scale values before displaying 
      heatMapValues[cnt,] <- scale(swaVal[feature, cellLineOrder])
      nonMutFea <- c(nonMutFea, cnt)
    }
    
  	if(startsWith(feature, "mda")) {
  		heatMapValues[cnt,] <- mdaVal[feature, cellLineOrder]
  		nonMutFea <- c(nonMutFea, cnt)
  	}
	
    if(verbose) {
      cat("Fea0: ", feature, 
          " Max: ", round(max(heatMapValues[cnt,]), 2), 
          " Min: ", round(min(heatMapValues[cnt,]), 2), "\n")
    }
    
    cnt <- cnt + 1
  }
  
  # Threshold the extreme limits of the heatmap values 
  if(thresholdValues) {
    heatMapValues[which(heatMapValues > 3)] <- 3 
    heatMapValues[which(heatMapValues < -3)] <- -3 
  }
  
  # The if statement checks whether all shown features are mutations
  if(length(setdiff(1:numFeatures, nonMutFea)) == numFeatures) {
    
    # Find the extremes of the mutation values in the heatmap
    maxHeatMapValues <- round(max(heatMapValues), 2) 
    minHeatMapValues <- round(min(heatMapValues), 2)
    midHeatMapValues <- round((max(heatMapValues) + min(heatMapValues)) / 2, 2)
  } else {
    
    # Find the extremes of the expression and copy number values in the heatmap
    maxHeatMapValues <- round(max(heatMapValues[nonMutFea,]), 2) 
    minHeatMapValues <- round(min(heatMapValues[nonMutFea,]), 2)
    midHeatMapValues <- round((max(heatMapValues[nonMutFea,]) + min(heatMapValues[nonMutFea,])) / 2, 2)
  }
  
  if(verbose) {
    cat("maxHeatMapValues: ", maxHeatMapValues, 
        " minHeatMapValues: ", minHeatMapValues, "\n")
  }
  
  # Set mutation values to expression or copy number extremes
  cnt <- 1
  for(feature in names(weights)) {
    if(verbose) {
      cat("Fea0: ", feature, "\n")
    }
    
    if(isMutBin && startsWith(feature, "mut")) {
      if(verbose) {
        cat("Mut: ", feature, 
            " Max: ", maxHeatMapValues, 
            " Min: ", minHeatMapValues, "\n")
      }
      
      # Final the middle value for each mutation. This value will be used to binarize.
      tmpMax <- max(mutFea[feature, cellLineOrder])
      tmpMin <- min(mutFea[feature, cellLineOrder])
      tmpMid <- (tmpMax + tmpMin) / 2
      
      heatMapValues[cnt, 
                    which(mutFea[feature, 
                                 cellLineOrder] >= tmpMid)] <- maxHeatMapValues
      
      heatMapValues[cnt, 
                    which(mutFea[feature, 
                                 cellLineOrder] < tmpMid)] <- minHeatMapValues
    } 
    
    cnt <- cnt + 1
  }
  
  if(!is.null(pdfFilename)) {
    pdf(pdfFilename, width=pdfWidth, height=7)
  }
  
  ### Initiate Graphic  
  # NOTE: To make longer gene names fit use cex.axis 
  # in Heatmap section to scale down font or MAYBE increase the 
  # left margin c(BOTTOM, LEFT, TOP, RIGHT) in Heatmap and GI50
  m <- cbind(
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3),
    
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3),
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    c(1,1,1,1,1,1,1,1,1,1,1,1,3,3,3,3), 
    c(2,2,2,2,2,2,2,2,2,2,2,2,4,4,5,5), 
    
    c(2,2,2,2,2,2,2,2,2,2,2,2,4,4,5,5)
  )
  layout(t(m))
  
  # Test layout
  #layout.show(5)
  
  # Set outer margin values 
  # Resize depending on whether annotation is present
  if(!is.null(featureAnnotations)) {
    par(las=1, oma=c(1,6,1,1))
  } else {
    par(las=1, oma=c(1,1,1,1))
  }
  
  heatmap.colors <- colorRampPalette(c("blue", "white", "red"))
  drugact.colors <- colorRampPalette(c("cyan", "white", "magenta"))
  
  ### Heatmap
  #par(mar=c(3,5,4,1))
  par(mar=c(3,leftMarHeatmap,4,1))
  #mat <- matrix(sample(-10:10, 200, replace=T), ncol=20)
  #labels <- rep("XXXXX_MUT", 10)
  
  image(t(heatMapValues), col=heatmap.colors(12), axes=FALSE)
  title(main=paste("Drug Compound: ", drugName, sep="")) 
  
  labels <- names(weights[1:numFeatures])
  
  # Add annotations to labels
  # Resize label depending on whether annotation is present
  if(!is.null(featureAnnotations)) {
    labels <- paste(labels, " (", featureAnnotations, ")", sep="")
    axisSize <- 0.65 # for 10 features, 0.5 for 38 features    
  } else {
    axisSize <- 0.8
  }
  
  if(!is.null(axisScaling)) {
  	axisSize <- axisScaling  
  }
  
  if(numFeatures == 1) {
    axis(2, at=0, labels=labels, cex.axis=axisSize)  
  } else {
    axis(2, at=seq(0,1,1/(numFeatures-1)), labels=labels, cex.axis=axisSize)    
  }
  
  ### GI50
  par(mar=c(8,leftMarHeatmap,0,1))
  #mat <- matrix(sample(-10:10, numCellLines, replace=T), ncol=numCellLines)
  #labels <- rep("CELLLINE_A", 20)
  
  # NOTE: If the drug data is not transposed properly the drug activity bar may look odd.
  image(as.matrix(drugActValues), col=drugact.colors(12), axes=FALSE)
  axis(2, at=0.1, labels=responseLabel, cex.axis=0.75) 
  
  if(showCellLineLabels) {
    par(las=3)
    labels <- cellLineOrder
    axis(1, at=seq(0,1,1/(numCellLines-1)), labels=labels, cex.axis=0.60)
  }
  
  ### Weights
  weightColors <- rep(NA, length(weights[1:numFeatures]))
  weightColors[which(weights[1:numFeatures] > 0)] <- "green"
  weightColors[which(weights[1:numFeatures] < 0)] <- "red"
  
  #par(las=1, mar=c(2.75,0,4,2)) (10 features)
  #bottomMarWt <- 2.75
  #topMarWt <- 4

  #   par(las=1, mar=c(1.6,0,2.75,2)) (38 features)
  par(las=1, mar=c(bottomMarWt,0,topMarWt,2))  
  barplot2(weights[1:numFeatures], xlab="", ylab="", 
           horiz=TRUE, plot.grid=TRUE, space=0, main="Weights", 
           col=weightColors, names.arg=NULL)
  
  ##### SCALE BARS #####
  scaleBarSize = numCellLines
  
  ### Scale bar / Heatmap
  par(las=1, mar=c(4,2,1,2))
  #mat <- matrix(sample(-10:10, 20, replace=T), ncol=1)
  
  image(t(seq(0,1,1/(scaleBarSize-1))), col=heatmap.colors(scaleBarSize), axes=FALSE)
  axis(2, at=seq(0,1,1/2), labels=c(minHeatMapValues,midHeatMapValues,maxHeatMapValues), cex.axis=0.75)
  mtext(side=1, text="Values", line=1, cex=0.75)
  
  ### Scale bar / GI50
  par(las=1, mar=c(4,2,1,2))
  #mat <- matrix(sample(-10:10, 20, replace=T), ncol=1)
  
  image(t(seq(0,1,1/(scaleBarSize-1))), col=drugact.colors(scaleBarSize), axes=FALSE)
  axis(2, at=seq(0,1,1/2), labels=c(minCellLineValue,midCellLineValue,maxCellLineValue), cex.axis=0.75)
  mtext(side=1, text=responseLabel, line=1, cex=0.75)
  
  if(!is.null(pdfFilename)) {
    dev.off()
  }
  
  if(verbose) {
    enPlotValues <- list(drugActValues=drugActValues, heatMapValues=t(heatMapValues), 
                         cellLineOrder=cellLineOrder)
    return(enPlotValues)
  }
}

