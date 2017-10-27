library(parallel)
library(jsonlite)

# PARAMETERS ----
if(Sys.getenv("NSCS_FILE") == "") {
  #nscs <- c("609699", "757804", "119875", "747856", "761431", "718781")
  #nscs <- c("747856", "609699")
  nscs <- c("761431", "764134", "752", "34462", "3088", "330507", "409962", "756717", "616348", "63878")
} else {
  tmp <- read.table(Sys.getenv("NSCS_FILE"), header=TRUE, stringsAsFactors=FALSE)
  nscs <- tmp[, 1]
}

if(Sys.getenv("RENDER_RMD") == "") {
  Sys.setenv("RENDER_RMD"="TRUE")
}

if(Sys.getenv("SWATH_OUTPUT") == "") {
  #Sys.setenv("SWATH_OUTPUT"="/")
  Sys.setenv("SWATH_OUTPUT"="/Users/rajapaksevn/Dropbox/TMP/swath_en_results/test_070616/")
}

if(Sys.getenv("TEMPLATE_FILE") == "") {
  Sys.setenv("TEMPLATE_FILE"=system.file("brew", "templateElasticNetBatchAnalysisBasicBrew.Rmd", package="rcellminerElasticNet"))
}

if(Sys.getenv("NUM_CORES") == "") {
  numCores <- detectCores()
} else {
  numCores <- as.numeric(Sys.getenv("NUM_CORES"))
}

if(Sys.getenv("CONFIG_FILE") == "") {
  Sys.setenv("CONFIG_FILE"=system.file("config", "fdaClin", "config.json", package="rcellminerElasticNet"))
}

# Setup Cluster
cl <- makeCluster(numCores)

configFile <- Sys.getenv("CONFIG_FILE")
config <- jsonlite::fromJSON(configFile)

# GENERATE COMBINATIONS TO BE TESTED ----
# NOTE: Necessary to convert a vector into a string representation of the call
# to make the vector
tmp <- lapply(config$datasets, function(x) {
  return(deparse(do.call(call, as.list(c("c", x)))))
})

includeDatasets <- unlist(tmp)

reports <- expand.grid(nscs, includeDatasets, stringsAsFactors=FALSE)
colnames(reports) <- c("nsc", "includeDatasets")

# NOTE: Use parSapply to render over multiple nodes; use sapply to just generate Rmd (easier to debug)
#parSapply(cl, 1:nrow(reports), function(x, reports) {
sapply(1:nrow(reports), function(x, reports) {
  library(brew)
  library(knitr)
  library(rmarkdown)

  # Function to generate report and render output
  create.report <- function(reports, prepend="", workingDir=workingDir, outputDir=getwd(), renderRmd=TRUE) {
    setwd(workingDir)

    # Set parameters
    timeStamp <- format(Sys.time(), "%m%d%y")

    nsc <- reports$nsc
    includeDatasets <- reports$includeDatasets

    tmp <- gsub("\"", "", reports$includeDatasets)
    tmp <- gsub(" ", "", tmp)
    tmp <- gsub(",", "_", tmp)
    tmp <- gsub("^c", "", tmp)
    data <- gsub("[\\(\\)]", "", tmp)

    filePrefix <- paste0(prepend, paste(c("DATA", data, "NSC", nsc), collapse="_"))

    rmdFile <- file.path(outputDir, paste0(filePrefix, ".Rmd"))
    errorFile <- file.path(outputDir, paste0("invalid_", filePrefix, ".txt"))
    outputFile <- paste0(filePrefix, ".html")

    tmplFile <- Sys.getenv("TEMPLATE_FILE")

    #DEBUG
    cat("tmplFile: ", tmplFile, " rmdFile: ", rmdFile, "\n")

    brew(tmplFile, rmdFile)

    # Generate report
    if(renderRmd) {
      tryCatch({
        render(rmdFile, output_file=outputFile, output_dir=outputDir)
      }, error=function(e) {
        #write.table(1, file=errorFile)
        cat(paste0("ERROR: ", e), file=errorFile, append=TRUE, sep = "\n")
      })
    }
  }

  # SETUP FIXED PARAMETERS ----
  # Load Configuration
  configFile <- Sys.getenv("CONFIG_FILE")
  config <- jsonlite::fromJSON(configFile)

  tmp <- config$excludeCellLines
  excludeCellLines <- deparse(do.call(call, as.list(c("c", tmp))))

  corThreshold <- config$corThreshold
  predCorFeatureThreshold <- config$predCorFeatureThreshold
  minNumResponsiveLines <- config$minNumResponsiveLines
  responseThreshold <- config$responseThreshold
  minNumMutations <- config$minNumMutations
  maxPlotPredictors <- config$maxPlotPredictors
  alphaValMin <- config$alphaValMin
  alphaValMax <- config$alphaValMax
  verbose <- config$verbose
  useCommonFeaturesInMolDB <- config$useCommonFeaturesInMolDB
  includeFeaturesPath <- config$includeFeaturesPath
  useModelYIntercept <- config$useModelYIntercept

  workingDir <- Sys.getenv("SWATH_OUTPUT")
  outputDir <- Sys.getenv("SWATH_OUTPUT")

  # Set FALSE for just testing
  renderRmd <- as.logical(Sys.getenv("RENDER_RMD"))

  create.report(reports=reports[x, , drop=FALSE], prepend=NULL, workingDir=workingDir, outputDir=outputDir, renderRmd=renderRmd)
}, reports)

stopCluster(cl)
