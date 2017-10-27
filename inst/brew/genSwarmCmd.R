#' Generate a swarm filePrefix
#'
#' @param ids IDs in the form DATA_X_NSC_Y
#' @param workDir the directory with the Rmds
#' @param outputFile a path for the swarm file
#'
#' @return nothing
#'
#' @notes
#' # Rmd IDs using: ls -1 *.Rmd | cut -d '.' -f 1
#' #workDir <- "/data/augustin/swathOutput_022416_common"
#' #outputFile <- "swathCmd_common.swarm"
#' #idsFile <- "commonFeaturesEn/rmdIdsCommon.txt"
#' #ids <- read.table(idsFile, header=FALSE, stringsAsFactors=FALSE)
#' #ids <- ids$V1
genSwarmCmd <- function(ids=NULL, workDir=Sys.getenv("SWATH_OUTPUT"), outputFile=Sys.getenv("SWARM_OUTPUT")) {
  if(is.null(ids)) {
    files <- dir(workDir, pattern="Rmd")
    tmp <- strsplit(files, "\\.Rmd")
    ids <- unlist(tmp)
  }

  cmds <- NULL
  for(id in ids) {
    #cmd <- paste0("cd ", workDir, " && test ! -e ", id, ".Rdata && Rscript -e \"library(knitr); knit2html('", id, ".Rmd')\"")
    #cmd <- paste0("cd ", workDir, " && test ! -e ", id, ".Rdata && mkdir ", id, " && cp ", id, ".Rmd ", id, "/ && cd ", id, " && Rscript -e \"library(knitr); knit2html('", id, ".Rmd')\" && cp *.{html,Rdata,txt} ../")
    #cmd <- paste0("cd ", workDir, " ; mkdir -p ", id, " ; cp ", id, ".Rmd ", id, "/ ; cd ", id, " ; Rscript -e \"library(knitr); knit2html('", id, ".Rmd')\"")

    #cmd <- paste0("echo \"START_TIME: `date`\" ; cd ", workDir, " ; mkdir -p ", id, " ; cp ", id, ".Rmd ", id, "/ ; cd ", id, " ; Rscript -e \"library(knitr); purl('", id, ".Rmd'); source('", id, ".R')\" ; cp *.{html,Rdata,txt} ../ ; echo \"END_TIME: `date`\"")
    cmd <- paste0("echo \"START_TIME: `date`\" ; cd ", workDir, " ; mkdir -p ", id, " ; cp ", id, ".Rmd ", id, "/ ; cd ", id, " ; Rscript -e \"library(knitr); knit2html('", id, ".Rmd', force_v1=TRUE)\" ; cp *.{html,Rdata,txt} ../ ; echo \"END_TIME: `date`\"")

    cmds <- c(cmds, cmd)
  }

  write.table(cmds, file=outputFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
