#' Find unfinished Swarm commands
#' inputFile <- "/data/augustin/swathCmd_common.swarm"
#' ouputFile <- "/data/augustin/swathCmdUnfinished.swarm"
#' workDir <- "/data/augustin/swathOutput_022416_common"
#'
#' @export
findUnfinished <- function(inputFile, outputFile, workDir) {
  swarmFile <- readLines(inputFile)

  rmds <- dir(workDir, ".Rmd")
  rdatas <- dir(workDir, ".Rdata")

  tmpRmd <- gsub('.*(DATA.*)\\..*','\\1', rmds)
  tmpRdata <- gsub('.*(DATA.*)\\..*','\\1', rdatas)

  toMatch <- setdiff(tmpRmd, tmpRdata)

  matches <- unique(grep(paste(toMatch,collapse="|"), swarmFile, value=TRUE))

  writeLines(matches, ouputFile)
}
