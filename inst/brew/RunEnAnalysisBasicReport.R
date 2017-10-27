# DEBUG
Sys.setenv("NSCS_FILE"="")
Sys.setenv("RENDER_RMD"="TRUE")
Sys.setenv("NUM_CORES"="1")

Sys.setenv("TEMPLATE_FILE"=system.file("brew", "templateBasicBrew.Rmd", package="rcellminerElasticNet"))
Sys.setenv("CONFIG_FILE"=system.file("brew", "config.json", package="rcellminerElasticNet"))

Sys.setenv("SWATH_OUTPUT"="/Users/cannin/Downloads/delSO")

genScript <- system.file("brew", "genEnAnalysisBasicReport.R", package="rcellminerElasticNet")
source(genScript)

# BIOWULF (ALL)

gsplit --numeric-suffixes=1 -l 1000 allNscsIds.txt allNscsIds

Sys.setenv("NSCS_FILE"="/data/augustin/allNscsIds_1001_2000.txt")
Sys.setenv("RENDER_RMD"="FALSE")
Sys.setenv("NUM_CORES"="1")

Sys.setenv("TEMPLATE_FILE"=system.file("brew", "templateElasticNetBatchAnalysisBasicBrew", package="rcellminerElasticNet"))
Sys.setenv("CONFIG_FILE"=system.file("brew", "config.json", package="rcellminerElasticNet"))

Sys.setenv("SWATH_OUTPUT"="/data/augustin/swathOutput_020916_1001_2000")

genScript <- system.file("brew", "genEnAnalysisBasicReport.R", package="rcellminerElasticNet")
source(genScript)

# BIOWULF (FDA/CLIN/COMMON)

Sys.setenv("NSCS_FILE"="/data/augustin/fdaClinNscs.txt")
Sys.setenv("RENDER_RMD"="FALSE")
Sys.setenv("NUM_CORES"="1")

Sys.setenv("TEMPLATE_FILE"=system.file("brew", "templateElasticNetBatchAnalysisBasicBrew", package="rcellminerElasticNet"))
Sys.setenv("CONFIG_FILE"="/data/augustin/commonFeaEnConfig.json")

Sys.setenv("SWATH_OUTPUT"="/data/augustin/swathOutput_022416_common")
Sys.setenv("SWARM_OUTPUT"="/data/augustin/swathOutput_022416_common.swarm")

genScript <- system.file("brew", "genEnAnalysisBasicReport.R", package="rcellminerElasticNet")
source(genScript)

genSwarmScript <- system.file("brew", "genSwarmCmd.R", package="rcellminerElasticNet")
source(genSwarmScript)

workDir <- Sys.getenv("SWATH_OUTPUT")
swarmFile <- Sys.getenv("SWARM_OUTPUT")

files <- dir(".", pattern="Rmd")
tmp <- strsplit(files, "\\.Rmd")
ids <- unlist(tmp)

genSwarmCmd(ids, workDir, outputFile)

swarmCmd <- paste0("swarm -f ", swarmFile, " -g 4 --time=48:00:00 --module R --verbose 5")
#system(swarmCmd)
