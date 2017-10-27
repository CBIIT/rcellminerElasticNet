#rsync -av --exclude '.git' --exclude '.Rproj.user/' ~/default/workspace/nci60imsb augustin@helix.nih.gov:/data/augustin/rPackages/

#rsync -av --exclude '.git' --exclude '.Rproj.user/' ~/default/workspace/rcellminerElasticNet augustin@helix.nih.gov:/data/augustin/rPackages/

# BASIC SETUP
cd /Users/rajapaksevn/REPOS/

# INSTALL PACKAGES
Rscript -e 'library(devtools); install("rcellminer")'
Rscript -e 'library(devtools); install("rcellminerData")'
Rscript -e 'library(devtools); install("rcellminerElasticNet")'

# Get list of FDA/CT NSCs
Rscript -e 'library(rcellminerElasticNet); getFdaClinicalTrialNscs("/Users/rajapaksevn/REPOS/rcellminerElasticNet/inst/config/fdaClinNscs.txt")'

# SET ENVIRONMENT VARIABLES
## NOTE: USE FULL PATHS FOR ENVIRONMENT VARIABLES
export DATE=$(date +"%m%d%y")
export NSCS_FILE=/Users/rajapaksevn/REPOS/rcellminerElasticNet/inst/config/fdaClinNscs.txt
export RENDER_RMD=FALSE
export NUM_CORES=1
export TEMPLATE_FILE=/Users/rajapaksevn/REPOS/rcellminerElasticNet/inst/brew/templateElasticNetBatchAnalysisBasicBrew.Rmd
export CONFIG_FILE=/Users/rajapaksevn/REPOS/rcellminerElasticNet/inst/config/test/config.json
export SWATH_OUTPUT=/Users/rajapaksevn/Dropbox/TMP/swath_en_results/TEST/swathOutput_$DATE
#export SWARM_OUTPUT=/data/augustin/swathOutput_$DATE.swarm
#printenv | sort # View variables

# MAKE OUTPUT DIRECTORY
mkdir /Users/rajapaksevn/Dropbox/TMP/swath_en_results/TEST/swathOutput_$DATE

# GENERATE RMD
Rscript -e 'source("/Users/rajapaksevn/REPOS/rcellminerElasticNet/inst/brew/genEnAnalysisBasicReport.R")'

# GENERATE SWARM THAT RENDERS EACH RMD FILE
#Rscript -e 'source("rPackages/rcellminerElasticNet/inst/brew/genSwarmCmd.R"); genSwarmCmd()'

# RUN SWARM
#swarm -f $SWARM_OUTPUT -g 4 --time=48:00:00 --module R --verbose 5
