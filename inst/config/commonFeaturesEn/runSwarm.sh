# LOAD R ON NIH BIOWULF
#sinteractive
#module load R/3.2.3_gcc-4.9.1

# BASIC SETUP
cd /data/augustin

# DOWNLOAD PACKAGES
git clone https://github.com/cannin/rcellminerData.git rPackages/rcellminerData
git clone https://github.com/cannin/rcellminer.git rPackages/rcellminer
git clone https://bitbucket.org/cbio_mskcc/rcellminerelasticnet.git rPackages/rcellminerElasticNet

# INSTALL PACKAGES
Rscript -e 'library(devtools); install("rPackages/rcellminerData")'
Rscript -e 'library(devtools); install("rPackages/rcellminer")'
Rscript -e 'library(devtools); install("rPackages/rcellminerElasticNet")'

# Get list of FDA/CT NSCs
Rscript -e 'library(rcellminerElasticNet); getFdaClinicalTrialNscs("/data/augustin/rPackages/rcellminerElasticNet/inst/config/fdaClinNscs.txt")'

# SET ENVIRONMENT VARIABLES
## NOTE: USE FULL PATHS FOR ENVIRONMENT VARIABLES
export DATE=$(date +"%m%d%y")
export NSCS_FILE=/data/augustin/rPackages/rcellminerElasticNet/inst/config/fdaClinNscs.txt
export RENDER_RMD=FALSE
export NUM_CORES=1
export TEMPLATE_FILE=/data/augustin/rPackages/rcellminerElasticNet/inst/brew/templateElasticNetBatchAnalysisBasicBrew.Rmd
export CONFIG_FILE=/data/augustin/rPackages/rcellminerElasticNet/inst/config/commonFeaturesEn/config.json
export SWATH_OUTPUT=/data/augustin/swathOutput_${DATE}_common
export SWARM_OUTPUT=/data/augustin/swathOutput_${DATE}_common.swarm
#printenv | sort # View variables

# MAKE OUTPUT DIRECTORY
mkdir swathOutput_${DATE}_common

# GENERATE RMD
Rscript -e 'source("rPackages/rcellminerElasticNet/inst/brew/genEnAnalysisBasicReport.R")'

# GENERATE SWARM THAT RENDERS EACH RMD FILE
Rscript -e 'source("rPackages/rcellminerElasticNet/inst/brew/genSwarmCmd.R"); genSwarmCmd()'

# RUN SWARM
swarm -f $SWARM_OUTPUT -g 4 --time=48:00:00 --module R/3.2.3_gcc-4.9.1 --verbose 5
