# install RCellMiner in R 3.2.*
source("https://bioconductor.org/biocLite.R")
biocLite("rcellminer")
biocLite("rcellminerData")

# Open Project

Open Rproj file in RStudio 

# Install Dependencies 
    if (!"devtools" %in% installed.packages()) install.packages("devtools")
    
    setRepositories(ind=1:6)
    library(devtools)

    install_bitbucket("cbio_mskcc/rcellminerelasticnet",
      auth_user="discoverUser",
      password="discoverUserPassword",
      build_vignettes=FALSE,
      dependencies=TRUE,
      args="--no-multiarch")
    
    # Custom pheatmap version  
    install_github("cannin/pheatmap", dependencies=TRUE)

# Run Vignette

Open vignette (Rmd files) from vignettes/ folder and use "Knit HTML" button to generate the HTML file.

# Download RData files

wget -r -nH -N -R index.html http://cmuser:eeX3aishoo@sanderlab.org/swathAnalysis/

# Convert Rmd to R
library(knitr); purl(file.path("vignettes", "checkEnData.Rmd"), output=file.path("vignettes", "checkEnData.R"))
library(knitr); purl(file.path("vignettes", "makeEnResultsHeatmap.Rmd"), output=file.path("vignettes", "makeEnResultsHeatmap.R"))
library(knitr); purl(file.path("vignettes", "makeEnResultsTable.Rmd"), output=file.path("vignettes", "makeEnResultsTable.R"))


# Command to obtain elastic net errors
grep -h  "## Error" *.html | sort -u 

## Errors due to empty elastic net feature set
<pre><code>## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
<pre><code>## Error in apply(cvFoldEnPredictorWts, MARGIN = 1, FUN = function(x) {: dim(X) must have a positive length
<pre><code>## Error in kable_markdown(x = structure(character(0), .Dim = c(0L, 0L), .Dimnames = list(: the table must have a header (column names)
<pre><code>## Error in sqrt(elNetResults$fullEnCvResults$cvPredRsqared): non-numeric argument to mathematical function
<pre><code>## Error: identical(names(reducedFeatureSetAvgWts), names(reducedFeatureSetNonzeroWtPct)) is not TRUE

## Errors from missing data
<pre><code>## Error in molDB[[molDBType]][paste0(molDBType, includeFeatures), ]: subscript out of bounds