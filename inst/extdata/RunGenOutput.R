library(rmarkdown)

dataDir <- "~/Dropbox/TMP/swath_en_results/jun2016_biowulf/"
outputDir <- "~/Dropbox/TMP/swath_en_results/jun2016_biowulf/output/"
allDataDir <- "swathOutput_062316"
commonDataDir <- "swathOutput_062316_common"  

# HEATMAP
rmarkdown::render(file.path("vignettes", "makeEnResultsHeatmap.Rmd"), 
                  params = list(commonFlag = TRUE, partialFlag = TRUE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "heatmap_ct_pt.html"), 
                  output_dir = outputDir)

rmarkdown::render(file.path("vignettes", "makeEnResultsHeatmap.Rmd"), 
                  params = list(commonFlag = FALSE, partialFlag = FALSE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "heatmap_cf_pf.html"), 
                  output_dir = outputDir)

rmarkdown::render(file.path("vignettes", "makeEnResultsHeatmap.Rmd"), 
                  params = list(commonFlag = TRUE, partialFlag = FALSE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "heatmap_ct_pf.html"), 
                  output_dir = outputDir)

rmarkdown::render(file.path("vignettes", "makeEnResultsHeatmap.Rmd"), 
                  params = list(commonFlag = FALSE, partialFlag = TRUE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "heatmap_cf_pt.html"), 
                  output_dir = outputDir)

# PREDICTOR TABLE
rmarkdown::render(file.path("vignettes", "makeEnResultsTable.Rmd"), 
                  params = list(commonFlag = TRUE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "table_ct.html"), 
                  output_dir = outputDir)

rmarkdown::render(file.path("vignettes", "makeEnResultsTable.Rmd"), 
                  params = list(commonFlag = FALSE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "table_cf.html"), 
                  output_dir = outputDir)

# PREDICTOR ANNOTATION TABLE 
rmarkdown::render(file.path("vignettes", "makePredictorAnnotTable.Rmd"), 
                  params = list(commonFlag = TRUE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "annot_table_ct.html"), 
                  output_dir = outputDir)

rmarkdown::render(file.path("vignettes", "makePredictorAnnotTable.Rmd"), 
                  params = list(commonFlag = FALSE, dataDir = dataDir, 
                                outputDir = outputDir, allDataDir=allDataDir, commonDataDir=commonDataDir), 
                  output_file = file.path(outputDir, "annot_table_cf.html"), 
                  output_dir = outputDir)
