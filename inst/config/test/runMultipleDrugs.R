library(rmarkdown)

outputDir <- "/Users/rajapaksevn/Dropbox/TMP/swath_en_results/TEST/swathOutput_062116"


rmdFiles <- dir(outputDir, pattern = ".Rmd", full.names = TRUE)

for (rmdFile in rmdFiles){
  rmarkdown::render(input = rmdFile, output_dir = outputDir)
}

