#--------------------------------------------------------------------------------------------------
# SET INPUTS, ETC.
#--------------------------------------------------------------------------------------------------
library(rcellminerElasticNet)
library(rcellminer)
library(nci60imsb)
library(glmnet)
library(impute)

nci60Swath <- getAllFeatureData(nci60imsb::molData)[["swathms_avg"]]
nci60Exp <- getAllFeatureData(rcellminerData::molData)[["exp"]]
nci60Act <- exprs(getAct(rcellminerData::drugData))


# Exclude cell line with no data.
nci60Swath <- nci60Swath[, -which(colnames(nci60Swath) == "ME:MDA_N")]
nci60Exp <- nci60Exp[, colnames(nci60Swath)]
nci60Exp <- impute.knn(nci60Exp)$data
nci60Act <- nci60Act[, colnames(nci60Swath)]

drugNsc <- "757027"
rVec <- na.exclude(nci60Act[drugNsc, , drop=TRUE])
#fMat <- nci60Exp[, names(rVec)]
fMat <- nci60Swath[, names(rVec)]

rVecZ <- setNames(as.numeric(scale(rVec)), names(rVec))
fMatZ <- t(scale(t(fMat)))

set.seed(1)
source("R/getCvFoldIds.R")
cvFoldIds <- getCvFoldIds(nObs = length(rVec), nFolds = 10, nRepeats = 1)
#--------------------------------------------------------------------------------------------------
# ELASTIC NET, CV: ORIGINAL DATA (NON-ZSCORE), ALL FEATURES (NO 1 SE RULE), GLMNET STANDARDIZE
#--------------------------------------------------------------------------------------------------
set.seed(1)
enResults <- elasticNet(featureMat=fMat, responseVec=rVec, standardize = TRUE,
                        nFoldsForParamSelection=10, useOneStdErrRule = FALSE)

set.seed(1)
cvEnResults <- cvElasticNet(cvFoldIds=cvFoldIds, keepEnResults=TRUE,
                            featureMat=fMat, responseVec=rVec, standardize = TRUE,
                            nFoldsForParamSelection=10, useOneStdErrRule = FALSE)

save(cvEnResults, file = "inst/extdata/cvEnResults.RData")

#--------------------------------------------------------------------------------------------------
# ELASTIC NET, CV: ORIGINAL DATA (NON-ZSCORE), ALL FEATURES (NO 1 SE RULE), NO GLMNET STANDARDIZE
#--------------------------------------------------------------------------------------------------
set.seed(1)
enResults_no_glmnet_std <- elasticNet(featureMat=fMat, responseVec=rVec, standardize = FALSE,
                        nFoldsForParamSelection=10, useOneStdErrRule = FALSE)

set.seed(1)
cvEnResults_no_glmnet_std <- cvElasticNet(cvFoldIds=cvFoldIds, keepEnResults=TRUE,
                            featureMat=fMat, responseVec=rVec, standardize = FALSE,
                            nFoldsForParamSelection=10, useOneStdErrRule = FALSE)

save(cvEnResults_no_glmnet_std, file = "inst/extdata/cvEnResults_no_glmnet_std.RData")

#--------------------------------------------------------------------------------------------------
# ELASTIC NET, CV: ORIGINAL DATA (NON-ZSCORE), TRIM FEATURES (1 SE RULE), NO GLMNET STANDARDIZE
#--------------------------------------------------------------------------------------------------
set.seed(1)
enResults_no_glmnet_std_1se <- elasticNet(featureMat=fMat, responseVec=rVec, standardize = FALSE,
                                      nFoldsForParamSelection=10, useOneStdErrRule = TRUE)

set.seed(1)
cvEnResults_no_glmnet_std_1se <- cvElasticNet(cvFoldIds=cvFoldIds, keepEnResults=TRUE,
                                          featureMat=fMat, responseVec=rVec, standardize = FALSE,
                                          nFoldsForParamSelection=10, useOneStdErrRule = TRUE)

save(cvEnResults_no_glmnet_std_1se, file = "inst/extdata/cvEnResults_no_glmnet_std_1se.RData")

#--------------------------------------------------------------------------------------------------
# ELASTIC NET, CV: ORIGINAL DATA (NON-ZSCORE), TRIM FEATURES (1 SE RULE), NO GLMNET STANDARDIZE
#--------------------------------------------------------------------------------------------------
# testing leave out out cross-validation.

set.seed(1)
cvFoldIds <- getCvFoldIds(nObs = length(rVec), nFolds = length(rVec), nRepeats = 1)
cvEnResults_no_glmnet_std_1se <- cvElasticNet(minFeatureResponseAbsCor = 0.5, nTrainingRuns = 20, # For testing. 
                                              cvFoldIds=cvFoldIds, keepEnResults=TRUE,
                                              featureMat=fMat, responseVec=rVec, standardize = FALSE,
                                              nFoldsForParamSelection=10, useOneStdErrRule = TRUE)

#--------------------------------------------------------------------------------------------------
# DEBUGGGING
#--------------------------------------------------------------------------------------------------

set.seed(1)
enResults <- elasticNet(featureMat=fMatZ, responseVec=rVec, standardize = FALSE,
                        nFoldsForParamSelection=10, useOneStdErrRule = TRUE)

#--------------------------------------------------------------------------------------------------