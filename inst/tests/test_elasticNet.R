test_that("expected elasticNet() output is produced", {
  require(rcellminerElasticNet)
  require(rcellminer)
  require(impute)
  require(testthat)
  
  drugActData <- exprs(getAct(rcellminerData::drugData))
  drugActProfile <- setNames(as.numeric(drugActData["609699", ]), colnames(drugActData))
  
  set.seed(1)
  expData <- getAllFeatureData(rcellminerData::molData)[["exp"]][, names(drugActProfile)]
  expData <- invisible(impute.knn(expData)$data)
  
  set.seed(1)
  featureData <- expData[c("SLFN11", "LOC353194", "LDB1"), ]
  randFeatureSet <- sample(setdiff(rownames(expData), rownames(featureData)), size = 96)
  featureData <- rbind(featureData, expData[randFeatureSet, ])
  
  set.seed(1)
  enResults_1se <- elasticNet(featureMat = featureData, responseVec = drugActProfile, 
                          id = "609699", alphaVals = seq(0.2, 1, length = 10), 
                          useLambdaMin = TRUE,  verbose = FALSE, standardize = FALSE, 
                          useOneStdErrRule = TRUE)
  
  # testing -------------------------------------------------------------------
  expect_equal(enResults_1se$predictorWts,
               c(SLFN11 = 0.6253101), tolerance = 10^-6)
  expect_equal(enResults_1se$predictorSelectionFreq,
               c(SLFN11 = 1), tolerance = 10^-6)
  
  expect_equal(enResults_1se$predictorWts_preModelSelection,
               c(SLFN11 = 0.6253101, LDB1=0.1367193, 
                 CD40=-0.1203740), 
               tolerance = 10^-6)
  expect_equal(enResults_1se$predictorSelectionFreq_preModelSelection,
               c(SLFN11 = 1, LDB1=1, CD40=0.99), 
               tolerance = 10^-6)
  #----------------------------------------------------------------------------
})