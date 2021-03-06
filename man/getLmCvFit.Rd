% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLmCvFit.R
\name{getLmCvFit}
\alias{getLmCvFit}
\title{Estimate the prediction error of a linear model using repeated k-fold
cross validation.}
\usage{
getLmCvFit(X, y, nFolds = 10, nRepeats = 10, cvFoldIds = NULL)
}
\arguments{
\item{X}{A data matrix with observation data along the rows and predictor data
along the columns.}

\item{y}{A response vector.}

\item{nFolds}{The number of folds used in cross-validation.}

\item{nRepeats}{The number of cross-validation repeats (with different fold
assignments).}

\item{cvFoldIds}{An nObs x nReplicates matrix with integer fold identifiers 
along the columns.}
}
\value{
A list with the following elements
\itemize{
\item{cvPred} {A vector of cross-validation predicted values (averaged over cross validation
 repeats).}
\item{cvPredR} {The Pearson's correlation between the above vector of predicted response 
 values and the actual response values.}
\item{cvPredRsqared} {The square of the Pearson's correlation between the above vector of
 predicted response values and the actual response values.}
\item{cvMeanSqErr} {The mean squared error (averaged over the results from nFolds x nRepeats
 sets of cross-validation predictions).}
\item{cvSdMeanSqErr} {The standard deviation of the above (nFolds x nRepeats) cross validation
 mean squared error values.}
\item{cvPredMat} {Matrix with predicted cross-validation predictions}
\item{cvMeanSqErrMat} {Matrix with the cross-validation mean squared errors}
}
}
\description{
Estimate the prediction error of a linear model using repeated k-fold
cross validation.
}
\concept{rcellminerElasticNet}
