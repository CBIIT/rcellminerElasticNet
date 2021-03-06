% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genFeatureAnnotations.R
\name{genFeatureAnnotations}
\alias{genFeatureAnnotations}
\title{Generate feature annotations}
\usage{
genFeatureAnnotations(molDb, elNetResults, responseVec, verbose = FALSE)
}
\arguments{
\item{molDb}{molecular data used for elastic net results}

\item{elNetResults}{elastic net results object}

\item{responseVec}{response vector used in the elastic net results}

\item{verbose}{a boolean, whether to display debugging information}
}
\value{
a vector of strings with feature annotations in the form: cr=CORRELATION 
  cm=CUMULATIVE_CORRELATION fq=FEATURE_FREQUENCY
  \itemize{
    \item{cr} {individual correlation of feature}
    \item{cm} {cumulative correlation of features up to this feature}
    \item{fq} {frequency of feature in multiple training runs}
  }
}
\description{
Generate feature annotations
}
\concept{rcellminerElasticNet}
