% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov.pace.R
\name{cov.pace}
\alias{cov.pace}
\title{PACE approach to covariance estimation}
\usage{
cov.pace(
  Lt,
  Ly,
  bw = NULL,
  newt = NULL,
  mu = NULL,
  tuning = "GCV",
  weig = "SUBJ",
  kernel = "epanechnikov",
  delta = NULL,
  binData = "AUTO",
  numBins = NULL,
  nRegGrid = 51
)
}
\arguments{
\item{Lt}{a list (for irregular design) or a vector (for regular design)}

\item{Ly}{a list (for irregular design) for a matrix (for regular design). If \code{Ly} is a matrix, then \code{ncol(Ly)} must be equal to \code{length(Lt)}}

\item{bw}{bandwidth or bandwidth candidates}

\item{delta}{the snippet parameter, only used for irregular design}
}
\description{
PACE approach to covariance estimation
}
\keyword{internal}
