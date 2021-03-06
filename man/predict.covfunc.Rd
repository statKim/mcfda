% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov.R
\name{predict.covfunc}
\alias{predict.covfunc}
\title{predict cov functions at new locations}
\usage{
\method{predict}{covfunc}(covobj, newt)
}
\arguments{
\item{covobj}{the object obtained by calling \code{covfunc}}

\item{newt}{a vector or a list of vectors of real numbers}
}
\value{
the estimated mean function evaluated at \code{newt}. It has the same format of \code{newt}
}
\description{
predict cov functions at new locations
}
