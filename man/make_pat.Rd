% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_matrix.R
\name{make_pat}
\alias{make_pat}
\title{Create a marked point pattern from coordinates and types}
\usage{
make_pat(x, y, type, intensity_factor = 0)
}
\arguments{
\item{x}{Numeric vector of x coordinates}

\item{y}{Numeric vector of y coordinates}

\item{type}{Factor vector indicating types of points}

\item{intensity_factor}{Numeric multiplier for intensity thresholding.
Set to 0 to disable intensity masking.}
}
\value{
A \code{ppp} object (spatstat point pattern) with optional masking based on intensity
}
\description{
Create a marked point pattern from coordinates and types
}
