% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{prepare_spatial_model_data}
\alias{prepare_spatial_model_data}
\title{Prepare Spatial Model Input Data}
\usage{
prepare_spatial_model_data(
  x,
  y,
  cell_type,
  image_id,
  patient_metadata,
  n_dummy = 1000,
  type_idx = 1,
  path = NULL,
  mean_alpha = -10,
  scale_sigmas = 5,
  scale_sigma_betas = NULL,
  scale_sigma_alpha = 5,
  n_basis_functions = 3,
  max_dist = 75,
  basis_function_sigma = 15
)
}
\arguments{
\item{x}{Numeric vector of x-coordinates.}

\item{y}{Numeric vector of y-coordinates.}

\item{cell_type}{Factor or character vector of cell types (filtered).}

\item{image_id}{Character or factor vector of image/sample IDs.}

\item{patient_metadata}{A data frame with columns: 'Spot', 'Patient' and'Group'.}

\item{n_dummy}{Number of dummy points for quadrature (default: 1000).}

\item{type_idx}{Index of the cell type to model (default: 1).}

\item{path}{Optional directory path to write output files. If NULL, returns data in memory.}

\item{mean_alpha}{Mean of prior on intercept (default: -10).}

\item{scale_sigmas}{Scale for prior on interaction strengths (default: 5).}

\item{scale_sigma_betas}{Scale for prior on feature coefficients (default: seq(5, 1, length.out = num_pot)).}

\item{scale_sigma_alpha}{Scale for intercept prior (default: 5).}

\item{n_basis_functions}{Number of RBF basis functions (default: 3).}

\item{max_dist}{Maximum distance for RBF support (default: 75).}

\item{basis_function_sigma}{Spread of RBF functions (default: 15).}
}
\value{
A list with elements `stan_data`, `sparse_matrix`, and `metadata` if `path` is NULL.
Otherwise, writes outputs and returns invisibly.
}
\description{
Prepares spatial point pattern data and design matrices for hierarchical Bayesian modeling
of cell interactions in tissues. Assumes input data has already been filtered to relevant cell types.
}
