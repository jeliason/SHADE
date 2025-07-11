% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{write_json_chunked}
\alias{write_json_chunked}
\title{Write Large JSON Object in Chunks}
\usage{
write_json_chunked(data, file_path, chunk_size = 1000, verbose = TRUE)
}
\arguments{
\item{data}{A named list containing elements to serialize. Each element must be a numeric scalar, numeric vector, or numeric matrix.}

\item{file_path}{Path to the output JSON file.}

\item{chunk_size}{Number of elements per write chunk when writing long vectors or matrix rows. Default is 1000.}

\item{verbose}{Logical; if TRUE (default), prints progress messages for each top-level element and chunk.}
}
\value{
Invisibly returns NULL. Writes the JSON file to disk.
}
\description{
Efficiently writes a list-based data object (used for Stan input, etc.) to disk in JSON format.
Supports numeric scalars, vectors, and matrices. Large vectors and matrices are written in chunks
to reduce memory usage and allow partial serialization of large arrays.
}
\details{
The function writes JSON line-by-line to avoid creating large in-memory representations of large vectors or matrices.
Each element in the input list becomes a top-level key in the resulting JSON object.
For vectors and matrices, only numeric types are supported.

Unsupported types (e.g., lists, data frames, arrays > 2D) will trigger an error.
}
\examples{
\dontrun{
stan_data <- list(
  scalar = 42,
  vector = rnorm(10000),
  matrix = matrix(runif(100), nrow = 10)
)
write_json_chunked(stan_data, "output.json")
}

}
