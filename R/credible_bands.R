#' Add Simultaneous Credible Bands to SIC Posterior Data
#'
#' Computes simultaneous credible bands that control the family-wise error rate
#' across all distances. This ensures the entire curve falls within the band with
#' (1-alpha)% probability, rather than just pointwise coverage.
#'
#' @param sic_data Data frame containing rvar columns (from posterior package)
#'   representing SIC posteriors. Must have a column identifying distance coordinates.
#' @param distance_col Name of the column containing distance values (default: "distance").
#' @param sic_cols Character vector of rvar column names to compute bands for.
#'   If NULL (default), processes all rvar columns except the distance column.
#' @param alpha Significance level (default 0.05 for 95% simultaneous bands).
#' @param keep_rvar Logical; if TRUE, keeps the original rvar column. Default FALSE.
#'
#' @return A data frame with the original data plus columns for each SIC:
#'   \item{<sic>_mean}{Posterior mean}
#'   \item{<sic>_lower}{Simultaneous lower bound}
#'   \item{<sic>_upper}{Simultaneous upper bound}
#'
#' @details
#' The algorithm:
#' 1. Extracts mean and SD at each distance for each rvar column
#' 2. Standardizes posterior deviations: (draw - mean) / SD
#' 3. Finds maximum absolute deviation across all distances for each draw
#' 4. Uses the (1-alpha) quantile of these maxima as the critical value
#' 5. Constructs bands: mean ± critical_value × SD
#'
#' Simultaneous bands are wider than pointwise bands but provide stronger guarantees:
#' the entire curve is contained with (1-alpha)% probability.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(posterior)
#'
#' # Compute SIC posteriors
#' sics <- compute_sic_posterior(fit, prep, level = "group")
#'
#' # Add simultaneous bands
#' sics_with_bands <- sics %>%
#'   add_simultaneous_bands(alpha = 0.1)
#'
#' # Plot
#' ggplot(sics_with_bands, aes(x = distance, y = sic_mean)) +
#'   geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper), alpha = 0.2) +
#'   geom_line() +
#'   facet_wrap(~source)
#' }
#'
#' @export
add_simultaneous_bands <- function(sic_data,
                                    distance_col = "distance",
                                    sic_cols = NULL,
                                    alpha = 0.05,
                                    keep_rvar = FALSE) {

  # Validate inputs
  if (!distance_col %in% names(sic_data)) {
    stop(sprintf("distance_col '%s' not found in sic_data", distance_col))
  }

  # Identify rvar columns to process
  if (is.null(sic_cols)) {
    sic_cols <- names(sic_data)[sapply(sic_data, posterior::is_rvar)]
    sic_cols <- setdiff(sic_cols, distance_col)
  }

  if (length(sic_cols) == 0) {
    stop("No rvar columns found in sic_data")
  }

  # Validate that specified columns exist and are rvars
  missing_cols <- setdiff(sic_cols, names(sic_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
  }

  non_rvar_cols <- sic_cols[!sapply(sic_data[sic_cols], posterior::is_rvar)]
  if (length(non_rvar_cols) > 0) {
    stop(sprintf("Columns are not rvar type: %s", paste(non_rvar_cols, collapse = ", ")))
  }

  # Convert to data frame for processing
  result <- as.data.frame(sic_data)

  # Step 1: Extract mean and SD at each distance
  for (col in sic_cols) {
    result[[paste0(col, "_mean")]] <- as.vector(posterior::E(result[[col]]))
    result[[paste0(col, "_sd")]] <- as.vector(posterior::sd(result[[col]]))
  }

  # Step 2: Standardize residuals
  # Extract just the rvar columns into a separate data frame
  lp_data <- result[sic_cols]

  # Standardize each column
  lp_std <- lapply(sic_cols, function(col) {
    rvar_col <- lp_data[[col]]
    # Compute mean and sd for this column using posterior package functions
    col_mean <- mean(rvar_col)
    col_sd <- posterior::sd(rvar_col)  # Use posterior::sd for rvar
    # Standardize
    posterior::as_rvar((rvar_col - col_mean) / col_sd)
  })
  names(lp_std) <- sic_cols

  # Step 3: Find maximum deviation across distances for each posterior draw
  max_dev <- sapply(lp_std, function(col) posterior::rvar_max(abs(col)))

  # Step 4: Find critical value for simultaneous coverage
  z_score_band <- sapply(max_dev, quantile, probs = 1 - alpha, na.rm = TRUE)

  # Step 5: Construct simultaneous lower and upper bands
  for (col in sic_cols) {
    result[[paste0(col, "_lower")]] <- result[[paste0(col, "_mean")]] -
      z_score_band[[col]] * result[[paste0(col, "_sd")]]
    result[[paste0(col, "_upper")]] <- result[[paste0(col, "_mean")]] +
      z_score_band[[col]] * result[[paste0(col, "_sd")]]

    # Remove rvar column if requested
    if (!keep_rvar) {
      result[[col]] <- NULL
    }

    # Remove SD column (internal use only)
    result[[paste0(col, "_sd")]] <- NULL
  }

  return(result)
}


#' Add Pointwise Credible Bands to SIC Posterior Data
#'
#' Computes pointwise credible bands using simple quantiles at each distance.
#' Unlike simultaneous bands, these only guarantee pointwise coverage - each
#' individual point has (1-alpha) probability of being in the band, but the
#' entire curve may not be contained.
#'
#' @param sic_data Data frame containing rvar columns (from posterior package)
#'   representing SIC posteriors. Must have a column identifying distance coordinates.
#' @param distance_col Name of the column containing distance values (default: "distance").
#' @param sic_cols Character vector of rvar column names to compute bands for.
#'   If NULL (default), processes all rvar columns except the distance column.
#' @param alpha Significance level (default 0.05 for 95% pointwise bands).
#' @param keep_rvar Logical; if TRUE, keeps the original rvar column. Default FALSE.
#'
#' @return A data frame with the original data plus columns for each SIC:
#'   \item{<sic>_mean}{Posterior mean}
#'   \item{<sic>_lower}{Pointwise lower bound}
#'   \item{<sic>_upper}{Pointwise upper bound}
#'
#' @details
#' The algorithm:
#' 1. Extracts mean at each distance for each rvar column
#' 2. Computes (alpha/2) and (1-alpha/2) quantiles at each distance independently
#' 3. Returns mean with quantile-based bounds
#'
#' Pointwise bands are narrower than simultaneous bands but provide weaker guarantees:
#' each individual distance has (1-alpha)% coverage, but multiple testing issues arise.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(posterior)
#'
#' # Compute SIC posteriors
#' sics <- compute_sic_posterior(fit, prep, level = "group")
#'
#' # Add pointwise bands
#' sics_with_bands <- sics %>%
#'   add_pointwise_bands(alpha = 0.05)
#'
#' # Compare with simultaneous bands
#' sics_simul <- add_simultaneous_bands(sics, alpha = 0.05)
#' }
#'
#' @export
add_pointwise_bands <- function(sic_data,
                                 distance_col = "distance",
                                 sic_cols = NULL,
                                 alpha = 0.05,
                                 keep_rvar = FALSE) {

  # Validate inputs
  if (!distance_col %in% names(sic_data)) {
    stop(sprintf("distance_col '%s' not found in sic_data", distance_col))
  }

  # Identify rvar columns to process
  if (is.null(sic_cols)) {
    sic_cols <- names(sic_data)[sapply(sic_data, posterior::is_rvar)]
    sic_cols <- setdiff(sic_cols, distance_col)
  }

  if (length(sic_cols) == 0) {
    stop("No rvar columns found in sic_data")
  }

  # Validate that specified columns exist and are rvars
  missing_cols <- setdiff(sic_cols, names(sic_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Columns not found: %s", paste(missing_cols, collapse = ", ")))
  }

  non_rvar_cols <- sic_cols[!sapply(sic_data[sic_cols], posterior::is_rvar)]
  if (length(non_rvar_cols) > 0) {
    stop(sprintf("Columns are not rvar type: %s", paste(non_rvar_cols, collapse = ", ")))
  }

  # Process each SIC column
  result <- sic_data

  for (col in sic_cols) {
    rvar_col <- result[[col]]

    # Compute mean and quantiles at each distance
    result <- result %>%
      dplyr::mutate(
        !!paste0(col, "_mean") := as.vector(posterior::E(rvar_col)),
        !!paste0(col, "_lower") := sapply(rvar_col, quantile, probs = alpha / 2),
        !!paste0(col, "_upper") := sapply(rvar_col, quantile, probs = 1 - alpha / 2)
      )

    # Remove rvar column if requested
    if (!keep_rvar) {
      result[[col]] <- NULL
    }
  }

  return(result)
}
