#' Compute Spatial Interaction Curve Posteriors
#'
#' Extracts posterior distributions of Spatial Interaction Curves (SICs) from a
#' fitted SHADE model at the specified hierarchical level (image, patient, or group).
#' Returns posteriors as rvar objects for flexible downstream analysis.
#'
#' @param fit A CmdStan model fit object returned by `run_SHADE_model()`.
#' @param prep The preprocessing object returned by `prepare_spatial_model_data()`.
#' @param level Character string specifying hierarchical level: "image", "patient", or "group".
#' @param indices Optional integer vector specifying which images/patients/groups to extract.
#'   If NULL (default), extracts all available at the specified level.
#' @param distance_seq Numeric vector of distances at which to evaluate SICs.
#'   Default: seq(0, 100, by = 1).
#' @param sources Optional character vector of source cell types to extract.
#'   If NULL (default), extracts all source types.
#'
#' @return A tibble with columns:
#'   \item{distance}{Distance value}
#'   \item{source}{Source cell type name}
#'   \item{level_id}{ID at the specified level (image/patient/group number)}
#'   \item{level_name}{Name at the specified level (if available in metadata)}
#'   \item{sic}{rvar object containing posterior draws of the SIC}
#'
#' @details
#' The function extracts coefficients at the specified hierarchical level:
#' \itemize{
#'   \item \strong{image}: Uses `beta_local` - image-specific coefficients
#'   \item \strong{patient}: Uses `beta_indiv` - patient-specific coefficients
#'   \item \strong{group}: Uses `beta_global` - group-level coefficients
#' }
#'
#' SICs represent how target cell density varies with distance from source cells.
#' The returned rvar objects can be used with `posterior` package functions or
#' passed to `add_simultaneous_bands()` or `add_pointwise_bands()` for uncertainty
#' quantification.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(posterior)
#'
#' # Fit model
#' fit <- run_SHADE_model(prep$stan_data)
#'
#' # Extract group-level SICs
#' sics_group <- compute_sic_posterior(fit, prep, level = "group")
#'
#' # Extract patient-level SICs for specific patients
#' sics_patients <- compute_sic_posterior(
#'   fit, prep,
#'   level = "patient",
#'   indices = c(1, 3, 5)
#' )
#'
#' # Add credible bands and plot
#' sics_group %>%
#'   add_simultaneous_bands(alpha = 0.1) %>%
#'   ggplot(aes(x = distance, y = sic_mean, color = level_name)) +
#'     geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper, fill = level_name), alpha = 0.2) +
#'     geom_line() +
#'     facet_wrap(~source)
#' }
#'
#' @export
compute_sic_posterior <- function(fit,
                                   prep,
                                   level = c("group", "patient", "image"),
                                   indices = NULL,
                                   distance_seq = seq(0, 100, by = 1),
                                   sources = NULL) {

  # Validate inputs
  level <- match.arg(level)

  if (!inherits(fit, "CmdStanMCMC") && !inherits(fit, "CmdStanVB")) {
    stop("fit must be a CmdStanMCMC or CmdStanVB object from run_SHADE_model()")
  }

  if (!is.list(prep) || !all(c("stan_data", "metadata") %in% names(prep))) {
    stop("prep must be the preprocessing object from prepare_spatial_model_data()")
  }

  # Extract metadata
  metadata <- prep$metadata
  potentials <- metadata$potentials
  all_types <- metadata$types
  coef_names <- metadata$coef_names
  n_basis <- length(potentials)

  # Determine target and source types
  # Target is the focal type stored in metadata
  target_type <- metadata$focal_type
  if (is.null(target_type)) {
    stop("focal_type not found in metadata. Please re-run prepare_spatial_model_data() with updated package.")
  }
  source_types <- setdiff(all_types, target_type)

  # Filter sources if specified
  if (!is.null(sources)) {
    sources <- intersect(sources, source_types)
    if (length(sources) == 0) {
      stop("No valid source types specified")
    }
    source_types <- sources
  }

  # Extract posterior draws as rvar
  draws <- posterior::as_draws_rvars(fit$draws())

  # Create design matrix from potentials
  x_des <- lapply(potentials, function(pot) pot(distance_seq)) %>%
    do.call(cbind, .)

  # Extract SICs based on level
  if (level == "group") {
    result <- .extract_group_sics(draws, x_des, distance_seq, source_types, target_type,
                                   n_basis, coef_names, metadata, indices)
  } else if (level == "patient") {
    result <- .extract_patient_sics(draws, x_des, distance_seq, source_types, target_type,
                                     n_basis, coef_names, metadata, prep$stan_data, indices)
  } else if (level == "image") {
    result <- .extract_image_sics(draws, x_des, distance_seq, source_types, target_type,
                                   n_basis, coef_names, metadata, prep$stan_data, indices)
  }

  return(result)
}


# Internal function to extract group-level SICs
.extract_group_sics <- function(draws, x_des, distance_seq, source_types, target_type,
                                 n_basis, coef_names, metadata, indices) {

  if (!"beta_global" %in% names(draws)) {
    stop("Model does not have group-level effects (beta_global not found)")
  }

  beta_global <- draws$beta_global
  n_groups <- dim(beta_global)[2]

  # Validate and set indices
  if (is.null(indices)) {
    indices <- 1:n_groups
  } else {
    if (any(indices < 1 | indices > n_groups)) {
      stop(sprintf("indices must be between 1 and %d (number of groups)", n_groups))
    }
  }

  # Extract SICs for each group and source
  result_list <- list()

  for (group_i in indices) {
    for (source in source_types) {
      # Find coefficient indices for this target-source pair
      coef_pattern <- paste0("_", target_type, "_", source)
      coef_ix <- grep(coef_pattern, coef_names, fixed = TRUE)

      if (length(coef_ix) == 0) {
        warning(sprintf("No coefficients found for %s -> %s", source, target_type))
        next
      }

      # Should have n_basis coefficients per source
      if (length(coef_ix) != n_basis) {
        warning(sprintf("Expected %d coefficients for %s -> %s, found %d",
                        n_basis, source, target_type, length(coef_ix)))
      }

      # Extract coefficients for this group and source
      group_coeffs <- beta_global[coef_ix, group_i]

      # Compute linear predictor: X %*% beta
      lp <- x_des %*% group_coeffs

      # Create result tibble
      result_list[[length(result_list) + 1]] <- tibble::tibble(
        distance = distance_seq,
        source = source,
        level_id = group_i,
        level_name = paste0("Group ", group_i),
        sic = lp
      )
    }
  }

  dplyr::bind_rows(result_list)
}


# Internal function to extract patient-level SICs
.extract_patient_sics <- function(draws, x_des, distance_seq, source_types, target_type,
                                   n_basis, coef_names, metadata, stan_data, indices) {

  if (!"beta_indiv" %in% names(draws)) {
    stop("Model does not have patient-level effects (beta_indiv not found)")
  }

  beta_indiv <- draws$beta_indiv
  n_patients <- dim(beta_indiv)[2]

  # Validate and set indices
  if (is.null(indices)) {
    indices <- 1:n_patients
  } else {
    if (any(indices < 1 | indices > n_patients)) {
      stop(sprintf("indices must be between 1 and %d (number of patients)", n_patients))
    }
  }

  # Try to get patient names from metadata
  patient_names <- if (!is.null(metadata$patient_ids)) {
    metadata$patient_ids
  } else {
    paste0("Patient ", 1:n_patients)
  }

  # Extract SICs for each patient and source
  result_list <- list()

  for (patient_i in indices) {
    for (source in source_types) {
      # Find coefficient indices for this target-source pair
      coef_pattern <- paste0("_", target_type, "_", source)
      coef_ix <- grep(coef_pattern, coef_names, fixed = TRUE)

      if (length(coef_ix) == 0) {
        warning(sprintf("No coefficients found for %s -> %s", source, target_type))
        next
      }

      # Extract coefficients for this patient and source
      patient_coeffs <- beta_indiv[coef_ix, patient_i]

      # Compute linear predictor
      lp <- x_des %*% patient_coeffs

      # Create result tibble
      result_list[[length(result_list) + 1]] <- tibble::tibble(
        distance = distance_seq,
        source = source,
        level_id = patient_i,
        level_name = patient_names[patient_i],
        sic = lp
      )
    }
  }

  dplyr::bind_rows(result_list)
}


# Internal function to extract image-level SICs
.extract_image_sics <- function(draws, x_des, distance_seq, source_types, target_type,
                                 n_basis, coef_names, metadata, stan_data, indices) {

  if (!"beta_local" %in% names(draws)) {
    stop("Model does not have image-level effects (beta_local not found)")
  }

  beta_local <- draws$beta_local
  n_images <- length(beta_local)

  # Validate and set indices
  if (is.null(indices)) {
    indices <- 1:n_images
  } else {
    if (any(indices < 1 | indices > n_images)) {
      stop(sprintf("indices must be between 1 and %d (number of images)", n_images))
    }
  }

  # Get image names from metadata
  image_names <- if (!is.null(metadata$spots)) {
    as.character(metadata$spots)
  } else {
    paste0("Image ", 1:n_images)
  }

  # Extract SICs for each image and source
  result_list <- list()

  for (image_i in indices) {
    for (source in source_types) {
      # Find coefficient indices for this target-source pair
      coef_pattern <- paste0("_", target_type, "_", source)
      coef_ix <- grep(coef_pattern, coef_names, fixed = TRUE)

      if (length(coef_ix) == 0) {
        warning(sprintf("No coefficients found for %s -> %s", source, target_type))
        next
      }

      # Extract coefficients for this image and source
      # beta_local is an rvar of length n_images, each element is a vector of d_cells coefficients
      image_coeffs <- beta_local[[image_i]][coef_ix]

      # Compute linear predictor
      lp <- x_des %*% image_coeffs

      # Create result tibble
      result_list[[length(result_list) + 1]] <- tibble::tibble(
        distance = distance_seq,
        source = source,
        level_id = image_i,
        level_name = image_names[image_i],
        sic = lp
      )
    }
  }

  dplyr::bind_rows(result_list)
}


#' Extract Group-Level Spatial Interaction Curves
#'
#' Convenience wrapper around `compute_sic_posterior()` that extracts group-level
#' SICs and optionally adds credible bands in one call.
#'
#' @inheritParams compute_sic_posterior
#' @param bands Type of credible bands to add: "simultaneous", "pointwise", or "none".
#'   Default is "simultaneous".
#' @param alpha Significance level for credible bands (default: 0.05).
#' @param keep_rvar Logical; if TRUE and bands != "none", keeps the rvar column.
#'   Default FALSE.
#'
#' @return A tibble with distance, source, level information, and either:
#'   \itemize{
#'     \item If `bands = "none"`: rvar column containing posterior
#'     \item Otherwise: mean, lower, and upper bound columns
#'   }
#'
#' @examples
#' \dontrun{
#' # Simple extraction with simultaneous bands
#' sics <- extract_group_sics(fit, prep)
#'
#' # Custom distance range and pointwise bands
#' sics <- extract_group_sics(
#'   fit, prep,
#'   distance_seq = seq(30, 100, by = 2),
#'   bands = "pointwise",
#'   alpha = 0.1
#' )
#' }
#'
#' @export
extract_group_sics <- function(fit,
                                prep,
                                indices = NULL,
                                distance_seq = seq(0, 100, by = 1),
                                sources = NULL,
                                bands = c("simultaneous", "pointwise", "none"),
                                alpha = 0.05,
                                keep_rvar = FALSE) {

  bands <- match.arg(bands)

  # Compute posteriors
  sics <- compute_sic_posterior(
    fit = fit,
    prep = prep,
    level = "group",
    indices = indices,
    distance_seq = distance_seq,
    sources = sources
  )

  # Add bands if requested
  if (bands == "simultaneous") {
    sics <- add_simultaneous_bands(sics, distance_col = "distance",
                                    sic_cols = "sic", alpha = alpha,
                                    keep_rvar = keep_rvar)
  } else if (bands == "pointwise") {
    sics <- add_pointwise_bands(sics, distance_col = "distance",
                                 sic_cols = "sic", alpha = alpha,
                                 keep_rvar = keep_rvar)
  }

  return(sics)
}


#' Extract Patient-Level Spatial Interaction Curves
#'
#' Convenience wrapper around `compute_sic_posterior()` that extracts patient-level
#' SICs and optionally adds credible bands in one call.
#'
#' @inheritParams extract_group_sics
#'
#' @return A tibble with distance, source, level information, and either rvar or
#'   summarized columns depending on `bands` argument.
#'
#' @examples
#' \dontrun{
#' # Extract for all patients with simultaneous bands
#' sics <- extract_patient_sics(fit, prep)
#'
#' # Extract for specific patients only
#' sics <- extract_patient_sics(fit, prep, indices = c(1, 3, 5))
#' }
#'
#' @export
extract_patient_sics <- function(fit,
                                  prep,
                                  indices = NULL,
                                  distance_seq = seq(0, 100, by = 1),
                                  sources = NULL,
                                  bands = c("simultaneous", "pointwise", "none"),
                                  alpha = 0.05,
                                  keep_rvar = FALSE) {

  bands <- match.arg(bands)

  # Compute posteriors
  sics <- compute_sic_posterior(
    fit = fit,
    prep = prep,
    level = "patient",
    indices = indices,
    distance_seq = distance_seq,
    sources = sources
  )

  # Add bands if requested
  if (bands == "simultaneous") {
    sics <- add_simultaneous_bands(sics, distance_col = "distance",
                                    sic_cols = "sic", alpha = alpha,
                                    keep_rvar = keep_rvar)
  } else if (bands == "pointwise") {
    sics <- add_pointwise_bands(sics, distance_col = "distance",
                                 sic_cols = "sic", alpha = alpha,
                                 keep_rvar = keep_rvar)
  }

  return(sics)
}


#' Extract Image-Level Spatial Interaction Curves
#'
#' Convenience wrapper around `compute_sic_posterior()` that extracts image-level
#' SICs. Note that image-level estimates are typically noisy, so credible bands
#' default to "none".
#'
#' @inheritParams extract_group_sics
#'
#' @return A tibble with distance, source, level information, and either rvar or
#'   summarized columns depending on `bands` argument.
#'
#' @examples
#' \dontrun{
#' # Extract for all images (no bands by default due to noise)
#' sics <- extract_image_sics(fit, prep, bands = "none")
#'
#' # Extract for specific images
#' sics <- extract_image_sics(fit, prep, indices = c(1, 2, 3))
#' }
#'
#' @export
extract_image_sics <- function(fit,
                                prep,
                                indices = NULL,
                                distance_seq = seq(0, 100, by = 1),
                                sources = NULL,
                                bands = c("none", "pointwise", "simultaneous"),
                                alpha = 0.05,
                                keep_rvar = FALSE) {

  bands <- match.arg(bands)

  # Compute posteriors
  sics <- compute_sic_posterior(
    fit = fit,
    prep = prep,
    level = "image",
    indices = indices,
    distance_seq = distance_seq,
    sources = sources
  )

  # Add bands if requested
  if (bands == "simultaneous") {
    sics <- add_simultaneous_bands(sics, distance_col = "distance",
                                    sic_cols = "sic", alpha = alpha,
                                    keep_rvar = keep_rvar)
  } else if (bands == "pointwise") {
    sics <- add_pointwise_bands(sics, distance_col = "distance",
                                 sic_cols = "sic", alpha = alpha,
                                 keep_rvar = keep_rvar)
  }

  return(sics)
}
