#' Generate hierarchical parameters for SHADE simulations
#'
#' Generates a hierarchical set of coefficients for simulating spatial interactions.
#' Parameters are generated for global (group), individual (patient), and image (local) levels.
#'
#' @param mean_alpha Mean of the global intercept.
#' @param sigma_beta_global Standard deviation of global coefficients.
#' @param sigma_beta_indiv Standard deviation of individual-level deviations.
#' @param sigma_beta_local Standard deviation of image-level deviations.
#' @param scale_sigmas Not currently used.
#' @param num_pt_groups Number of patient groups.
#' @param num_types Number of cell types (last one is predicted).
#' @param num_combos Number of pairwise combinations (usually num_types - 1).
#' @param num_pot Number of basis functions.
#' @param indiv_to_group Integer vector mapping each individual to a group.
#' @param num_pts Number of individuals.
#' @param num_images Number of images.
#' @param sample_to_indiv Integer vector mapping each image to an individual.
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{betas_local}{Matrix of image-level coefficients.}
#'   \item{betas_indiv}{Matrix of individual-level coefficients.}
#'   \item{betas_global}{(Optional) Matrix of group-level coefficients.}
#' }
#' @export
make_simulation_parameters <- function(
    mean_alpha,
    sigma_beta_global,
    sigma_beta_indiv,
    sigma_beta_local,
    scale_sigmas,
    num_pt_groups,
    num_types,
    num_combos,
    num_pot,
    indiv_to_group,
    num_pts,
    num_images,
    sample_to_indiv,
    seed
) {
  set.seed(seed)
  num_params <- 1 + num_combos * num_pot  # 1 intercept + interaction terms
  
  # ---- Generate group-level coefficients ----
  if (num_pt_groups > 0) {
    betas_global <- vapply(seq_len(num_pt_groups), function(group_id) {
      beta0 <- rnorm(1, mean = mean_alpha, sd = sigma_beta_global)
      betas <- unlist(lapply(seq_len(num_combos), function(j) {
        b <- rnorm(num_pot, mean = 0, sd = sigma_beta_global)
        b[order(-abs(b))]  # enforce acyclicity
      }))
      c(beta0, betas)
    }, numeric(num_params))
    
    # ---- Individual-level coefficients ----
    betas_indiv <- vapply(seq_len(num_pts), function(i) {
      group_id <- indiv_to_group[i]
      rnorm(num_params, mean = betas_global[, group_id], sd = sigma_beta_indiv)
    }, numeric(num_params))
  } else {
    # No group-level variation
    betas_indiv <- vapply(seq_len(num_pts), function(i) {
      beta0 <- rnorm(1, mean = mean_alpha, sd = sigma_beta_indiv)
      betas <- unlist(lapply(seq_len(num_combos), function(j) {
        b <- rnorm(num_pot, mean = 0, sd = sigma_beta_indiv)
        b[order(-abs(b))]
      }))
      c(beta0, betas)
    }, numeric(num_params))
  }
  
  # ---- Image-level coefficients ----
  betas_local <- vapply(seq_len(num_images), function(i) {
    indiv_id <- sample_to_indiv[i]
    rnorm(num_params, mean = betas_indiv[, indiv_id], sd = sigma_beta_local)
  }, numeric(num_params))
  
  # ---- Output ----
  out <- list(
    betas_local = betas_local,
    betas_indiv = betas_indiv
  )
  if (num_pt_groups > 0) {
    out$betas_global <- betas_global
  }
  
  return(out)
}

#' Perform 2D FFT-Based Density Smoothing (Simple Version)
#'
#' @param X A point pattern (`ppp` object) or density grid (`im` object)
#' @param kernel_func A function that computes kernel values dynamically
#' @param resolution Resolution of the density grid if `X` is a `ppp` object
#' @return A smoothed density estimate as an `im` object
#' @export
smooth_density_fft <- function(X, kernel_func, resolution = 128) {
  # Convert point pattern to a density grid if needed
  if (spatstat.geom::is.ppp(X)) {
    X <- spatstat.geom::pixellate(X, dimyx = resolution, padzero = TRUE)
  }
  
  # Ensure input is an `im` object
  if (!spatstat.geom::is.im(X)) stop("X must be a `ppp` or `im` object")
  
  # Extract matrix and dimensions
  Y <- X$v
  nr <- nrow(Y)
  nc <- ncol(Y)
  
  # Pad the image to prevent wrap-around effects
  Ypad <- matrix(0, nrow = 2 * nr, ncol = 2 * nc)
  Ypad[1:nr, 1:nc] <- Y
  
  # Generate kernel values dynamically
  xcol.ker <- X$xstep * c(0:(nc - 1), -(nc:1))
  yrow.ker <- X$ystep * c(0:(nr - 1), -(nr:1))
  
  Kern <- outer(yrow.ker, xcol.ker, kernel_func)  # Compute dynamically
  
  # Compute FFT of image and kernel
  fft_Y <- stats::fft(Ypad)
  fft_Kern <- stats::fft(Kern)
  
  # Multiply in Fourier space and take inverse FFT
  smooth_Y <- Re(stats::fft(fft_Y * fft_Kern, inverse = TRUE)) / (4 * nc * nr)
  
  # Extract valid region
  smooth_Y <- smooth_Y[1:nr, 1:nc]
  
  # Convert back to `im` object
  smoothed_image <- spatstat.geom::im(smooth_Y, xcol = X$xcol, yrow = X$yrow, unitname = spatstat.geom::unitname(X))
  
  return(smoothed_image)
}