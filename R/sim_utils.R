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

#' Simulate Spatial Data with Directional Associations
#'
#' Generates simulated spatial point patterns with biologically plausible asymmetric
#' spatial associations between cell types. The function creates a hierarchical dataset
#' with images nested within patients and patients within groups.
#'
#' @param n_images Total number of images to simulate.
#' @param n_patients Number of patients.
#' @param n_groups Number of patient groups.
#' @param cell_types Character vector of cell type names. The last type is treated as the target.
#' @param target_type Character string specifying which cell type is the target (default: last type).
#' @param window_size Size of the square observation window (default: 1500).
#' @param n_cells_per_type Number of cells per type per image (default: 150).
#' @param n_basis_functions Number of RBF basis functions (default: 3).
#' @param max_dist Maximum distance for basis functions (default: 75).
#' @param basis_function_sigma Spread of basis functions (default: 15).
#' @param seed Random seed for reproducibility (default: 2025).
#' @param verbose logical
#'
#' @return A list with elements:
#' \describe{
#'   \item{data}{A tibble with columns: x, y, cell_type, image_id}
#'   \item{true_params}{A list containing the true simulation parameters for comparison}
#' }
#' @export
simulate_spatial_data <- function(
    n_images = 4,
    n_patients = 2,
    n_groups = 2,
    cell_types = c("tumor", "immune", "stroma"),
    target_type = NULL,
    window_size = 1500,
    sigma_beta_global = 0.5,
    sigma_beta_indiv = 0.1,
    sigma_beta_local = 0.1,
    scale_sigmas = 1,
    n_cells_per_type = 150,
    n_basis_functions = 3,
    max_dist = 75,
    basis_function_sigma = 15,
    seed = 2025,
    verbose = FALSE
) {
  set.seed(seed)
  
  # Default target type is the last one
  if (is.null(target_type)) {
    target_type <- cell_types[length(cell_types)]
  }
  
  # Validate inputs
  if (!target_type %in% cell_types) {
    stop("target_type must be one of the provided cell_types")
  }
  if (n_images %% n_patients != 0) {
    stop("n_images must be evenly divisible by n_patients")
  }
  if (n_patients %% n_groups != 0) {
    stop("n_patients must be evenly divisible by n_groups")
  }
  
  # Setup hierarchical structure
  images_per_patient <- n_images / n_patients
  patients_per_group <- n_patients / n_groups
  
  # Generate simulation parameters
  type_idx <- which(cell_types == target_type)
  num_types <- length(cell_types)
  
  # Create patient groupings
  indiv_to_group <- rep(1:n_groups, each = patients_per_group)
  sample_to_indiv <- rep(1:n_patients, each = images_per_patient)
  
  # Generate hierarchical coefficients
  params <- make_simulation_parameters(
    mean_alpha = log(n_cells_per_type / (window_size^2)),
    sigma_beta_global = sigma_beta_global,
    sigma_beta_indiv = sigma_beta_indiv,
    sigma_beta_local = sigma_beta_local,
    scale_sigmas = scale_sigmas,
    num_pt_groups = n_groups,
    num_types = num_types,
    num_combos = num_types - 1,
    num_pot = n_basis_functions,
    indiv_to_group = indiv_to_group,
    num_pts = n_patients,
    num_images = n_images,
    sample_to_indiv = sample_to_indiv,
    seed = seed
  )
  
  # Create RBF basis functions
  basis <- make_rbfs(
    n_basis_functions = n_basis_functions,
    max_dist = max_dist,
    basis_function_sigma = basis_function_sigma
  )
  
  # Observation window
  W <- spatstat.geom::owin(c(0, window_size), c(0, window_size))
  area <- window_size^2
  
  # Extract local coefficients (excluding intercept)
  betas_local <- params$betas_local[-1, ]
  
  # Simulate spatial patterns for each image
  sim_data <- lapply(1:n_images, function(i) {
    if(verbose) print(i)
    # Simulate background cell types
    if (num_types == 2) {
      pat <- spatstat.random::rpoispp(lambda = n_cells_per_type / area, win = W)
      spatstat.geom::marks(pat) <- factor(cell_types[1])
    } else {
      pat <- spatstat.random::rmpoispp(
        lambda = rep(n_cells_per_type / area, num_types - 1), 
        win = W
      )
      # Assign proper cell type names
      mark_counts <- as.vector(table(spatstat.geom::marks(pat)))
      spatstat.geom::marks(pat) <- factor(
        rep(cell_types[1:(num_types-1)], mark_counts),
        levels = cell_types[1:(num_types-1)]
      )
    }
    
    # Create intensity surfaces for target cell type
    dens_list <- lapply(1:(num_types - 1), function(j) {
      coeffs <- betas_local[((j - 1) * n_basis_functions + 1):(j * n_basis_functions), i]
      custom_kernel <- Vectorize(function(x, y) {
        d <- sqrt(x^2 + y^2)
        sum(sapply(seq_along(coeffs), function(k) coeffs[k] * basis[[k]](d)))
      })
      subs <- spatstat.geom::unmark(spatstat.geom::subset.ppp(pat, spatstat.geom::marks(pat) == cell_types[j]))
      smooth_density_fft(subs, custom_kernel, resolution = 128)
    })
    
    # Combine intensity surfaces
    int_dens <- Reduce("+", dens_list)
    lambda_integral <- sum(exp(int_dens$v)) * (int_dens$xstep * int_dens$ystep)
    beta0 <- log(n_cells_per_type / lambda_integral)
    
    # Simulate target cell type
    target_pat <- spatstat.random::rpoispp(lambda = exp(int_dens + beta0))
    spatstat.geom::marks(target_pat) <- factor(target_type, levels = target_type)
    
    # Combine all cell types
    final_pat <- spatstat.geom::superimpose(pat, target_pat)
    
    # Convert to data frame
    df <- as.data.frame(final_pat)
    df$image_id <- factor(paste0("img_", i))
    df$cell_type <- factor(df$marks, levels = cell_types)
    
    return(df[, c("x", "y", "cell_type", "image_id")])
  })
  
  # Combine all simulations
  result <- do.call(rbind, sim_data)
  
  # Return data with true parameters for comparison
  return(list(
    data = tibble::tibble(result),
    true_params = list(
      betas_global = if (n_groups > 0) params$betas_global else NULL,
      betas_indiv = params$betas_indiv,
      betas_local = params$betas_local,
      basis_functions = basis,
      n_basis_functions = n_basis_functions,
      max_dist = max_dist,
      basis_function_sigma = basis_function_sigma
    )
  ))
}