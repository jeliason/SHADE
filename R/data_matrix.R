#' Create a marked point pattern from coordinates and types
#'
#' @param x Numeric vector of x coordinates
#' @param y Numeric vector of y coordinates
#' @param type Factor vector indicating types of points
#' @param intensity_factor Numeric multiplier for intensity thresholding.
#'        Set to 0 to disable intensity masking.
#'
#' @return A \code{ppp} object (spatstat point pattern) with optional masking based on intensity
#' @export
make_pat <- function(x, y, type, intensity_factor = 0) {
  # Validate inputs
  stopifnot("type must be a factor!" = is.factor(type))
  stopifnot("x and y must be the same length" = length(x) == length(y))
  stopifnot("x and type must be the same length" = length(x) == length(type))
  
  # Define initial observation window based on x and y ranges
  window_ss <- spatstat.geom::owin(
    xrange = range(x),
    yrange = range(y)
  )
  
  # Create initial point pattern with marks
  pat <- spatstat.geom::ppp(x, y, window = window_ss, marks = type)
  
  # If intensity masking is enabled
  if (intensity_factor > 0) {
    # Remove marks to compute intensity over spatial locations
    unmarked_pat <- spatstat.geom::unmark(pat)
    
    # Estimate intensity using Diggle's bandwidth
    intensity_img <- spatstat.explore::density.ppp(unmarked_pat, sigma = spatstat.explore::bw.diggle)
    
    # Compute intensity threshold and convert mask to logical image
    intensity_threshold <- intensity_factor * spatstat.geom::intensity.ppp(unmarked_pat)
    mask_df <- spatstat.geom::as.data.frame.im(intensity_img > intensity_threshold)
    
    # Create new window from mask and update pattern
    window_mask <- spatstat.geom::owin(mask = mask_df)
    pat <- spatstat.geom::ppp(x, y, marks = type, window = window_mask)
  }
  
  return(pat)
}

#' Construct a quadrature scheme for spatial logistic regression
#'
#' Generates a quadrature scheme from a multitype spatial point pattern by adding 
#' dummy points of each type and pairing them with observed data points. This is 
#' used in approximating the conditional intensity via logistic regression in spatial models.
#'
#' The number of dummy points per cell type can be controlled either by a fixed value 
#' (`n_dummy`) or as a multiple of the observed count via a scaling factor (`scale_factor`).
#' Only one of `n_dummy` or `scale_factor` should be supplied.
#'
#' @param pp A \code{ppp} object (spatstat point pattern). Must be a multitype pattern (i.e., have factor-valued marks).
#' @param scale_factor Optional numeric value specifying the multiplier applied to the number of observed points 
#' for determining the number of dummy points per cell type. Must not be used with \code{n_dummy}.
#' @param n_dummy Optional integer specifying a fixed number of dummy points per cell type. Overrides \code{scale_factor}.
#' @param min_dummy Minimum number of dummy points per type (used only when \code{scale_factor} is set and observed count is low). Default is 100.
#' @param dist Character string indicating how to generate dummy points. Options are \code{"random"} (default) for homogeneous Poisson process, 
#' or \code{"grid"} for grid-based dummy points.
#'
#' @return A \code{quadscheme} object (class "quad") from \pkg{spatstat.geom}, with appropriately marked data and dummy points.
#' @export
make_quadrature <- function(pp, scale_factor = NULL, n_dummy = NULL, min_dummy = 100, dist = "random") {
  
  # Ensure mutual exclusivity of scale_factor and n_dummy
  if (!is.null(scale_factor) && !is.null(n_dummy)) {
    stop("Specify either `scale_factor` or `n_dummy`, not both.")
  }
  
  # Check if the pattern is multitype
  if (is.null(spatstat.geom::marks(pp))) {
    stop("The point pattern `pp` must have marks (be multitype) to use scale_factor or n_dummy.")
  }
  
  # Get unique types
  types <- levels(spatstat.geom::marks(pp))
  
  # Determine number of dummy points per type
  if (!is.null(scale_factor)) {
    type_counts <- table(spatstat.geom::marks(pp))
    type_counts[type_counts < min_dummy] <- min_dummy
    n_dummy_per_type <- ceiling(scale_factor * type_counts)
  } else if (!is.null(n_dummy)) {
    n_dummy_per_type <- setNames(rep(n_dummy, length(types)), types)  # Fixed value per type
  } else {
    stop("Either `scale_factor` or `n_dummy` must be specified.")
  }
  
  # Generate dummy points
  dummy_points <- list()
  for (type in types) {
    nd <- ceiling(sqrt(n_dummy_per_type[type]))  # Ensure it's a square number
    if(dist == "random") {
      dummy_points[[type]] <- spatstat.random::rpoispp(lambda = n_dummy_per_type[type] / spatstat.geom::area.owin(spatstat.geom::as.owin(pp)), 
                                                       win = spatstat.geom::as.owin(pp))
    } else if (dist == "grid") {
      dummy_points[[type]] <- spatstat.geom::default.dummy(pp,nd = n_dummy)
    }
    
    
    # Assign marks to dummy points
    spatstat.geom::marks(dummy_points[[type]]) <- factor(rep(type, spatstat.geom::npoints(dummy_points[[type]])))           
  }
  
  # Combine dummy points into a single point pattern
  dummy <- Reduce(spatstat.geom::superimpose, dummy_points)
  
  # Create quadrature scheme
  Q <- spatstat.geom::quadscheme.logi(pp, dummy)
  
  # Ensure factor levels match pp
  levels(spatstat.geom::marks(Q$data)) <- levels(spatstat.geom::marks(pp))
  levels(spatstat.geom::marks(Q$dummy)) <- levels(spatstat.geom::marks(pp))
  
  return(Q)
}


#' Make radial basis functions
#'
#' @param max_dist maximum distance to consider
#' @param n_basis_functions number of basis functions
#' @param basis_function_sigma standard deviation for all basis functions
#'
#' @return list of named radial basis functions
#' @export
make_rbfs <- function(max_dist,
                      n_basis_functions=6,
                      basis_function_sigma=8) {
  gaussian_rbf <- function(x, mu, sigma) {
    exp(-(x - mu)^2 / (2 * sigma^2))
  }
  
  basis_function_centers <- seq(0, max_dist, length.out = n_basis_functions)  # Equally spaced centers
  
  rbfs <- lapply(basis_function_centers,\(mu) {
    function(x) {
      gaussian_rbf(x,mu,basis_function_sigma)
    }
  })
  
  names(rbfs) <- paste0("rbf",1:n_basis_functions)
  
  rbfs
}

#' Make regression dataframe from quadrature
#'
#' @param Q Logistic regression quadrature scheme, as generated by make_quadrature
#' @param potentials Named list of potential functions, each taking a distance as input
#' @param focal_cell Focal cell type (character) for which the regression is computed
#' @param verbose Logical flag for verbosity (default: FALSE)
#'
#' @return A list containing the design matrix (`data`) with an intercept and dispersion features,
#'         and the response vector (`response`) indicating data (1) vs. dummy (0) points.
#' @export
make_data <- function(Q,
                      potentials,
                      focal_cell = NULL,
                      verbose = FALSE) {
  
  # Ensure a focal cell is specified
  if (is.null(focal_cell)) {
    stop("focal_cell must be specified")
  }
  
  # Compute dispersion features using the provided potential functions
  dispersions <- make_dispersions(Q, focal_cell, potentials, verbose)
  
  # Subset quadrature scheme to the focal cell points using spatstat.core's quadscheme.logi
  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == focal_cell],
    Q$dummy[spatstat.geom::marks(Q$dummy) == focal_cell]
  )
  
  # Drop unused factor levels in marks for both data and dummy points
  spatstat.geom::marks(Q_focal$data)  <- droplevels(spatstat.geom::marks(Q_focal$data))
  spatstat.geom::marks(Q_focal$dummy) <- droplevels(spatstat.geom::marks(Q_focal$dummy))
  
  # Construct the response vector: 1 for data points, 0 for dummy points
  response <- as.numeric(spatstat.geom::is.data(Q_focal))
  
  # Build the design matrix by combining an intercept with the dispersion features
  design_matrix <- cbind(beta0 = 1, dispersions)
  
  return(list(data = design_matrix, response = response))
}

#' Make dispersions dataframe
#'
#' @param Q Logistic regression quadrature scheme, as generated by make_quadrature
#' @param focal_cell Focal cell type (character) whose interactions we're modeling
#' @param potentials Named list of potential functions, each taking a distance and returning a scalar
#' @param verbose Logical flag for verbosity
#'
#' @return Sparse matrix of dispersion features
#' @export
make_dispersions <- function(Q,
                             focal_cell,
                             potentials,
                             verbose = FALSE) {
  
  # Check input validity
  if (!is.factor(spatstat.geom::marks(Q$data))) {
    stop("The marks on the points in Q must be factors.")
  }
  
  # Determine all interaction types (excluding focal)
  m <- spatstat.geom::marks(Q)
  types <- setdiff(levels(m), focal_cell)
  types_grid <- tidyr::expand_grid(t1 = focal_cell, t2 = types)
  
  # Subset Q for focal cell
  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == focal_cell],
    Q$dummy[spatstat.geom::marks(Q$dummy) == focal_cell]
  )
  spatstat.geom::marks(Q_focal$data) <- droplevels(spatstat.geom::marks(Q_focal$data))
  spatstat.geom::marks(Q_focal$dummy) <- droplevels(spatstat.geom::marks(Q_focal$dummy))
  
  # Remaining (non-focal) data points
  pat_rest <- Q$data[spatstat.geom::marks(Q$data) != focal_cell]
  spatstat.geom::marks(pat_rest) <- droplevels(spatstat.geom::marks(pat_rest))
  
  # Compute pairwise distances
  print_vb("Calculating pairwise distances...", verbose)
  pd <- spatstat.geom::crossdist(spatstat.geom::union.quad(Q_focal), pat_rest)
  
  # Extract marks for grouping
  m_data <- spatstat.geom::marks(pat_rest)
  m_all <- spatstat.geom::marks(Q_focal)
  
  p <- progressr::progressor(along = potentials)
  
  data <- lapply(seq_along(potentials), function(i) {
    print_vb(paste0("Starting potential: ", i), verbose)
    
    potential <- potentials[[i]]
    pot_name <- names(potentials)[i]
    
    # Apply potential function to distance matrix
    pot <- matrix(potential(pd), nrow = nrow(pd))
    diag(pot) <- 0  # remove self-interaction
    
    # Precompute interaction matrices for each type pair
    pot_ix_jx <- lapply(seq_len(nrow(types_grid)), function(j) {
      t1 <- types_grid$t1[j]
      t2 <- types_grid$t2[j]
      ix <- which(m_all == t1)
      jx <- which(m_data == t2)
      as.matrix(pot[ix, jx])
    })
    print_vb(paste0("Calculated pot_ix_jx: ", i), verbose)
    
    # Compute dispersion vectors
    dispersions <- lapply(seq_len(nrow(types_grid)), function(j) {
      t1 <- types_grid$t1[j]
      ix <- which(m_all == t1)
      vec <- numeric(nrow(pot))
      vec[ix] <- rowSums(pot_ix_jx[[j]])
      Matrix::Matrix(vec, sparse = TRUE)
    }) %>% do.call(cbind, .)
    
    colnames(dispersions) <- paste0(pot_name, "_", types_grid$t1, "_", types_grid$t2)
    p()
    dispersions
  }) %>% do.call(cbind, .)
  
  # Order columns by cell type interaction (optional but consistent)
  col_names <- colnames(data)
  interaction_ids <- stringr::str_extract(col_names, "_[A-Za-z0-9+ ]+_[A-Za-z0-9+ ]+")
  source_types <- stringr::str_replace(interaction_ids,paste0("_",focal_cell,"_"),"")
  lvs <- levels(m)
  lvs <- lvs[which(lvs != focal_cell)]
  col_order <- order(match(source_types,lvs))
  data <- data[, col_order]
  
  return(data)
}
