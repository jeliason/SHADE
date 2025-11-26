#' Prepare Spatial Model Input Data
#'
#' Prepares spatial point pattern data and design matrices for hierarchical Bayesian modeling
#' of cell interactions in tissues. Assumes input data has already been filtered to relevant cell types.
#'
#' @param x Numeric vector of x-coordinates.
#' @param y Numeric vector of y-coordinates.
#' @param cell_type Factor or character vector of cell types (filtered).
#' @param image_id Character or factor vector of image/sample IDs.
#' @param patient_metadata A data frame with columns: 'Spot', 'Patient' and'Group'.
#' @param n_dummy Number of dummy points for quadrature (default: 1000).
#' @param type_idx Index of the cell type to model (default: 1).
#' @param covariate_list Optional named list of covariate matrices, one per image. Each matrix should have
#'        one row per quadrature point (data + dummy points for the focal cell type) and one column per covariate.
#'        List names must match the levels of `image_id`. If NULL (default), no covariates are included.
#' @param quadrature_list Optional named list of pre-made quadrature schemes (from make_quadrature), one per image.
#'        If provided, these will be used instead of generating new quadrature schemes internally.
#'        List names must match the levels of `image_id`. Useful when you need to create covariates that
#'        match specific quadrature points.
#' @param potentials Optional list of custom potential (basis) functions. Each function should take a distance
#'        and return a scalar value. If NULL (default), creates standard RBFs using make_rbfs() with the
#'        n_basis_functions, max_dist, and basis_function_sigma parameters. Useful for custom basis functions
#'        such as truncated RBFs or alternative functional forms.
#' @param path Optional directory path to write output files. If NULL, returns data in memory.
#' @param mean_alpha Mean of prior on intercept (default: -10).
#' @param scale_sigmas Scale for prior on interaction strengths (default: 5).
#' @param scale_sigma_betas Scale for prior on feature coefficients (default: seq(5, 1, length.out = num_pot)).
#' @param scale_sigma_alpha Scale for intercept prior (default: 5).
#' @param n_basis_functions Number of RBF basis functions (default: 3). Only used if potentials is NULL.
#' @param max_dist Maximum distance for RBF support (default: 75). Only used if potentials is NULL.
#' @param basis_function_sigma Spread of RBF functions (default: 15). Only used if potentials is NULL.
#'
#' @return A list with elements `stan_data`, `sparse_matrix`, and `metadata` if `path` is NULL.
#' Otherwise, writes outputs and returns invisibly.
#' @export
prepare_spatial_model_data <- function(
    x, y, cell_type, image_id,
    patient_metadata,
    n_dummy = 1000,
    type_idx = 1,
    covariate_list = NULL,
    quadrature_list = NULL,
    potentials = NULL,
    path = NULL,
    mean_alpha = -10,
    scale_sigmas = 5,
    scale_sigma_betas = NULL,
    scale_sigma_alpha = 5,
    n_basis_functions = 3,
    max_dist = 75,
    basis_function_sigma = 15
) {
  stopifnot(length(x) == length(y),
            length(x) == length(cell_type),
            length(x) == length(image_id))
  
  required_cols <- c("Spot", "Patient", "Group")
  if (!is.data.frame(patient_metadata) || !all(required_cols %in% colnames(patient_metadata))) {
    stop("patient_metadata must be a data.frame with columns: Spot, Patient, Group")
  }
  
  if(!is.factor(image_id)) {
    stop("image_id must be a factor!")
  }
  
  spots_pt <- patient_metadata[patient_metadata$Spot %in% image_id, , drop = FALSE]
  patient_ids <- unique(spots_pt$Patient)
  sample_to_indiv <- match(spots_pt$Patient, patient_ids)
  num_indiv <- length(patient_ids)
  
  if(all(is.na(spots_pt$Group))) {
    indiv_to_group <- rep(0, length(patient_ids))
    num_pt_groups <- 0
  } else {
    patient_to_group <- setNames(spots_pt$Group[!duplicated(spots_pt$Patient)], 
                                 spots_pt$Patient[!duplicated(spots_pt$Patient)])
    group_ids <- unique(spots_pt$Group[!is.na(spots_pt$Group)])
    indiv_to_group <- match(patient_to_group[patient_ids], group_ids)
    indiv_to_group[is.na(indiv_to_group)] <- 0
    num_pt_groups <- length(group_ids)
  }
  
  df_all <- tibble::tibble(X = x, Y = y, type = cell_type, Spot = image_id)
  dats <- split(df_all, df_all$Spot)
  
  unique_types <- levels(droplevels(factor(cell_type)))
  dats <- lapply(dats, function(df) {
    df$type <- factor(df$type,levels = unique_types)
    df
  })
  
  pats <- lapply(dats, function(df) {
    tryCatch({
      pat <- make_pat(df$X, df$Y, df$type)
      win <- spatstat.geom::owin(xrange = range(df$X), yrange = range(df$Y))
      spatstat.geom::Window(pat) <- win
      pat
    }, error = function(e) {
      stop("Error constructing point pattern for a sample: ", conditionMessage(e))
    })
  })

  # Create potentials if not provided
  if (is.null(potentials)) {
    potentials <- make_rbfs(
      n_basis_functions = n_basis_functions,
      max_dist = max_dist,
      basis_function_sigma = basis_function_sigma
    )
  } else {
    # Validate custom potentials
    if (!is.list(potentials)) {
      stop("potentials must be a list of functions")
    }
    if (!all(sapply(potentials, is.function))) {
      stop("All elements of potentials must be functions")
    }
    # Update n_basis_functions to match provided potentials
    n_basis_functions <- length(potentials)
  }

  if (is.null(scale_sigma_betas)) {
    scale_sigma_betas <- seq(5, 1, length.out = n_basis_functions)
  }

  # Use provided quadrature schemes or create new ones
  if (!is.null(quadrature_list)) {
    # Validate quadrature_list
    if (!is.list(quadrature_list)) {
      stop("quadrature_list must be a list of quadrature schemes")
    }

    image_levels <- levels(image_id)
    if (is.null(names(quadrature_list))) {
      stop("quadrature_list must be a named list with names matching levels of image_id")
    }

    if (!all(names(quadrature_list) %in% image_levels)) {
      missing <- setdiff(names(quadrature_list), image_levels)
      stop(sprintf("quadrature_list contains names not in image_id levels: %s",
                   paste(missing, collapse = ", ")))
    }

    if (!all(image_levels %in% names(quadrature_list))) {
      missing <- setdiff(image_levels, names(quadrature_list))
      stop(sprintf("quadrature_list is missing quadrature schemes for image_id levels: %s",
                   paste(missing, collapse = ", ")))
    }

    # Verify all are valid quadrature objects
    if (!all(sapply(quadrature_list, inherits, "quad"))) {
      stop("All elements of quadrature_list must be quadrature schemes (class 'quad')")
    }

    Qs <- quadrature_list
  } else {
    Qs <- lapply(pats, make_quadrature, n_dummy = n_dummy)
  }

  type <- unique_types[type_idx]

  # Validate covariate_list if provided
  if (!is.null(covariate_list)) {
    if (!is.list(covariate_list)) {
      stop("covariate_list must be a list of matrices")
    }

    # Check that list names match image_id levels
    image_levels <- levels(image_id)
    if (is.null(names(covariate_list))) {
      stop("covariate_list must be a named list with names matching levels of image_id")
    }

    if (!all(names(covariate_list) %in% image_levels)) {
      missing <- setdiff(names(covariate_list), image_levels)
      stop(sprintf("covariate_list contains names not in image_id levels: %s",
                   paste(missing, collapse = ", ")))
    }

    if (!all(image_levels %in% names(covariate_list))) {
      missing <- setdiff(image_levels, names(covariate_list))
      stop(sprintf("covariate_list is missing matrices for image_id levels: %s",
                   paste(missing, collapse = ", ")))
    }

    # Check that all covariate matrices have the same number of columns
    n_cols <- sapply(covariate_list, ncol)
    if (length(unique(n_cols)) > 1) {
      stop("All covariate matrices must have the same number of columns (covariates)")
    }

    # Check that column names are consistent across images
    col_names_list <- lapply(covariate_list, colnames)
    first_names <- col_names_list[[1]]
    if (!all(sapply(col_names_list, function(x) identical(x, first_names)))) {
      warning("Covariate matrices have inconsistent column names across images")
    }
  }

  offset <- lapply(Qs, function(Q) {
    log(spatstat.geom::intensity(Q$dummy)) |>
      tibble::enframe() |>
      dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
      dplyr::filter(name == type) |>
      dplyr::pull(value)
  }) |> unlist()

  # Pass covariates to make_data if provided
  data_lists <- lapply(names(Qs), function(img_name) {
    Q <- Qs[[img_name]]
    cov <- if (!is.null(covariate_list)) covariate_list[[img_name]] else NULL
    make_data(Q, potentials, type, covariates = cov, verbose = FALSE)
  })
  
  x_cells <- do.call(rbind, lapply(data_lists, \(d) d$data))
  is_cell <- unlist(lapply(data_lists, \(d) d$response))
  sample_id <- unlist(lapply(seq_along(data_lists), function(i) rep(i, length(data_lists[[i]]$response))))
  
  start_stop_idx <- lapply(seq_along(data_lists), function(i) {
    idx <- which(sample_id == i)
    c(min(idx), max(idx))
  }) %>% do.call(rbind,.)
  
  x_cells_sparse <- as(x_cells, "RsparseMatrix")
  u <- as.integer(x_cells_sparse@p + 1)
  w <- x_cells_sparse@x
  v <- x_cells_sparse@j + 1
  
  data_start_stop <- lapply(seq_len(nrow(start_stop_idx)), function(i) {
    row_start <- start_stop_idx[i, 1]
    row_stop <- start_stop_idx[i, 2]
    c(u[row_start], u[row_stop + 1] - 1)
  }) %>% do.call(rbind,.)

  # Detect single-image patients for adaptive hierarchy
  images_per_patient <- table(sample_to_indiv)
  is_single_image_patient <- as.integer(images_per_patient[sample_to_indiv] == 1)

  data_stan <- list(
    num_indiv = num_indiv,
    num_types = length(unique_types),
    num_pot = n_basis_functions,
    num_pt_groups = num_pt_groups,
    n_cells = nrow(x_cells),
    d_cells = ncol(x_cells),
    is_cell = is_cell,
    oset = -offset,
    n_samples = length(data_lists),
    sample_to_indiv = sample_to_indiv,
    indiv_to_group = indiv_to_group,
    sample_id = sample_id,
    y_start_stop = start_stop_idx,
    data_start_stop = data_start_stop,
    n_nz = length(w),
    w = w,
    v = v,
    u = u,
    is_single_image_patient = is_single_image_patient,
    mean_alpha = mean_alpha,
    scale_sigmas = scale_sigmas,
    scale_sigma_betas = scale_sigma_betas,
    scale_sigma_alpha = scale_sigma_alpha,
    grainsize = 1
  )
  
  # Identify which columns are covariates
  n_dispersion_cols <- n_basis_functions * (length(unique_types) - 1)
  n_total_cols <- ncol(x_cells)
  has_covariates <- !is.null(covariate_list)
  n_covariates <- if (has_covariates) n_total_cols - n_dispersion_cols - 1 else 0  # -1 for intercept

  covariate_cols <- if (has_covariates) {
    # Column indices for covariates (after intercept and dispersions)
    (n_dispersion_cols + 2):n_total_cols
  } else {
    integer(0)
  }

  metadata <- list(
    coef_names = colnames(x_cells),
    potentials = potentials,
    spots = unique(spots_pt$Spot),
    types = unique_types,
    focal_type = type,
    n_basis_functions = n_basis_functions,
    sample_to_indiv = sample_to_indiv,
    indiv_to_group = indiv_to_group,
    n_covariates = n_covariates,
    covariate_cols = covariate_cols,
    covariate_names = if (has_covariates) colnames(x_cells)[covariate_cols] else character(0)
  )
  
  if (!is.null(path)) {
    type_safe <- make.names(type)
    file_json <- file.path(path, paste0("data_stan_CRC_type_", type_safe, ".json"))
    file_dlist <- file.path(path, paste0("data_lists_CRC_type_", type_safe, ".rds"))
    file_metadata <- file.path(path, paste0("metadata_CRC_type_", type_safe, ".rds"))
    
    write_json_chunked(data_stan, file_json, chunk_size = 1e6)
    saveRDS(data_lists, file_dlist)
    saveRDS(metadata, file_metadata)
    invisible(NULL)
  } else {
    return(list(
      stan_data = data_stan,
      sparse_matrix = x_cells_sparse,
      metadata = metadata
    ))
  }
}