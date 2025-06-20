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
#' @param path Optional directory path to write output files. If NULL, returns data in memory.
#' @param mean_alpha Mean of prior on intercept (default: -10).
#' @param scale_sigmas Scale for prior on interaction strengths (default: 5).
#' @param scale_sigma_betas Scale for prior on feature coefficients (default: seq(5, 1, length.out = num_pot)).
#' @param scale_sigma_alpha Scale for intercept prior (default: 5).
#' @param n_basis_functions Number of RBF basis functions (default: 3).
#' @param max_dist Maximum distance for RBF support (default: 75).
#' @param basis_function_sigma Spread of RBF functions (default: 15).
#'
#' @return A list with elements `stan_data`, `sparse_matrix`, and `metadata` if `path` is NULL.
#' Otherwise, writes outputs and returns invisibly.
#' @export
prepare_spatial_model_data <- function(
    x, y, cell_type, image_id,
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
  
  potentials <- make_rbfs(
    n_basis_functions = n_basis_functions,
    max_dist = max_dist,
    basis_function_sigma = basis_function_sigma
  )
  
  if (is.null(scale_sigma_betas)) {
    scale_sigma_betas <- seq(5, 1, length.out = n_basis_functions)
  }
  
  Qs <- lapply(pats, make_quadrature, n_dummy = n_dummy)
  
  type <- unique_types[type_idx]
  
  offset <- lapply(Qs, function(Q) {
    log(spatstat.geom::intensity(Q$dummy)) |>
      tibble::enframe() |>
      dplyr::right_join(tibble::tibble(name = spatstat.geom::marks(Q)), by = "name") |>
      dplyr::filter(name == type) |>
      dplyr::pull(value)
  }) |> unlist()
  
  data_lists <- lapply(Qs, function(Q) make_data(Q, potentials, type, verbose = FALSE))
  
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
    mean_alpha = mean_alpha,
    scale_sigmas = scale_sigmas,
    scale_sigma_betas = scale_sigma_betas,
    scale_sigma_alpha = scale_sigma_alpha,
    grainsize = 1
  )
  
  metadata <- list(
    coef_names = colnames(x_cells),
    potentials = potentials,
    spots = unique(spots_pt$Spot),
    types = unique_types,
    sample_to_indiv = sample_to_indiv,
    indiv_to_group = indiv_to_group
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