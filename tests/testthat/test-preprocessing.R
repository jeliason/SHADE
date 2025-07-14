test_that("prepare_spatial_model_data works with fake data", {
  # Generate fake coordinates and cell types
  set.seed(123)
  n_cells <- 100
  x <- runif(n_cells, 0, 100)
  y <- runif(n_cells, 0, 100)
  image_id <- factor(sample(paste0("img", 1:4), n_cells, replace = TRUE))
  cell_type <- sample(c("Tcell", "Bcell"), n_cells, replace = TRUE)
  
  # Create minimal patient metadata
  patient_metadata <- tibble::tibble(
    Spot = paste0("img", 1:4),
    Patient = c("P1", "P1", "P2", "P3"),
    Group = c("G1", "G1", "G2", "G2")
  )
  
  # Run the function (without writing files)
  result <- prepare_spatial_model_data(
    x = x,
    y = y,
    cell_type = cell_type,
    image_id = image_id,
    patient_metadata = patient_metadata,
    n_dummy = 10,
    type_idx = 1,
    path = NULL
  )
  
  # Check that output is a list with expected elements
  expect_type(result, "list")
  expect_true(all(c("stan_data", "sparse_matrix", "metadata") %in% names(result)))
  
  stan_data <- result$stan_data
  sparse_matrix <- result$sparse_matrix
  metadata <- result$metadata
  
  # Check basic structure
  expect_true(is.list(stan_data))
  expect_s4_class(sparse_matrix, 'RsparseMatrix')
  expect_true(is.list(metadata))
  
  # Check that data dimensions are consistent
  expect_equal(nrow(sparse_matrix), length(stan_data$is_cell))
  expect_equal(length(stan_data$w), stan_data$n_nz)
  expect_equal(length(stan_data$v), stan_data$n_nz)
  expect_equal(length(stan_data$u), nrow(sparse_matrix) + 1)
  
  # Check metadata content
  expect_true("coef_names" %in% names(metadata))
  expect_true("types" %in% names(metadata))
  expect_true(all(metadata$types %in% unique(cell_type)))
  
  # Check error when patient_metadata is missing a required column
  bad_metadata <- patient_metadata |> dplyr::select(-Group)
  expect_error(
    prepare_spatial_model_data(
      x = x,
      y = y,
      cell_type = cell_type,
      image_id = image_id,
      patient_metadata = bad_metadata
    ),
    "patient_metadata must be a data.frame with columns"
  )
  
  # Check error when input lengths do not match
  expect_error(
    prepare_spatial_model_data(
      x = x[1:50],
      y = y,
      cell_type = cell_type,
      image_id = image_id,
      patient_metadata = patient_metadata
    ),
    "length"
  )
})