# Example: Adding Spatial Covariates to SHADE
#
# This example demonstrates how to include spatial covariates in SHADE models.
# Covariates can capture location-specific effects like distance from tissue
# boundary, local cell density, tissue regions, etc.

library(SHADE)

# Simulate spatial data with 3 cell types across 2 images
set.seed(123)

# Image 1
n1 <- 200
x1 <- runif(n1, 0, 100)
y1 <- runif(n1, 0, 100)
type1 <- factor(sample(c("CD8", "CD4", "Tumor"), n1, replace = TRUE))

# Image 2
n2 <- 250
x2 <- runif(n2, 0, 100)
y2 <- runif(n2, 0, 100)
type2 <- factor(sample(c("CD8", "CD4", "Tumor"), n2, replace = TRUE))

# Combine data
x <- c(x1, x2)
y <- c(y1, y2)
cell_type <- factor(c(as.character(type1), as.character(type2)))
image_id <- factor(c(rep("Image_1", n1), rep("Image_2", n2)))

# Patient metadata
patient_metadata <- data.frame(
  Spot = c("Image_1", "Image_2"),
  Patient = c("P1", "P1"),
  Group = c("Control", "Control")
)

# === WITHOUT COVARIATES (standard SHADE) ===
cat("\n=== Running SHADE without covariates ===\n")

result_no_cov <- prepare_spatial_model_data(
  x = x,
  y = y,
  cell_type = cell_type,
  image_id = image_id,
  patient_metadata = patient_metadata,
  n_dummy = 500,
  type_idx = 1,  # Model CD4 as focal type
  n_basis_functions = 3,
  max_dist = 50
)

cat(sprintf("Design matrix dimensions: %d rows x %d columns\n",
            nrow(result_no_cov$sparse_matrix),
            ncol(result_no_cov$sparse_matrix)))
cat("Column names:", paste(result_no_cov$metadata$coef_names, collapse = ", "), "\n")
cat("Number of covariates:", result_no_cov$metadata$n_covariates, "\n")


# === WITH COVARIATES ===
cat("\n=== Running SHADE with covariates ===\n")

# Strategy: Create quadrature schemes first, then create covariates based on them,
# then pass both to prepare_spatial_model_data

# Step 1: Create point patterns (reuse from above)
pats <- list(
  Image_1 = make_pat(x1, y1, type1),
  Image_2 = make_pat(x2, y2, type2)
)

# Step 2: Create quadrature schemes
Qs <- lapply(pats, make_quadrature, n_dummy = 500)

# Step 3: Create covariate matrices based on the quadrature schemes
# Each covariate must have one row per quadrature point (data + dummy) for the FOCAL cell type
focal_type <- "CD4"

covariate_list <- lapply(names(Qs), function(img_name) {
  Q <- Qs[[img_name]]

  # Subset to focal cell type
  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == focal_type],
    Q$dummy[spatstat.geom::marks(Q$dummy) == focal_type]
  )

  # Drop unused levels
  spatstat.geom::marks(Q_focal$data)  <- droplevels(spatstat.geom::marks(Q_focal$data))
  spatstat.geom::marks(Q_focal$dummy) <- droplevels(spatstat.geom::marks(Q_focal$dummy))

  # Get the correct number of rows: data points + dummy points
  n_data <- spatstat.geom::npoints(Q_focal$data)
  n_dummy <- spatstat.geom::npoints(Q_focal$dummy)
  n_quad <- n_data + n_dummy

  # Get coordinates for data and dummy separately, then combine
  coords_data <- spatstat.geom::coords(Q_focal$data)
  coords_dummy <- spatstat.geom::coords(Q_focal$dummy)
  all_x <- c(coords_data$x, coords_dummy$x)
  all_y <- c(coords_data$y, coords_dummy$y)

  # Create example covariates:
  # 1. Distance from center (tissue region proxy)
  center_x <- mean(all_x)
  center_y <- mean(all_y)
  dist_from_center <- sqrt((all_x - center_x)^2 + (all_y - center_y)^2)

  # 2. X position (gradient effect)
  x_position <- all_x

  # 3. Random covariate (e.g., staining intensity, local density measure)
  random_effect <- rnorm(n_quad, mean = 10, sd = 2)

  # Combine into matrix
  cov_matrix <- cbind(
    dist_from_center = dist_from_center,
    x_position = x_position,
    intensity = random_effect
  )

  return(cov_matrix)
})

names(covariate_list) <- names(Qs)

# Step 4: Run prepare_spatial_model_data with both quadrature schemes and covariates
result_with_cov <- prepare_spatial_model_data(
  x = x,
  y = y,
  cell_type = cell_type,
  image_id = image_id,
  patient_metadata = patient_metadata,
  n_dummy = 500,  # Not used when quadrature_list is provided
  type_idx = 1,  # Model CD4 as focal type
  covariate_list = covariate_list,  # Pass covariates
  quadrature_list = Qs,  # Pass pre-made quadrature schemes
  n_basis_functions = 3,
  max_dist = 50
)

cat(sprintf("Design matrix dimensions: %d rows x %d columns\n",
            nrow(result_with_cov$sparse_matrix),
            ncol(result_with_cov$sparse_matrix)))
cat("Column names:", paste(result_with_cov$metadata$coef_names, collapse = ", "), "\n")
cat("Number of covariates:", result_with_cov$metadata$n_covariates, "\n")
cat("Covariate column indices:", paste(result_with_cov$metadata$covariate_cols, collapse = ", "), "\n")
cat("Covariate names:", paste(result_with_cov$metadata$covariate_names, collapse = ", "), "\n")

cat("\n=== Summary ===\n")
cat(sprintf("Without covariates: %d columns (1 intercept + %d dispersion features)\n",
            ncol(result_no_cov$sparse_matrix),
            ncol(result_no_cov$sparse_matrix) - 1))
cat(sprintf("With covariates: %d columns (1 intercept + %d dispersion features + %d covariates)\n",
            ncol(result_with_cov$sparse_matrix),
            ncol(result_with_cov$sparse_matrix) - result_with_cov$metadata$n_covariates - 1,
            result_with_cov$metadata$n_covariates))

cat("\nThe result can now be passed to run_SHADE_model() for fitting.\n")
cat("The Stan model will estimate hierarchical coefficients for all features,\n")
cat("including the spatial covariates.\n")
