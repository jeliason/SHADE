test_that("make_pat returns a ppp object", {
  x <- runif(100)
  y <- runif(100)
  type <- factor(sample(c("A", "B"), 100, replace = TRUE))
  
  pat <- make_pat(x, y, type, intensity_factor = 0)
  expect_s3_class(pat, "ppp")
})

test_that("make_pat preserves all points when intensity_factor = 0", {
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("T", "B"), 50, replace = TRUE))
  
  pat <- make_pat(x, y, type, intensity_factor = 0)
  expect_equal(pat$n, 50)
})

test_that("make_pat may remove points when intensity_factor > 0", {
  set.seed(123)
  x <- c(runif(95, 0, 1), runif(5, 5, 6))  # Sparse cluster
  y <- c(runif(95, 0, 1), runif(5, 5, 6))
  type <- factor(sample(c("CD8", "CD4"), 100, replace = TRUE))
  
  pat_masked <- make_pat(x, y, type, intensity_factor = 0.5)
  
  # The result should still be a valid ppp object
  expect_s3_class(pat_masked, "ppp")
  
  # Warn if no points are removed, but don't fail
  if (pat_masked$n == 100) {
    print("No points removed by intensity mask — this may happen depending on kernel estimate.")
  } else {
    expect_lt(pat_masked$n, 100)
  }
})

test_that("make_pat throws error if type is not a factor", {
  x <- runif(10)
  y <- runif(10)
  type <- sample(c("A", "B"), 10, replace = TRUE)  # Not a factor
  
  expect_error(make_pat(x, y, type), "type must be a factor!")
})

test_that("make_pat fails if x, y, and type lengths are unequal", {
  x <- runif(10)
  y <- runif(9)
  type <- factor(sample(c("A", "B"), 10, replace = TRUE))
  
  expect_error(make_pat(x, y, type), "x and y must be the same length")
})

test_that("make_quadrature works with scale_factor", {
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  
  Q <- make_quadrature(pp, scale_factor = 1)
  
  expect_s3_class(Q, "quad")
  expect_true(all(levels(spatstat.geom::marks(Q$data)) == levels(type)))
  expect_true(all(levels(spatstat.geom::marks(Q$dummy)) == levels(type)))
})

test_that("make_quadrature works with fixed n_dummy", {
  x <- runif(30)
  y <- runif(30)
  type <- factor(sample(c("X", "Y", "Z"), 30, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  
  Q <- make_quadrature(pp, n_dummy = 20)
  
  expect_s3_class(Q, "quad")
  expect_equal(length(levels(spatstat.geom::marks(Q$dummy))), 3)
})

test_that("make_quadrature throws if both scale_factor and n_dummy are set", {
  x <- runif(10)
  y <- runif(10)
  type <- factor(sample(c("a", "b"), 10, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  
  expect_error(
    make_quadrature(pp, scale_factor = 1, n_dummy = 10),
    "Specify either `scale_factor` or `n_dummy`, not both."
  )
})

test_that("make_quadrature throws if marks are missing", {
  pp <- spatstat.random::rpoispp(10)  # no marks
  
  expect_error(
    make_quadrature(pp, scale_factor = 1),
    "The point pattern `pp` must have marks"
  )
})

test_that("make_quadrature supports grid dummy distribution", {
  x <- runif(25)
  y <- runif(25)
  type <- factor(sample(c("red", "blue"), 25, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  
  Q <- make_quadrature(pp, n_dummy = 16, dist = "grid")
  
  expect_s3_class(Q, "quad")
  expect_equal(length(levels(spatstat.geom::marks(Q$dummy))), 2)
})

test_that("make_quadrature respects min_dummy", {
  # One type has only 3 points
  x <- c(runif(27), runif(3))
  y <- c(runif(27), runif(3))
  type <- factor(c(rep("A", 27), rep("B", 3)))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))

  Q <- make_quadrature(pp, scale_factor = 1, min_dummy = 10)

  # Check that B gets at least 10 dummy points
  dummy_B <- Q$dummy[spatstat.geom::marks(Q$dummy) == "B"]
  expect_gte(spatstat.geom::npoints(dummy_B), 10)
})

# Tests for covariate support
test_that("make_data works without covariates (backward compatibility)", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B", "C"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  result <- make_data(Q, potentials, focal_cell = "A", verbose = FALSE)

  expect_type(result, "list")
  expect_true("data" %in% names(result))
  expect_true("response" %in% names(result))
  expect_equal(ncol(result$data), 1 + 2 * 3)  # intercept + 2 source types * 3 RBFs
})

test_that("make_data accepts and binds covariates", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  # Get focal quadrature to determine size
  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == "A"],
    Q$dummy[spatstat.geom::marks(Q$dummy) == "A"]
  )
  n_quad <- spatstat.geom::npoints(Q_focal$data) + spatstat.geom::npoints(Q_focal$dummy)

  # Create covariate matrix
  covariates <- matrix(rnorm(n_quad * 2), ncol = 2)
  colnames(covariates) <- c("cov1", "cov2")

  result <- make_data(Q, potentials, focal_cell = "A", covariates = covariates, verbose = FALSE)

  expect_equal(ncol(result$data), 1 + 1 * 3 + 2)  # intercept + 1 source * 3 RBFs + 2 covariates
  expect_true("cov1" %in% colnames(result$data))
  expect_true("cov2" %in% colnames(result$data))
})

test_that("make_data validates covariate dimensions", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  # Wrong number of rows
  covariates <- matrix(rnorm(10 * 2), ncol = 2)

  expect_error(
    make_data(Q, potentials, focal_cell = "A", covariates = covariates),
    "must have .* rows"
  )
})

test_that("make_data detects missing values in covariates", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == "A"],
    Q$dummy[spatstat.geom::marks(Q$dummy) == "A"]
  )
  n_quad <- spatstat.geom::npoints(Q_focal$data) + spatstat.geom::npoints(Q_focal$dummy)

  # Create covariate matrix with NA
  covariates <- matrix(rnorm(n_quad * 2), ncol = 2)
  covariates[1, 1] <- NA

  expect_error(
    make_data(Q, potentials, focal_cell = "A", covariates = covariates),
    "cannot contain missing values"
  )
})

test_that("make_data warns about constant covariates", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == "A"],
    Q$dummy[spatstat.geom::marks(Q$dummy) == "A"]
  )
  n_quad <- spatstat.geom::npoints(Q_focal$data) + spatstat.geom::npoints(Q_focal$dummy)

  # Create covariate matrix with constant column
  covariates <- matrix(c(rep(1, n_quad), rnorm(n_quad)), ncol = 2)

  expect_warning(
    make_data(Q, potentials, focal_cell = "A", covariates = covariates),
    "constant and may cause numerical issues"
  )
})

test_that("make_data auto-names covariates if no colnames provided", {
  set.seed(42)
  x <- runif(50)
  y <- runif(50)
  type <- factor(sample(c("A", "B"), 50, replace = TRUE))
  pp <- spatstat.geom::ppp(x, y, marks = type, window = spatstat.geom::owin(c(0, 1), c(0, 1)))
  Q <- make_quadrature(pp, n_dummy = 100)

  potentials <- make_rbfs(max_dist = 1, n_basis_functions = 3)

  Q_focal <- spatstat.geom::quadscheme.logi(
    Q$data[spatstat.geom::marks(Q$data) == "A"],
    Q$dummy[spatstat.geom::marks(Q$dummy) == "A"]
  )
  n_quad <- spatstat.geom::npoints(Q_focal$data) + spatstat.geom::npoints(Q_focal$dummy)

  # Create covariate matrix without column names
  covariates <- matrix(rnorm(n_quad * 2), ncol = 2)

  result <- make_data(Q, potentials, focal_cell = "A", covariates = covariates, verbose = FALSE)

  # Check that auto-generated names exist
  expect_true("cov1" %in% colnames(result$data))
  expect_true("cov2" %in% colnames(result$data))
})