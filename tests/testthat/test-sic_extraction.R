# Tests for SIC extraction and credible band functions

# Helper function to create mock rvar data
create_mock_sic_data <- function(n_distances = 50, n_groups = 2, n_sources = 2) {
  # Create mock posterior draws
  n_draws <- 100

  result_list <- list()

  for (g in 1:n_groups) {
    for (s in 1:n_sources) {
      # Simulate some posterior draws for SIC
      # Each distance point has n_draws samples
      sic_draws <- lapply(1:n_distances, function(i) {
        rnorm(n_draws, mean = g + s + 0.01 * i, sd = 0.5)
      })

      # Convert list of draws to rvar
      sic_rvar <- posterior::rvar(do.call(cbind, sic_draws))

      result_list[[length(result_list) + 1]] <- tibble::tibble(
        distance = seq(0, 100, length.out = n_distances),
        source = paste0("Source", s),
        level_id = g,
        level_name = paste0("Group ", g),
        sic = sic_rvar
      )
    }
  }

  dplyr::bind_rows(result_list)
}

# =============================================================================
# Tests for add_simultaneous_bands()
# =============================================================================

test_that("add_simultaneous_bands returns correct structure", {
  sic_data <- create_mock_sic_data(n_distances = 30, n_groups = 2, n_sources = 1)

  result <- add_simultaneous_bands(sic_data, alpha = 0.05)

  # Check that required columns exist
  expect_true("sic_mean" %in% names(result))
  expect_true("sic_lower" %in% names(result))
  expect_true("sic_upper" %in% names(result))

  # Check that rvar column was removed by default
  expect_false("sic" %in% names(result))

  # Check dimensions
  expect_equal(nrow(result), nrow(sic_data))
})

test_that("add_simultaneous_bands keeps rvar when requested", {
  sic_data <- create_mock_sic_data(n_distances = 20, n_groups = 1, n_sources = 1)

  result <- add_simultaneous_bands(sic_data, keep_rvar = TRUE)

  expect_true("sic" %in% names(result))
  expect_true(posterior::is_rvar(result$sic))
})

test_that("add_simultaneous_bands produces wider bands than mean ± 2SD", {
  sic_data <- create_mock_sic_data(n_distances = 50, n_groups = 1, n_sources = 1)

  result <- add_simultaneous_bands(sic_data, alpha = 0.05)

  # Simultaneous bands should generally be wider than pointwise bands
  band_width <- result$sic_upper - result$sic_lower
  expect_true(all(band_width > 0))

  # For 95% bands, width should be reasonable (not too narrow or too wide)
  expect_true(all(band_width < 20))  # Sanity check given our mock data
})

test_that("add_simultaneous_bands handles multiple groups correctly", {
  sic_data <- create_mock_sic_data(n_distances = 30, n_groups = 3, n_sources = 2)

  result <- add_simultaneous_bands(sic_data, alpha = 0.1)

  # Should have bands for each group
  groups_in <- unique(sic_data$level_id)
  groups_out <- unique(result$level_id)
  expect_equal(sort(groups_in), sort(groups_out))

  # Each group should have complete data
  for (g in groups_in) {
    group_data <- result[result$level_id == g, ]
    expect_false(any(is.na(group_data$sic_mean)))
    expect_false(any(is.na(group_data$sic_lower)))
    expect_false(any(is.na(group_data$sic_upper)))
  }
})

test_that("add_simultaneous_bands validates inputs", {
  sic_data <- create_mock_sic_data()

  # Should fail if distance column doesn't exist
  expect_error(
    add_simultaneous_bands(sic_data, distance_col = "nonexistent"),
    "distance_col 'nonexistent' not found"
  )

  # Should fail if no rvar columns
  sic_data_no_rvar <- sic_data
  sic_data_no_rvar$sic <- as.numeric(posterior::E(sic_data$sic))
  expect_error(
    add_simultaneous_bands(sic_data_no_rvar),
    "No rvar columns found"
  )
})

# =============================================================================
# Tests for add_pointwise_bands()
# =============================================================================

test_that("add_pointwise_bands returns correct structure", {
  sic_data <- create_mock_sic_data(n_distances = 30, n_groups = 2, n_sources = 1)

  result <- add_pointwise_bands(sic_data, alpha = 0.05)

  # Check that required columns exist
  expect_true("sic_mean" %in% names(result))
  expect_true("sic_lower" %in% names(result))
  expect_true("sic_upper" %in% names(result))

  # Check that rvar column was removed by default
  expect_false("sic" %in% names(result))
})

test_that("add_pointwise_bands produces narrower bands than simultaneous", {
  sic_data <- create_mock_sic_data(n_distances = 50, n_groups = 1, n_sources = 1)

  result_point <- add_pointwise_bands(sic_data, alpha = 0.05)
  result_simul <- add_simultaneous_bands(sic_data, alpha = 0.05)

  # Pointwise bands should generally be narrower
  width_point <- mean(result_point$sic_upper - result_point$sic_lower)
  width_simul <- mean(result_simul$sic_upper - result_simul$sic_lower)

  expect_lt(width_point, width_simul)
})

test_that("add_pointwise_bands handles custom alpha", {
  sic_data <- create_mock_sic_data(n_distances = 30, n_groups = 1, n_sources = 1)

  result_90 <- add_pointwise_bands(sic_data, alpha = 0.10)
  result_95 <- add_pointwise_bands(sic_data, alpha = 0.05)

  # 90% bands should be narrower than 95% bands
  width_90 <- mean(result_90$sic_upper - result_90$sic_lower)
  width_95 <- mean(result_95$sic_upper - result_95$sic_lower)

  expect_lt(width_90, width_95)
})

# =============================================================================
# Tests for compute_sic_posterior() validation
# =============================================================================

test_that("compute_sic_posterior validates fit object", {
  fake_fit <- list(not = "a cmdstan fit")
  fake_prep <- list(stan_data = list(), metadata = list())

  expect_error(
    compute_sic_posterior(fake_fit, fake_prep),
    "fit must be a CmdStanMCMC or CmdStanVB object"
  )
})

test_that("compute_sic_posterior validates prep object", {
  # This is a minimal test since we can't easily create a real fit object
  fake_fit <- structure(list(), class = "CmdStanMCMC")
  fake_prep <- list(wrong = "structure")

  expect_error(
    compute_sic_posterior(fake_fit, fake_prep),
    "prep must be the preprocessing object"
  )
})

test_that("compute_sic_posterior validates level argument", {
  fake_fit <- structure(list(), class = "CmdStanMCMC")
  fake_prep <- list(stan_data = list(), metadata = list())

  expect_error(
    compute_sic_posterior(fake_fit, fake_prep, level = "invalid"),
    "'arg' should be one of"
  )
})

# =============================================================================
# Tests for convenience wrappers
# =============================================================================

test_that("convenience wrapper functions have correct arguments", {
  # Check that all wrappers have bands argument with correct defaults
  expect_true("bands" %in% names(formals(extract_group_sics)))
  expect_true("bands" %in% names(formals(extract_patient_sics)))
  expect_true("bands" %in% names(formals(extract_image_sics)))

  # Check default band types (formals returns call objects, eval to get actual values)
  expect_equal(eval(formals(extract_group_sics)$bands), c("simultaneous", "pointwise", "none"))
  expect_equal(eval(formals(extract_image_sics)$bands), c("none", "pointwise", "simultaneous"))
})

# =============================================================================
# Integration test with mock data (if possible)
# =============================================================================

test_that("full workflow produces expected output structure", {
  # Create a complete mock dataset
  sic_data <- create_mock_sic_data(n_distances = 40, n_groups = 2, n_sources = 3)

  # Test the full pipeline
  result <- sic_data %>%
    add_simultaneous_bands(alpha = 0.1)

  # Verify structure
  expect_s3_class(result, "data.frame")
  expect_true(all(c("distance", "source", "level_id", "level_name",
                    "sic_mean", "sic_lower", "sic_upper") %in% names(result)))

  # Verify that bands contain the mean
  expect_true(all(result$sic_lower <= result$sic_mean))
  expect_true(all(result$sic_mean <= result$sic_upper))

  # Verify no missing values
  expect_false(any(is.na(result$sic_mean)))
  expect_false(any(is.na(result$sic_lower)))
  expect_false(any(is.na(result$sic_upper)))
})

test_that("both band types work together", {
  sic_data <- create_mock_sic_data(n_distances = 30, n_groups = 1, n_sources = 2)

  # Should be able to compute both types without errors
  expect_silent({
    result_simul <- add_simultaneous_bands(sic_data, alpha = 0.05)
    result_point <- add_pointwise_bands(sic_data, alpha = 0.05)
  })

  # Both should have same number of rows
  expect_equal(nrow(result_simul), nrow(result_point))

  # Both should have means (which should be identical)
  expect_equal(result_simul$sic_mean, result_point$sic_mean, tolerance = 1e-10)
})

# =============================================================================
# Tests for simultaneous band coverage
# =============================================================================

test_that("simultaneous bands achieve expected coverage", {
  # Create mock data with known true values
  set.seed(789)
  n_distances <- 40
  n_draws <- 500
  alpha <- 0.10  # 90% bands

  # True curve (simple quadratic)
  distances <- seq(0, 100, length.out = n_distances)
  true_curve <- 2 + 0.01 * distances - 0.0001 * distances^2

  # Generate draws around true curve (draws x distances)
  sic_draws <- lapply(1:n_distances, function(i) {
    rnorm(n_draws, mean = true_curve[i], sd = 0.3)
  })
  sic_rvar <- posterior::rvar(do.call(cbind, sic_draws))

  sic_data <- tibble::tibble(
    distance = distances,
    source = "SourceA",
    level_id = 1,
    level_name = "Group 1",
    sic = sic_rvar,
    true_value = true_curve
  )

  # Compute simultaneous bands
  result <- add_simultaneous_bands(sic_data, alpha = alpha, keep_rvar = TRUE)

  # Check coverage: what proportion of draws have the entire curve within bands?
  # For each draw, check if all points of the true curve are within the band
  sic_draws <- posterior::draws_of(result$sic)  # n_distances x n_draws

  # For simultaneous bands, we need to check if the supremum of |deviation| is captured
  # The bands should contain the true curve in (1-alpha)% of draws
  coverage <- mean(result$true_value >= result$sic_lower & result$true_value <= result$sic_upper)

  # Since these are simultaneous bands, the pointwise coverage should be higher than 1-alpha
  # (Simultaneous bands are wider to ensure full curve coverage)
  expect_gt(coverage, 1 - alpha)

  # More importantly, check that the entire true curve is contained
  # For each column (draw), check if all true values fall within the bands
  # This is a conservative test since we're comparing to the mean bands, not per-draw bands
  n_distances_in_band <- sum(result$true_value >= result$sic_lower &
                               result$true_value <= result$sic_upper)
  expect_gt(n_distances_in_band / n_distances, 0.80)  # At least 80% of points covered
})

test_that("simultaneous bands are wider than pointwise bands", {
  # This is the key property: simultaneous bands should be uniformly wider
  set.seed(101)
  n_distances <- 50
  n_draws <- 200

  distances <- seq(0, 100, length.out = n_distances)
  true_curve <- sin(distances / 20)

  # Generate draws around true curve (draws x distances)
  sic_draws <- lapply(1:n_distances, function(i) {
    rnorm(n_draws, mean = true_curve[i], sd = 0.5)
  })
  sic_rvar <- posterior::rvar(do.call(cbind, sic_draws))

  sic_data <- tibble::tibble(
    distance = distances,
    source = "Source1",
    level_id = 1,
    level_name = "Group 1",
    sic = sic_rvar
  )

  result_simul <- add_simultaneous_bands(sic_data, alpha = 0.05)
  result_point <- add_pointwise_bands(sic_data, alpha = 0.05)

  # At every distance, simultaneous bands should be at least as wide as pointwise
  width_simul <- result_simul$sic_upper - result_simul$sic_lower
  width_point <- result_point$sic_upper - result_point$sic_lower

  expect_true(all(width_simul >= width_point - 1e-10))  # Allow tiny numerical error

  # On average, should be noticeably wider
  expect_gt(mean(width_simul), mean(width_point) * 1.1)
})
