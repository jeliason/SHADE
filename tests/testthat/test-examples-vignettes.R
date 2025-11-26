# Tests for example scripts and vignettes
#
# These are integration tests that verify the example scripts and vignettes
# run without error. They are skipped on CRAN due to resource requirements.

# Helper function to check if CmdStan is available
cmdstan_available <- function() {
  tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
}

# Helper function to check if quarto is available
quarto_available <- function() {
  tryCatch({
    result <- system2("quarto", "--version", stdout = TRUE, stderr = TRUE)
    length(result) > 0
  }, error = function(e) FALSE)
}

# Helper function to check if Stan models are compiled
stan_models_available <- function() {
  tryCatch({
    model <- instantiate::stan_package_model(name = "SHADE", package = "SHADE")
    !is.null(model)
  }, error = function(e) FALSE)
}

# Helper to find example files (works both during development and after install)
find_example_file <- function(filename) {
  # Try package inst directory first (for installed package)
  path <- system.file("examples", filename, package = "SHADE")
  if (path != "" && file.exists(path)) return(path)

  # Try relative path from test directory (for development)
  path <- file.path(getwd(), "..", "..", "examples", filename)
  if (file.exists(path)) return(normalizePath(path))

  # Try from package root
  path <- file.path(getwd(), "..", "..", "inst", "examples", filename)
  if (file.exists(path)) return(normalizePath(path))

  NULL
}

# Helper to find vignette files
find_vignette_file <- function(filename) {
  # Try relative path from test directory (for development)
  path <- file.path(getwd(), "..", "..", "vignettes", filename)
  if (file.exists(path)) return(normalizePath(path))

  NULL
}

test_that("covariate_example.R runs without error", {
  skip_on_cran()
  skip_if_not(cmdstan_available(), "CmdStan not available")

  example_path <- find_example_file("covariate_example.R")
  skip_if(is.null(example_path), "covariate_example.R not found")

  # Capture output to suppress cat() statements
  expect_no_error({
    suppressMessages(
      capture.output(
        source(example_path, local = TRUE),
        type = "output"
      )
    )
  })
})

test_that("sic_extraction_example.R runs without error", {
  skip_on_cran()
  skip_if_not(cmdstan_available(), "CmdStan not available")
  skip_if_not(stan_models_available(), "Stan models not compiled")

  example_path <- find_example_file("sic_extraction_example.R")
  skip_if(is.null(example_path), "sic_extraction_example.R not found")

  # This test takes longer as it fits a model with variational inference
  expect_no_error({
    suppressMessages(
      suppressWarnings(
        capture.output(
          source(example_path, local = TRUE),
          type = "output"
        )
      )
    )
  })
})

test_that("Introduction vignette renders successfully", {
  skip_on_cran()
  skip_if_not(cmdstan_available(), "CmdStan not available")
  skip_if_not(stan_models_available(), "Stan models not compiled")
  skip_if_not(quarto_available(), "Quarto not available")

  vignette_path <- find_vignette_file("Introduction.qmd")
  skip_if(is.null(vignette_path), "Introduction.qmd not found")

  vignette_dir <- dirname(vignette_path)

  # Build vignette in its own directory
  # Use --cache-refresh to avoid stale cache paths from different locations
  withr::with_dir(vignette_dir, {
    # Use output-dir to place output in a temp directory
    output_dir <- tempdir()

    result <- system2(
      "quarto",
      args = c(
        "render", "Introduction.qmd",
        "--to", "html",
        "--output-dir", output_dir,
        "--cache-refresh"
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    exit_status <- attr(result, "status")
    output_file <- file.path(output_dir, "Introduction.html")

    # Check that rendering succeeded
    expect_true(
      is.null(exit_status) || exit_status == 0,
      info = paste("Quarto render failed with output:", paste(result, collapse = "\n"))
    )

    # Check that output file was created
    expect_true(
      file.exists(output_file),
      info = paste("Quarto output:", paste(result, collapse = "\n"))
    )

    # Clean up
    unlink(output_file)
  })
})
