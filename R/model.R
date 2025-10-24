#' @title Fit the SHADE model.
#' @export
#' @family models
#' @description Fit the SHADE Stan model and return fit
#' @return Model fit
#' @param data_stan list of arguments needed for SHADE model
#' @param method Character string specifying the inference method: "sample" for MCMC sampling (default) or "variational" for variational inference
#' @param parameterization Character string specifying the parameterization: "centered" for centered parameterization (default) or "non_centered" for non-centered parameterization. Non-centered parameterization can improve sampling efficiency when hierarchical variance parameters are small.
#' @param ... Named arguments to the `sample()` or `variational()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_SHADE_model <- function(data_stan, method = "sample", parameterization = "centered", ...) {
  # Select model based on parameterization
  if (parameterization == "centered") {
    model <- instantiate::stan_package_model(name = "SHADE", package = "SHADE")
  } else if (parameterization == "non_centered") {
    model <- instantiate::stan_package_model(name = "SHADE_ncp", package = "SHADE")
  } else {
    stop("parameterization must be either 'centered' or 'non_centered'")
  }

  if (method == "sample") {
    model$sample(data = data_stan, ...)
  } else if (method == "variational") {
    model$variational(data = data_stan, ...)
  } else {
    stop("Method must be either 'sample' or 'variational'")
  }
}

#' @title Generate predictions from fitted SHADE model.
#' @export
#' @family models
#' @description Generate predictions
#' @return Model fit
#' @param data_gq list of arguments needed for SHADE generated-quantities
#' @param draws_fit Model fit from `run_SHADE_model`
#' @param ... Named arguments to the `generate_quantities()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_SHADE_gq <- function(data_gq, draws_fit, ...) {
  model <- instantiate::stan_package_model(name = "SHADE_gq", package = "SHADE")
  
  fit_gq <- model$generate_quantities(draws_fit, data = data_gq, ...)
}
