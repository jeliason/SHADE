#' @title Fit the SHADE model.
#' @export
#' @family models
#' @description Fit the SHADE Stan model and return fit
#' @return Model fit
#' @param data_stan list of arguments needed for SHADE model
#' @param ... Named arguments to the `sample()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
run_SHADE_model <- function(data_stan, ...) {
  model <- instantiate::stan_package_model(name = "SHADE", package = "SHADE")
  model$sample(data = data_stan, ...)
}