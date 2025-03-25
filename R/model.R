#' @title Fit the SHADE model.
#' @export
#' @family models
#' @description Fit the SHADE Stan model and return posterior summaries.
#' @return A data frame of posterior summaries.
#' @param data_stan list of arguments needed for SHADE model
#' @param ... Named arguments to the `sample()` method of CmdStan model
#'   objects: <https://mc-stan.org/cmdstanr/reference/model-method-sample.html>
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   run_bernoulli_model(y = c(1, 0, 1, 0, 1, 0, 0, 0, 0, 0))
#' }
run_SHADE_model <- function(data_stan, ...) {
  model <- instantiate::stan_package_model(
    name = "SHADE",
    package = "SHADE"
  )
  fit <- model$sample(data = data_stan, ...)
}
