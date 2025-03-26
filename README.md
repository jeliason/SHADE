# SHADE: Spatial Hierarchical Asymmetry via Directional Estimation

SHADE is an R package for modeling asymmetric spatial associations between cell types in tissue images using a hierarchical Bayesian framework. It is designed for high-resolution spatial data from technologies such as multiplexed immunofluorescence (mIF) and CODEX.

SHADE captures directional spatial interactions using nonparametric spatial interaction curves (SICs) and models variation across biological levels, including images, patients, and cohorts. The model is implemented in Stan and fit using the logistic regression approximation for conditional intensity functions.

## Installation

To install the development version locally:

```r
# Install dependencies
install.packages(c("devtools", "cmdstanr"))

# Install SHADE from local source or GitHub
devtools::install_local("path/to/SHADE")  # or install_github("jeliason/SHADE")
```

During development, load with:

```r
devtools::load_all()
```

## Basic Usage

```r
library(SHADE)

# Input data must include:
# - x, y: spatial coordinates of cells
# - cell_type: a factor or character vector of cell type labels
# - image_id: identifier for each image or sample
# - patient_metadata: a data.frame with columns 'Spot', 'Patient', and 'Group'

# Example: preprocess spatial data into SHADE model input
prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$type,
  image_id = coords$image_id,
  patient_metadata = metadata,
  type_idx = 3  # Index of target cell type
)

# Fit the SHADE model using Stan
fit <- run_SHADE_model(
  prep$stan_data,
  chains = 4,
  iter_warmup = 500,
  iter_sampling = 1000
)

# Summarize posterior estimates
library(posterior)
rvars <- as_draws_rvars(fit$draws())
summary <- summarise_draws(rvars$beta_global)
```

For a complete end-to-end example, see the file `vignettes/Introduction.Rmd`.

## Package Overview

SHADE includes tools for:

- Preparing spatial point pattern data for hierarchical analysis
- Defining flexible spatial interaction features via radial basis functions
- Fitting multilevel spatial point process models using Stan
- Summarizing posterior distributions of interaction curves
- Comparing results across images, patients, and groups

## Citation and References

If you use SHADE in your work, please cite the accompanying paper (forthcoming). Key references include:

- Baddeley et al. (2014), "Logistic regression for spatial Gibbs point processes"
- Grabarnik and Särkkä (2009), "Modelling the spatial structure of forest stands by multivariate point processes with hierarchical interactions"

## License

SHADE is released under the MIT license.