# SHADE: Spatial Hierarchical Asymmetry via Directional Estimation

<!-- badges: start -->

[![tests](https://github.com/jeliason/SHADE/actions/workflows/tests.yaml/badge.svg)](https://github.com/jeliason/SHADE/actions/workflows/tests.yaml)

<!-- badges: end -->

SHADE is an R package for modeling asymmetric spatial associations between cell types in tissue images using a multilevel Bayesian framework.

Spatial relationships in tissue microenvironments are often directional—for example, immune cells may localize near tumor cells as part of immune surveillance, while tumor cell positions may be governed more by vasculature and tissue structure. These associations are not necessarily reciprocal, and conventional spatial models that assume symmetry fail to capture this critical asymmetry.

SHADE addresses this by modeling conditional spatial intensity—that is, how the presence of one or more *source* cell types predict the spatial distribution of a *target* cell type. These relationships are summarized using Spatial Interaction Curves (SICs), which quantify how the expected density of the target cell type varies with distance from the source. SICs offer a smooth, interpretable, and biologically grounded representation of directional interactions across spatial scales.

Additionally, biological data from spatial imaging studies are often hierarchically structured, with tissue images nested within patients and patients within cohorts. SHADE explicitly models this structure, enabling partial pooling and uncertainty quantification across levels to improve inference in sparse or heterogeneous datasets.

The model is implemented in Stan and uses a logistic regression approximation to estimate conditional intensity functions efficiently, even for large high-resolution datasets.

For more technical details and case studies, please see [our preprint](https://doi.org/10.1101/2025.06.24.661393) and [our repo](https://github.com/jeliason/shade_paper_code) for reproducing analyses and figures from that manuscript.

## Installation

To install from Github, first install CmdStanr and CmdStan:

``` r
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()
```

For problems installing CmdStanr and/or CmdStan, please see their documentation: <https://mc-stan.org/cmdstanr/>.

Then install SHADE:

``` r
devtools::install_github("jeliason/SHADE")
```

## Basic Usage

``` r
library(SHADE)
library(dplyr)
library(ggplot2)
library(posterior)

# Input data must include:
# - x, y: spatial coordinates of cells
# - cell_type: a factor or character vector of cell type labels
# - image_id: identifier for each image or sample
# - patient_metadata: a data.frame with columns 'Spot', 'Patient', and 'Group'

# Example: simulate spatial data with directional associations
# (In practice, you would load your own spatial imaging data)

set.seed(2025)
out <- simulate_spatial_data(
  n_images = 8,
  n_patients = 4,
  n_groups = 2,
  cell_types = c("tumor", "immune", "stroma"),
  target_type = "immune"  # immune cells attracted to tumor cells
)

coords <- out$data
# Create patient metadata for hierarchical modeling
patient_metadata <- data.frame(
  Spot = unique(coords$image_id),
  Patient = rep(c("pt_1", "pt_2"), each = 2),
  Group = rep(c("group_A", "group_B"), each = 2)
)

# Visualize example spatial pattern
coords %>%
  filter(image_id == "img_1") %>%
  ggplot(aes(x, y, color = cell_type)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Example Spatial Pattern")

# Prepare data for SHADE model
prep <- prepare_spatial_model_data(
  x = coords$x,
  y = coords$y,
  cell_type = coords$cell_type,
  image_id = coords$image_id,
  patient_metadata = patient_metadata,
  type_idx = 2,  # Index of target cell type (immune)
  n_basis_functions = 3
)

# Fit the SHADE model using Stan
fit <- run_SHADE_model(
  prep$stan_data,
  method = "variational",
  draws = 1e3,
  threads = 2
)

# Extract and summarize posterior estimates
rvars <- as_draws_rvars(fit$draws())
summary <- summarise_draws(rvars$beta_global)
print(summary)

# Plot Spatial Interaction Curves (SICs)
plot_spatial_interaction_curves(fit, prep, distance_range = c(0, 100))
```

For a complete end-to-end example, see the file `vignettes/Introduction.qmd`.

## Package Overview

SHADE includes tools for:

-   Preparing spatial point pattern data for hierarchical analysis
-   Defining flexible spatial interaction features via basis functions
-   Fitting multilevel spatial point process models using Stan
-   Summarizing posterior distributions of interaction curves
-   Comparing results across images, patients, and groups

## Citation and References

If you use SHADE in your work, please cite the accompanying [preprint](https://doi.org/10.1101/2025.06.24.661393):

``` bibtex
@misc{
  title = {{{SHADE}}: {{A Multilevel Bayesian Approach}} to {{Modeling Directional Spatial Associations}} in {{Tissues}}},
  author = {Eliason, Joel and Peruzzi, Michele and Rao, Arvind},
  year = {2025},
  month = jun,
  publisher = {Cold Spring Harbor Laboratory},
  doi = {10.1101/2025.06.24.661393},
}
```

## License

SHADE is released under the MIT license.

## Contact

Joel Eliason

[joeleliason.com](https://joeleliason.com)
