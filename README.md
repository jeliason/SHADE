# SHADE: Spatial Hierarchical Asymmetry via Directional Estimation

<!-- badges: start -->
[![tests](https://github.com/jeliason/SHADE/actions/workflows/tests.yaml/badge.svg)](https://github.com/jeliason/SHADE/actions/workflows/tests.yaml)
<!-- badges: end -->

SHADE is an R package for modeling asymmetric spatial associations between cell types in tissue images using a multilevel Bayesian framework.

Spatial relationships in tissue microenvironments are often directional—for example, immune cells may localize near tumor cells as part of immune surveillance, while tumor cell positions may be governed more by vasculature and tissue structure. These associations are not necessarily reciprocal, and conventional spatial models that assume symmetry or analyze each image independently fail to capture this critical asymmetry.

SHADE addresses this by modeling conditional spatial intensity—that is, how the presence of one or more *source* cell types predict the spatial distribution of a *target* cell type. These relationships are summarized using Spatial Interaction Curves (SICs), which quantify how the expected density of the target cell type varies with distance from the source. SICs offer a smooth, interpretable, and biologically grounded representation of directional interactions across spatial scales.

Additionally, biological data from spatial imaging studies are often hierarchically structured, with tissue images nested within patients and patients within cohorts. SHADE explicitly models this structure, enabling partial pooling and uncertainty quantification across levels to improve inference in sparse or heterogeneous datasets.

The model is implemented in Stan and uses a logistic regression approximation to estimate conditional intensity functions efficiently, even for large high-resolution datasets.

For more technical details and case studies, please see [our preprint](https://doi.org/10.1101/2025.06.24.661393).

## Installation

To install from Github:

```r
devtools::install_github("jeliason/SHADE")
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
- Defining flexible spatial interaction features via basis functions
- Fitting multilevel spatial point process models using Stan
- Summarizing posterior distributions of interaction curves
- Comparing results across images, patients, and groups

## Citation and References

If you use SHADE in your work, please cite the accompanying [preprint](https://doi.org/10.1101/2025.06.24.661393):

```bibtex
@misc{,
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