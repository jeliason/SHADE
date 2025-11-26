# Example: Extracting and Analyzing Spatial Interaction Curves (SICs)
#
# This example demonstrates the new SIC extraction API that allows you to:
# 1. Extract SICs at different hierarchical levels (image, patient, group)
# 2. Compute credible bands (simultaneous or pointwise)
# 3. Perform custom analyses with posterior distributions

library(SHADE)
library(dplyr)
library(ggplot2)
library(posterior)  # For working with rvar objects

# ==============================================================================
# 1. SIMULATE DATA AND FIT MODEL
# ==============================================================================

cat("Simulating spatial data...\n")

# Simulate spatial data with 3 cell types across 4 images (2 patients, 2 groups)
set.seed(456)

# Create data for 4 images
n_per_image <- 150
n_images <- 4

x <- c(
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100)
)

y <- c(
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100),
  runif(n_per_image, 0, 100)
)

cell_type <- factor(c(
  sample(c("CD8", "CD4", "Tumor"), n_per_image, replace = TRUE),
  sample(c("CD8", "CD4", "Tumor"), n_per_image, replace = TRUE),
  sample(c("CD8", "CD4", "Tumor"), n_per_image, replace = TRUE),
  sample(c("CD8", "CD4", "Tumor"), n_per_image, replace = TRUE)
))

image_id <- factor(c(
  rep("Image_1", n_per_image),
  rep("Image_2", n_per_image),
  rep("Image_3", n_per_image),
  rep("Image_4", n_per_image)
))

# Patient metadata: 2 patients, 2 groups
patient_metadata <- data.frame(
  Spot = c("Image_1", "Image_2", "Image_3", "Image_4"),
  Patient = c("P1", "P1", "P2", "P2"),
  Group = c("Control", "Control", "Treatment", "Treatment")
)

cat("Preparing spatial model data...\n")

# Prepare data
prep <- prepare_spatial_model_data(
  x = x,
  y = y,
  cell_type = cell_type,
  image_id = image_id,
  patient_metadata = patient_metadata,
  n_dummy = 500,
  type_idx = 1,  # Model CD8 as target
  n_basis_functions = 4,
  max_dist = 75
)

cat("\nFitting SHADE model (this may take a minute)...\n")

# Fit model (use variational inference for speed in example)
fit <- run_SHADE_model(
  prep$stan_data,
  method = "variational",
  seed = 123
)

cat("Model fitting complete!\n\n")

# ==============================================================================
# 2. BASIC USAGE: CONVENIENCE WRAPPERS
# ==============================================================================

cat("======================================================\n")
cat("Example 1: Simple Extraction with Convenience Wrappers\n")
cat("======================================================\n\n")

# Extract group-level SICs with simultaneous bands (most common use case)
cat("Extracting group-level SICs with 95% simultaneous bands...\n")
sics_group <- extract_group_sics(
  fit = fit,
  prep = prep,
  bands = "simultaneous",
  alpha = 0.05
)

print(head(sics_group, 10))
cat("\n")

# Plot
cat("Creating plot...\n")
p1 <- ggplot(sics_group, aes(x = distance, y = sic_mean, color = level_name, fill = level_name)) +
  geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~source, scales = "free_y") +
  labs(
    title = "Group-Level SICs with 95% Simultaneous Bands",
    x = "Distance (μm)",
    y = "Spatial Interaction Curve",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)

# ==============================================================================
# 3. COMPARING BAND TYPES
# ==============================================================================

cat("\n======================================================\n")
cat("Example 2: Comparing Simultaneous vs Pointwise Bands\n")
cat("======================================================\n\n")

# Extract with both band types
sics_simultaneous <- extract_group_sics(fit, prep, bands = "simultaneous", alpha = 0.1)
sics_pointwise <- extract_group_sics(fit, prep, bands = "pointwise", alpha = 0.1)

# Combine for comparison
sics_comparison <- bind_rows(
  sics_simultaneous %>% mutate(band_type = "Simultaneous (90%)"),
  sics_pointwise %>% mutate(band_type = "Pointwise (90%)")
)

# Plot comparison
p2 <- ggplot(sics_comparison, aes(x = distance, y = sic_mean, color = band_type, fill = band_type)) +
  geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_grid(level_name ~ source, scales = "free_y") +
  scale_color_manual(values = c("Simultaneous (90%)" = "#D95F02", "Pointwise (90%)" = "#1B9E77")) +
  scale_fill_manual(values = c("Simultaneous (90%)" = "#D95F02", "Pointwise (90%)" = "#1B9E77")) +
  labs(
    title = "Simultaneous vs Pointwise Credible Bands",
    subtitle = "Simultaneous bands are wider but provide stronger guarantees",
    x = "Distance (μm)",
    y = "SIC",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)

# Compute band width comparison
band_widths <- sics_comparison %>%
  mutate(width = sic_upper - sic_lower) %>%
  group_by(source, level_name, band_type) %>%
  summarize(
    mean_width = mean(width),
    max_width = max(width),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = band_type,
    values_from = c(mean_width, max_width)
  )

cat("\nBand width comparison:\n")
print(band_widths)

# ==============================================================================
# 4. ADVANCED: WORKING WITH POSTERIOR DISTRIBUTIONS
# ==============================================================================

cat("\n\n======================================================\n")
cat("Example 3: Advanced Analysis with Posterior Draws\n")
cat("======================================================\n\n")

# Extract raw posteriors (as rvar objects) for custom analysis
cat("Extracting raw posterior distributions...\n")
sics_posterior <- compute_sic_posterior(
  fit = fit,
  prep = prep,
  level = "group",
  distance_seq = seq(0, 100, by = 2)  # Coarser grid for efficiency
)

cat("Structure of posterior data:\n")
print(head(sics_posterior, 3))
cat("\n")

# Example custom analysis: Compute probability that SIC > 0 at each distance
cat("Computing probability that SIC > 0 at each distance...\n")
prob_positive <- sics_posterior %>%
  mutate(
    prob_positive = sapply(sic, function(rv) mean(rv > 0))
  ) %>%
  select(distance, source, level_name, prob_positive)

cat("\nProbability that SIC > 0:\n")
print(head(prob_positive, 10))

# Plot probability of positive interaction
p3 <- ggplot(prob_positive, aes(x = distance, y = prob_positive, color = level_name)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = c(0.025, 0.975), linetype = "dotted", alpha = 0.3) +
  facet_wrap(~source) +
  labs(
    title = "Probability of Positive Spatial Interaction",
    x = "Distance (μm)",
    y = "P(SIC > 0)",
    color = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p3)

# Example: Compute difference between groups
cat("\nComputing group differences...\n")
group_comparison <- sics_posterior %>%
  dplyr::select(distance, source, level_name, sic) %>%
  tidyr::pivot_wider(names_from = level_name, values_from = sic) %>%
  mutate(
    diff = `Group 2` - `Group 1`,
    diff_mean = E(diff),
    prob_greater = sapply(diff, function(rv) mean(rv > 0))
  )

cat("Posterior probability that Treatment > Control:\n")
print(group_comparison %>% dplyr::select(distance, source, prob_greater) %>% head(10))

# ==============================================================================
# 5. PATIENT-LEVEL ANALYSIS
# ==============================================================================

cat("\n\n======================================================\n")
cat("Example 4: Patient-Level SICs\n")
cat("======================================================\n\n")

# Extract patient-level SICs
cat("Extracting patient-level SICs...\n")
sics_patient <- extract_patient_sics(
  fit = fit,
  prep = prep,
  bands = "pointwise",  # Pointwise bands for patient level
  alpha = 0.1
)

# Plot patient-level SICs
p4 <- ggplot(sics_patient, aes(x = distance, y = sic_mean, color = level_name, fill = level_name)) +
  geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~source, scales = "free_y") +
  labs(
    title = "Patient-Level SICs with 90% Pointwise Bands",
    x = "Distance (μm)",
    y = "SIC",
    color = "Patient",
    fill = "Patient"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p4)

# ==============================================================================
# 6. CUSTOM DISTANCE RANGES
# ==============================================================================

cat("\n\n======================================================\n")
cat("Example 5: Truncated Distance Range\n")
cat("======================================================\n\n")

# Often we want to exclude very short distances (cell overlap artifacts)
# and limit to biologically relevant distances

MIN_DISTANCE <- 30  # microns
MAX_DISTANCE <- 100

cat(sprintf("Extracting SICs for distances %d-%d μm...\n", MIN_DISTANCE, MAX_DISTANCE))

sics_truncated <- extract_group_sics(
  fit = fit,
  prep = prep,
  distance_seq = seq(MIN_DISTANCE, MAX_DISTANCE, by = 1),
  bands = "simultaneous",
  alpha = 0.1
)

p5 <- ggplot(sics_truncated, aes(x = distance, y = sic_mean, color = level_name, fill = level_name)) +
  geom_ribbon(aes(ymin = sic_lower, ymax = sic_upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~source, scales = "free_y") +
  labs(
    title = sprintf("Group-Level SICs (Truncated: %d-%d μm)", MIN_DISTANCE, MAX_DISTANCE),
    subtitle = "90% simultaneous credible bands",
    x = "Distance (μm)",
    y = "SIC",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p5)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n\n======================================================\n")
cat("Summary of SIC Extraction API\n")
cat("======================================================\n\n")

cat("Three main approaches:\n\n")

cat("1. CONVENIENCE WRAPPERS (recommended for most users):\n")
cat("   - extract_group_sics()   : Group-level SICs\n")
cat("   - extract_patient_sics() : Patient-level SICs\n")
cat("   - extract_image_sics()   : Image-level SICs\n")
cat("   Each automatically adds credible bands\n\n")

cat("2. CORE FUNCTION (for custom analyses):\n")
cat("   - compute_sic_posterior() : Returns rvar posteriors\n")
cat("   Allows custom summarization and analysis\n\n")

cat("3. BAND FUNCTIONS (for manual pipeline):\n")
cat("   - add_simultaneous_bands() : Simultaneous credible bands\n")
cat("   - add_pointwise_bands()    : Pointwise credible bands\n")
cat("   Can be applied to any rvar data\n\n")

cat("Documentation:\n")
cat("   ?compute_sic_posterior\n")
cat("   ?extract_group_sics\n")
cat("   ?add_simultaneous_bands\n")
cat("   ?add_pointwise_bands\n\n")
