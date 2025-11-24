print_vb <- function(str,verbose) {
  if(verbose) print(str)
}

#' Write Large JSON Object in Chunks
#'
#' Efficiently writes a list-based data object (used for Stan input, etc.) to disk in JSON format.
#' Supports numeric scalars, vectors, and matrices. Large vectors and matrices are written in chunks
#' to reduce memory usage and allow partial serialization of large arrays.
#'
#' @param data A named list containing elements to serialize. Each element must be a numeric scalar, numeric vector, or numeric matrix.
#' @param file_path Path to the output JSON file.
#' @param chunk_size Number of elements per write chunk when writing long vectors or matrix rows. Default is 1000.
#' @param verbose Logical; if TRUE (default), prints progress messages for each top-level element and chunk.
#'
#' @return Invisibly returns NULL. Writes the JSON file to disk.
#'
#' @details
#' The function writes JSON line-by-line to avoid creating large in-memory representations of large vectors or matrices.
#' Each element in the input list becomes a top-level key in the resulting JSON object.
#' For vectors and matrices, only numeric types are supported.
#'
#' Unsupported types (e.g., lists, data frames, arrays > 2D) will trigger an error.
#'
#' @examples
#' \dontrun{
#' stan_data <- list(
#'   scalar = 42,
#'   vector = rnorm(10000),
#'   matrix = matrix(runif(100), nrow = 10)
#' )
#' write_json_chunked(stan_data, "output.json")
#' }
#'
#' @export
write_json_chunked <- function(data, file_path, chunk_size = 1000,verbose=TRUE) {
  # Open a connection to the file
  file_conn <- file(file_path, "w")

  # Write the opening bracket
  writeLines("{", file_conn)

  # Loop over each element in the list
  elements <- names(data)
  for (i in seq_along(elements)) {
    element_name <- elements[i]
    element_value <- data[[element_name]]

    if(verbose) print(paste0("Element: ", element_name))

    # Write the element name as a JSON key
    writeLines(paste0('\t"', element_name, '": '), file_conn, sep = "")

    # Check the type of element and write accordingly
    if (is.numeric(element_value) && length(element_value) == 1) {
      # Scalar: preserve full precision
      writeLines(jsonlite::toJSON(element_value, auto_unbox = TRUE, digits = NA), file_conn, sep = "")

    } else if (is.vector(element_value)) {
      # Vector: Write in chunks with full precision, no extra brackets
      writeLines("[", file_conn, sep = "")
      for (j in seq(1, length(element_value), by = chunk_size)) {
        print(paste0("Chunk: ",j))
        chunk <- element_value[j:min(j + chunk_size - 1, length(element_value))]
        if (j + chunk_size - 1 < length(element_value)) {
          chunk_json <- paste0("\t",paste0(chunk,collapse = ","),",")
          writeLines(chunk_json, file_conn)
        } else {
          chunk_json <- paste0(chunk,collapse = ",")
          writeLines(chunk_json, file_conn,sep="")
        }
      }
      writeLines("]", file_conn,sep="")

    } else if (is.matrix(element_value)) {
      # Matrix: Write each row in chunks with full precision
      writeLines("[", file_conn, sep = "\n")
      for (j in 1:nrow(element_value)) {
        row_json <- jsonlite::toJSON(as.vector(element_value[j, ]), auto_unbox = TRUE, digits = NA)

        # Write each row with comma if not the last row
        if (j < nrow(element_value)) {
          writeLines(paste0("\t\t",row_json, ","), file_conn)
        } else {
          writeLines(paste0("\t\t",row_json), file_conn)
        }
      }
      writeLines("\t]", file_conn,sep="")

    } else {
      stop("Unsupported data type. Only scalars, vectors, and matrices are supported.")
    }

    # Add a comma and newline after each element except the last one
    if (i < length(elements)) {
      writeLines(",", file_conn)
    } else {
      writeLines("",file_conn)
    }
  }

  # Close the JSON structure and the file connection
  writeLines("}", file_conn)
  close(file_conn)
}

#' Plot Group-Level Spatial Interaction Curves (SICs)
#'
#' Creates a faceted plot showing the estimated group-level Spatial Interaction Curves (SICs) from a SHADE model fit.
#' The SICs represent how the expected density of the target cell type varies with distance
#' from each source cell type, across different groups if present. Each facet shows the SIC for
#' a specific source-target cell type interaction. Can optionally compare with true curves
#' if simulation parameters are provided.
#'
#' @param fit A CmdStan model fit object returned by `run_SHADE_model()`.
#' @param prep The preprocessing object returned by `prepare_spatial_model_data()`.
#' @param true_params Optional list of true simulation parameters (from `simulate_spatial_data()`) for comparison.
#' @param distance_range Numeric vector of length 2 specifying the range of distances to plot (default: c(0, 100)).
#' @param alpha Significance level for credible intervals (default: 0.05 for 95% intervals).
#' @param resolution Number of distance points to evaluate (default: 81).
#'
#' @return A ggplot object showing the SICs with credible intervals, faceted by source-target interactions.
#' @export
plot_spatial_interaction_curves <- function(fit, prep, true_params = NULL, distance_range = c(0, 100), alpha = 0.05, resolution = 81) {

  # Create distance sequence
  distance_seq <- seq(distance_range[1], distance_range[2], length.out = resolution)

  # Extract group-level SICs with simultaneous bands using new API
  plot_data <- extract_group_sics(
    fit = fit,
    prep = prep,
    distance_seq = distance_seq,
    bands = "simultaneous",
    alpha = alpha,
    keep_rvar = FALSE
  )

  # Get cell type information for labeling
  all_types <- prep$metadata$types
  target_type <- all_types[prep$stan_data$num_types]

  # Rename columns for compatibility with plotting code
  plot_data <- plot_data %>%
    dplyr::mutate(
      x = distance,
      mean_val = sic_mean,
      lower = sic_lower,
      upper = sic_upper,
      group = level_name,
      source_type = source,
      interaction = paste0(source, " → ", target_type)
    )
  
  # Add true curves if provided
  if (!is.null(true_params) && !is.null(true_params$betas_global)) {
    # Create design matrix using true basis functions
    true_basis <- true_params$basis_functions
    x_des_true <- lapply(true_basis, function(pot) pot(x_seq)) %>% 
      do.call(cbind, .)
    
    # Compute true curves for each group and source type
    true_data <- lapply(1:n_groups, function(group_i) {
      lapply(1:n_source_types, function(source_i) {
        coeff_start <- (source_i - 1) * n_basis + 1
        coeff_end <- source_i * n_basis
        coeff_indices <- coeff_start:coeff_end + 1  # +1 to skip intercept
        
        beta_true <- true_params$betas_global[coeff_indices, group_i]
        lp_true <- as.vector(x_des_true %*% beta_true)
        
        data.frame(
          x = x_seq,
          true_val = lp_true,
          group = paste0("Group ", group_i),
          source_type = source_types[source_i],
          interaction = paste0(source_types[source_i], " → ", target_type)
        )
      }) %>% 
        do.call(rbind, .)
    }) %>% 
      do.call(rbind, .)
    
    # Combine estimated and true data for plotting
    plot_data_lines <- plot_data %>%
      dplyr::select(x, group, source_type, interaction, estimated = mean_val) %>%
      dplyr::left_join(true_data, by = c("x", "group", "source_type", "interaction")) %>%
      tidyr::pivot_longer(cols = c(estimated, true_val), 
                         names_to = "type", values_to = "value") %>%
      dplyr::mutate(type = dplyr::case_when(
        type == "estimated" ~ "Estimated",
        type == "true_val" ~ "True",
        TRUE ~ type
      ))
    
    # Create plot with both estimated and true curves
    p <- ggplot2::ggplot() +
      ggplot2::geom_ribbon(data = plot_data,
                          ggplot2::aes(x = x, ymin = lower, ymax = upper, fill = group),
                          alpha = 0.2) +
      ggplot2::geom_line(data = plot_data_lines,
                        ggplot2::aes(x = x, y = value, color = group, linetype = type),
                        size = 1) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dotted", alpha = 0.5) +
      ggplot2::facet_wrap(~ interaction, scales = "free_y") +
      ggplot2::scale_color_manual(
        values = c("Group 1" = "#e6194b", "Group 2" = "#3cb44b", "Group 3" = "#3cb4f4")
      ) +
      ggplot2::scale_fill_manual(
        values = c("Group 1" = "#e6194b", "Group 2" = "#3cb44b", "Group 3" = "#3cb4f4")
      ) +
      ggplot2::scale_linetype_manual(
        values = c("Estimated" = "solid", "True" = "dashed")
      ) +
      ggplot2::labs(
        x = "Distance (microns)",
        y = "Spatial Interaction Curve",
        color = "Group",
        fill = "Group",
        linetype = "Curve Type",
        title = "Spatial Interaction Curves (SICs): Estimated vs True"
      ) +
      ggplot2::guides(fill = "none") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5),
        strip.text = ggplot2::element_text(size = 10, face = "bold")
      )
    
  } else {
    # Create plot with estimated curves only
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = mean_val, color = group, fill = group)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
      ggplot2::facet_wrap(~ interaction, scales = "free_y") +
      ggplot2::scale_color_manual(
        values = c("Group 1" = "#e6194b", "Group 2" = "#3cb44b", "Group 3" = "#3cb4f4")
      ) +
      ggplot2::scale_fill_manual(
        values = c("Group 1" = "#e6194b", "Group 2" = "#3cb44b", "Group 3" = "#3cb4f4")
      ) +
      ggplot2::labs(
        x = "Distance (microns)",
        y = "Spatial Interaction Curve",
        color = "Group",
        fill = "Group",
        title = "Spatial Interaction Curves (SICs)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5),
        strip.text = ggplot2::element_text(size = 10, face = "bold")
      )
  }
  
  return(p)
}
