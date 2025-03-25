print_vb <- function(str,verbose) {
  if(verbose) print(str)
}

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
