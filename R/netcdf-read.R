# Reading values from a netCDF variable: `start` and `count` are built in the
# variable's own dimension order from the indices chosen on each axis, and
# ncdf4 reads exactly that block. Shared by point and buffer extraction.

.build_start_count <- function(var_dims, index_list,
                               singleton_dims = character(0)) {
  # Build start/count in the exact dimension order of the variable.
  # Every dimension must be explicitly indexed, except known singletons.

  start <- integer(length(var_dims))
  count <- integer(length(var_dims))

  for (i in seq_along(var_dims)) {
    dim_name <- var_dims[i]
    idx <- index_list[[dim_name]]

    if (!is.null(idx)) {
      start[i] <- as.integer(idx$start)
      count[i] <- as.integer(idx$count)
    } else if (dim_name %in% singleton_dims) {
      start[i] <- 1L
      count[i] <- 1L
    } else {
      stop(
        "No extraction index was assigned to dimension `",
        dim_name, "`. Refusing to read the full dimension ",
        "because this could return an incorrect value.",
        call. = FALSE
      )
    }
  }

  list(start = start, count = count)
}
