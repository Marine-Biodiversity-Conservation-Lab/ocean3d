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

.extract_array <- function(nc,
                           var,
                           dims,
                           lon_idx,
                           lat_idx,
                           time_idx = NULL,
                           depth_start = NULL,
                           depth_count = NULL) {
  index_list <- list()
  index_list[[dims$lon]] <- list(start = lon_idx, count = 1L)
  index_list[[dims$lat]] <- list(start = lat_idx, count = 1L)

  if (!is.null(dims$time) && !is.null(time_idx)) {
    index_list[[dims$time]] <- list(start = time_idx, count = 1L)
  }

  if (!is.null(dims$depth)) {
    if (!is.null(depth_start)) {
      # Explicit depth selection for 3D extraction methods.
      index_list[[dims$depth]] <- list(
        start = depth_start,
        count = depth_count
      )
    } else if (nc$dim[[dims$depth]]$len == 1L) {
      # A single depth level can be selected without ambiguity.
      index_list[[dims$depth]] <- list(
        start = 1L,
        count = 1L
      )
    } else {
      stop(
        "Variable `", var, "` has multiple depth levels ",
        "but no depth index was provided.",
        call. = FALSE
      )
    }
  }

  sc <- .build_start_count(
    dims$var_dims, index_list,
    singleton_dims = dims$singleton_dims
  )

  ncdf4::ncvar_get(
    nc = nc,
    varid = var,
    start = sc$start,
    count = sc$count
  )
}
