# Reads one index window of a SpatProxy's variable with a single ncdf4 call.
# `window` is a named list with up to one element per axis (`lon`, `lat`,
# `depth`, `time`): the first and last index to read, in storage order. An
# axis left out is read whole. Returns the values, unpacked by ncdf4 (scale,
# offset, and fill values as NA), as a 4D array in the fixed order lon, lat,
# depth, time, with length 1 for an axis the variable doesn't have.
.read_proxy_window <- function(x, window = list()) {
  variable <- x@variable
  axes <- c("lon", "lat", "depth", "time")
  unknown <- setdiff(names(window), axes)
  if (length(unknown)) {
    rlang::abort(glue::glue(
      "Unknown axis in the read window: {toString(unknown)}."
    ), class = "ocean3d_error_read_window", call = NULL)
  }
  index_list <- list()
  shape <- stats::setNames(rep(1L, 4L), axes)
  for (axis in axes) {
    dim <- variable$axis_dims[[axis]]
    if (is.null(dim)) next
    n <- variable$dim_lengths[[dim]]
    span <- if (is.null(window[[axis]])) c(1L, n) else window[[axis]]
    if (!is.numeric(span) || length(span) != 2L || anyNA(span) ||
          span[1L] < 1L || span[2L] > n || span[1L] > span[2L]) {
      rlang::abort(glue::glue(
        "The read window for `{axis}` must be a first and last index ",
        "between 1 and {n}."
      ), class = "ocean3d_error_read_window", call = NULL)
    }
    shape[[axis]] <- as.integer(span[2L] - span[1L] + 1L)
    index_list[[dim]] <- list(start = span[1L], count = shape[[axis]])
  }
  sc <- .build_start_count(variable$dims, index_list, variable$singleton_dims)

  opened <- .open_source(list(source = x@source$path, type = "path"))
  on.exit(.close_source(opened), add = TRUE)
  values <- ncdf4::ncvar_get(opened$nc, variable$name, start = sc$start,
                             count = sc$count, collapse_degen = FALSE)
  # Put the axes in fixed order; other dimensions have length 1.
  present <- unlist(variable$axis_dims[axes])
  dim(values) <- sc$count
  values <- aperm(values, c(match(present, variable$dims),
                            which(!variable$dims %in% present)))
  dim(values) <- unname(shape)
  values
}
