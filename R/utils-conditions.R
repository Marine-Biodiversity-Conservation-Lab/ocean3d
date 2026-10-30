# R/utils-conditions.R
#
# Custom conditions shared by the exported extraction functions.
#
# Each helper signals the error or warning for one state of the input to
# extract_to_point() and extract_to_buffer(). Keeping them here gives both
# functions the same wording. Condition classes are internal: they let tests
# catch a specific condition without matching its wording, and are not
# documented for users. Nothing here is exported.

# Error: a selected variable has a time dimension, but no date input was
# supplied (no `date` argument, or `date_col` is not a column of `data`).
# An explicit NA date is allowed; only a missing input is an error.
# Called once per source by extract_to_point() and extract_to_buffer().
.require_observation_time <- function(data, date_col, var, source_file,
                                      direct_input = FALSE) {
  named_column <- is.character(date_col) &&
    length(date_col) == 1L &&
    !is.na(date_col) &&
    nzchar(date_col)
  if (named_column && date_col %in% names(data)) return(invisible(NULL))

  correction <- if (direct_input) {
    "Supply the `date` argument."
  } else if (named_column) {
    glue::glue(
      "Column `{date_col}` is missing from `data`. ",
      "Supply that column or set `date_col` to an existing date column."
    )
  } else {
    "Set `date_col` to an existing date column in `data`."
  }
  file <- basename(source_file)

  rlang::abort(glue::glue(
    "netCDF variable `{var}` has a time dimension.\n",
    "File: {file}\n",
    "{correction}\n",
    "An explicit NA date is allowed; a missing date input is not."
  ), call = NULL)
}

# Warning, class `ocean3d_timezone_warning`: some text date-times have no time
# zone (no trailing Z) and are being read as UTC. `prepared` is the result of
# .prepare_time_input(). Raised once per extraction, even when
# `time_report = FALSE`.
.warn_time_zone <- function(prepared, input_name) {
  n_unzoned <- length(prepared$unzoned_rows)
  if (!n_unzoned) return(invisible(NULL))

  observations <- if (n_unzoned == 1L) {
    "observation has"
  } else {
    "observations have"
  }
  rlang::warn(glue::glue(
    "Time zone warning: {n_unzoned} {observations} times without a ",
    "time zone.\n",
    "Input: {input_name}. These times are being interpreted as UTC.\n",
    "If they represent local time, supply POSIXct/POSIXlt values ",
    "constructed with the correct time zone.\n",
    "If they are UTC, add Z: \"2020-01-02T19:54Z\"."
  ), class = "ocean3d_timezone_warning")
}
