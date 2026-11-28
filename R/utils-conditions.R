# R/utils-conditions.R
#
# Custom conditions shared by several functions.
#
# Each helper signals the error or warning for one state of the input.
# Keeping them here gives them the same wording. Condition classes are
# internal: they let tests catch a specific condition without matching its
# wording, and are not documented for users. Nothing here is exported.


# Warning, class `ocean3d_timezone_warning`: some text date-times have no time
# zone (no trailing Z) and are being read as UTC. `prepared` is the result of
# .prepare_time_input(). Raised once per call.
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
