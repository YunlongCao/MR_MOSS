#' Load Built-In Formatted Example Data
#'
#' Loads the packaged real-data example used in the README:
#'
#' - exposure: `Smoking_initiation`
#' - outcomes: `AMD`, `AMD_dry`, `AMD_wet`
#'
#' @return Named list of formatted trait data.tables.
#' @export
mrmoss_example_formatted_data <- function() {
  formatted_dir <- system.file("extdata", "formatted_example", package = "MRMOSS")
  if (!nzchar(formatted_dir)) {
    stop("Cannot find built-in formatted example data in this MRMOSS installation.", call. = FALSE)
  }

  traits <- c("Smoking_initiation", "AMD", "AMD_dry", "AMD_wet")
  mrmoss_load_formatted_traits(formatted_dir, traits)
}
