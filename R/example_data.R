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

#' Get Built-In Example Reference Prefix
#'
#' Returns the prefix (without extension) for the packaged lightweight
#' reference panel used by the built-in example.
#'
#' @return Character scalar path to bundled example reference prefix.
#' @export
mrmoss_example_reference <- function() {
  formatted_dir <- system.file("extdata", "formatted_example", package = "MRMOSS")
  if (!nzchar(formatted_dir)) {
    stop("Cannot find built-in formatted example data in this MRMOSS installation.", call. = FALSE)
  }

  ref <- file.path(formatted_dir, "EUR_example")
  if (!file.exists(paste0(ref, ".bed"))) {
    stop(sprintf("Bundled example reference is missing: %s.bed", ref), call. = FALSE)
  }
  ref
}

#' Run Built-In Smoking->AMD Example
#'
#' Runs MR-MOSS on the built-in formatted example data with one exposure
#' (`Smoking_initiation`) and three outcomes (`AMD`, `AMD_dry`, `AMD_wet`).
#'
#' @param output_dir Output directory.
#' @param output_prefix Output basename (without extension).
#' @param use_online_clump Whether to use online OpenGWAS clumping instead of
#'   bundled local reference.
#' @param ... Extra arguments passed to `mrmoss_run_analysis`.
#' @return Result from `mrmoss_run_analysis`.
#' @export
mrmoss_run_example <- function(output_dir = file.path(getwd(), "results", "smoking_amd_example"),
                               output_prefix = "smoking_to_amd",
                               use_online_clump = FALSE,
                               ...) {
  dat <- mrmoss_example_formatted_data()
  ref <- if (isTRUE(use_online_clump)) NULL else mrmoss_example_reference()

  mrmoss_run_analysis(
    exposures = "Smoking_initiation",
    outcomes = c("AMD", "AMD_dry", "AMD_wet"),
    formatted_data = dat,
    reference_prefix = ref,
    output_dir = output_dir,
    output_prefix = output_prefix,
    iv_thresholds = c(5e-7),
    rd = 1.2,
    include_other_methods = FALSE,
    ...
  )
}
