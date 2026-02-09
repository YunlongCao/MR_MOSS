#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MRMOSS))

as_flag <- function(x, default = FALSE) {
  if (!nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

project_root <- normalizePath(Sys.getenv("MRMOSS_PROJECT_ROOT", "."), winslash = "/", mustWork = FALSE)

manifest <- Sys.getenv("MRMOSS_MANIFEST", unset = "")
if (!nzchar(manifest)) {
  manifest <- system.file("config/raw_data_manifest.csv", package = "MRMOSS")
}
if (!nzchar(manifest)) {
  manifest <- file.path(project_root, "inst", "config", "raw_data_manifest.csv")
}

formatted_dir <- Sys.getenv(
  "MRMOSS_FORMATTED_DIR",
  unset = file.path(project_root, "data", "formatted_summary_data_quickstart")
)
output_dir <- Sys.getenv(
  "MRMOSS_OUTPUT_DIR",
  unset = file.path(project_root, "results", "quickstart_minimal")
)
reference_prefix <- Sys.getenv(
  "MRMOSS_REFERENCE_PREFIX",
  unset = file.path(project_root, "data", "reference", "EUR")
)
plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}

overwrite_format <- as_flag(Sys.getenv("MRMOSS_OVERWRITE_FORMAT", unset = ""), default = FALSE)
include_other_methods <- as_flag(Sys.getenv("MRMOSS_INCLUDE_OTHER_METHODS", unset = ""), default = FALSE)

traits <- c(
  "Ground_coffee_consumption",
  "Emotional_neglect",
  "Physical_abuse",
  "Sexual_abuse"
)
outcomes <- c("Emotional_neglect", "Physical_abuse", "Sexual_abuse")

cat("Manifest:", manifest, "\n")
cat("Formatted dir:", formatted_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Reference prefix:", reference_prefix, "\n")
cat("Include other methods:", include_other_methods, "\n")

dir.create(formatted_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nStep 1/2: format raw inputs...\n")
format_summary <- mrmoss_batch_format_manifest(
  manifest_csv = manifest,
  output_dir = formatted_dir,
  traits = traits,
  root_dir = project_root,
  overwrite = overwrite_format,
  verbose = TRUE
)
format_summary_file <- file.path(output_dir, "format_summary.tsv")
data.table::fwrite(format_summary, format_summary_file, sep = "\t")

cat("\nStep 2/2: run MR-MOSS...\n")
res <- mrmoss_run_analysis(
  exposures = "Ground_coffee_consumption",
  outcomes = outcomes,
  formatted_dir = formatted_dir,
  reference_prefix = reference_prefix,
  output_dir = output_dir,
  output_prefix = "quickstart_groundcoffee_negctrl",
  iv_thresholds = c(5e-07),
  rd = 1.2,
  plink_bin = plink_bin,
  include_other_methods = include_other_methods,
  cache_dir = file.path(output_dir, "MRdat_cache"),
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(res$files)
