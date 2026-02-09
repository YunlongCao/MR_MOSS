#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop(
    "Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(MRMOSS))

asset_dir <- Sys.getenv("MRMOSS_QUICKSTART_ASSET_DIR", unset = file.path(getwd(), "quickstart_real_assets"))
formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = file.path(getwd(), "results", "quickstart_real", "formatted"))
output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "quickstart_real"))

quickstart <- mrmoss_prepare_quickstart_real_data(
  output_dir = asset_dir,
  overwrite = FALSE,
  verbose = TRUE
)
manifest <- quickstart$manifest_csv
reference_prefix <- quickstart$reference_prefix

dir.create(formatted_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Quickstart root:", quickstart$quickstart_root, "\n")
cat("Manifest:", manifest, "\n")
cat("Reference:", reference_prefix, "\n")
cat("Formatted dir:", formatted_dir, "\n")
cat("Output dir:", output_dir, "\n")

format_summary <- mrmoss_batch_format_manifest(
  manifest_csv = manifest,
  output_dir = formatted_dir,
  traits = c(quickstart$exposure, quickstart$outcomes),
  root_dir = quickstart$quickstart_root,
  overwrite = TRUE,
  verbose = TRUE
)
data.table::fwrite(format_summary, file.path(output_dir, "format_summary.tsv"), sep = "\t")

res <- mrmoss_run_analysis(
  exposures = quickstart$exposure,
  outcomes = quickstart$outcomes,
  formatted_dir = formatted_dir,
  reference_prefix = reference_prefix,
  output_dir = output_dir,
  output_prefix = "quickstart_real_result",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  cache_dir = file.path(output_dir, "MRdat_cache"),
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(res$files)
