#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop(
    "Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(MRMOSS))

as_flag <- function(x, default = FALSE) {
  if (!nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

demo_root <- system.file("extdata/quickstart_demo", package = "MRMOSS")
if (!nzchar(demo_root)) {
  stop("Cannot find packaged quickstart demo data.", call. = FALSE)
}

manifest <- file.path(demo_root, "manifest_quickstart_demo.csv")
reference_prefix <- file.path(demo_root, "reference", "EUR")
formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = file.path(getwd(), "results", "quickstart_demo", "formatted"))
output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "quickstart_demo"))
plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}
include_other_methods <- as_flag(Sys.getenv("MRMOSS_INCLUDE_OTHER_METHODS", unset = ""), default = FALSE)

dir.create(formatted_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Demo root:", demo_root, "\n")
cat("Manifest:", manifest, "\n")
cat("Reference prefix:", reference_prefix, "\n")
cat("Formatted dir:", formatted_dir, "\n")
cat("Output dir:", output_dir, "\n")
cat("Include other methods:", include_other_methods, "\n")

cat("\nStep 1/2: format demo raw files...\n")
format_summary <- mrmoss_batch_format_manifest(
  manifest_csv = manifest,
  output_dir = formatted_dir,
  traits = c("Ground_coffee_consumption", "Emotional_neglect", "Physical_abuse", "Sexual_abuse"),
  root_dir = demo_root,
  overwrite = TRUE,
  verbose = TRUE
)
data.table::fwrite(format_summary, file.path(output_dir, "format_summary.tsv"), sep = "\t")

cat("\nStep 2/2: run MR-MOSS on demo data...\n")
res <- mrmoss_run_analysis(
  exposures = "Ground_coffee_consumption",
  outcomes = c("Emotional_neglect", "Physical_abuse", "Sexual_abuse"),
  formatted_dir = formatted_dir,
  reference_prefix = reference_prefix,
  output_dir = output_dir,
  output_prefix = "quickstart_demo_result",
  iv_thresholds = c(5e-07),
  rd = 1.2,
  plink_bin = plink_bin,
  include_other_methods = include_other_methods,
  cache_dir = file.path(output_dir, "MRdat_cache"),
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(res$files)
