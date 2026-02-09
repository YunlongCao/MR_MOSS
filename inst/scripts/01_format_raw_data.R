#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MRMOSS))

project_root <- normalizePath(Sys.getenv("MRMOSS_PROJECT_ROOT", "."), winslash = "/", mustWork = FALSE)
manifest <- Sys.getenv("MRMOSS_MANIFEST", unset = "")
if (!nzchar(manifest)) {
  manifest <- system.file("config/raw_data_manifest.csv", package = "MRMOSS")
}
if (!nzchar(manifest)) {
  manifest <- file.path(project_root, "inst", "config", "raw_data_manifest.csv")
}

formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = file.path(project_root, "data", "formatted_summary_data"))

cat("Formatting with manifest:", manifest, "\n")
cat("Output directory:", formatted_dir, "\n")

summary_tab <- mrmoss_batch_format_manifest(
  manifest_csv = manifest,
  output_dir = formatted_dir,
  root_dir = project_root,
  overwrite = FALSE,
  verbose = TRUE
)

out_summary <- file.path(formatted_dir, "_format_summary.tsv")
data.table::fwrite(summary_tab, out_summary, sep = "\t")
cat("\nFormatting summary written to:", out_summary, "\n")
