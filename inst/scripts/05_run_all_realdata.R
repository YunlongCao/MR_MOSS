#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MRMOSS))

project_root <- normalizePath(Sys.getenv("MRMOSS_PROJECT_ROOT", "."), winslash = "/", mustWork = FALSE)
manifest <- Sys.getenv("MRMOSS_MANIFEST", unset = system.file("config/raw_data_manifest.csv", package = "MRMOSS"))
if (!nzchar(manifest)) {
  manifest <- file.path(project_root, "inst", "config", "raw_data_manifest.csv")
}

formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = file.path(project_root, "data", "formatted_summary_data"))
reference_prefix <- Sys.getenv("MRMOSS_REFERENCE_PREFIX", unset = file.path(project_root, "data", "reference", "EUR"))
plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}
include_other_methods <- tolower(Sys.getenv("MRMOSS_INCLUDE_OTHER_METHODS", unset = "false")) %in%
  c("1", "true", "t", "yes", "y")

profiles <- c("negative_control", "positive_control", "amd_application")

for (p in profiles) {
  cat("\n=== Running profile:", p, "===\n")
  out <- file.path(project_root, "results", p)
  res <- mrmoss_run_profile(
    profile = p,
    formatted_dir = formatted_dir,
    reference_prefix = reference_prefix,
    output_dir = out,
    manifest_csv = manifest,
    root_dir = project_root,
    overwrite_format = FALSE,
    plink_bin = plink_bin,
    include_other_methods = include_other_methods,
    cache_dir = file.path(out, "MRdat_cache")
  )
  print(res$files)
}
