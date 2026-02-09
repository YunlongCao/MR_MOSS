#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop("Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').", call. = FALSE)
}
suppressPackageStartupMessages(library(MRMOSS))

example_data <- mrmoss_example_formatted_data()

exposures <- "Smoking_initiation"
outcomes <- c("AMD", "AMD_dry", "AMD_wet")

reference_prefix <- Sys.getenv("MRMOSS_REFERENCE_PREFIX", unset = "")
if (!nzchar(reference_prefix)) {
  reference_prefix <- mrmoss_example_reference()
}
if (toupper(reference_prefix) == "ONLINE") {
  reference_prefix <- NULL
}
plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}
pop <- Sys.getenv("MRMOSS_POP", unset = "EUR")

output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "smoking_amd_example"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Exposure:", exposures, "\n")
cat("Outcomes:", paste(outcomes, collapse = ", "), "\n")
cat("Output dir:", output_dir, "\n")
if (is.null(reference_prefix)) {
  cat("Clumping mode: online OpenGWAS (ieugwasr)\n")
  if (!nzchar(Sys.getenv("OPENGWAS_JWT", unset = ""))) {
    cat("Note: OPENGWAS_JWT is not set; OpenGWAS may return 401.\n")
  }
} else {
  cat("Clumping mode: local PLINK reference:", reference_prefix, "\n")
  if (!file.exists(paste0(reference_prefix, ".bed"))) {
    stop(sprintf("Reference file not found: %s.bed", reference_prefix), call. = FALSE)
  }
}

res <- mrmoss_run_analysis(
  exposures = exposures,
  outcomes = outcomes,
  formatted_data = example_data,
  reference_prefix = reference_prefix,
  output_dir = output_dir,
  output_prefix = "smoking_to_amd",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  plink_bin = plink_bin,
  pop = pop,
  include_other_methods = FALSE,
  cache_dir = file.path(output_dir, "MRdat_cache"),
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(res$files)
