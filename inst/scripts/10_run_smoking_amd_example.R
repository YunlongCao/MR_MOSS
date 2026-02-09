#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop("Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').", call. = FALSE)
}
suppressPackageStartupMessages(library(MRMOSS))

formatted_data <- mrmoss_example_formatted_data()
exposure_name <- "Smoking_initiation"
outcome_names <- c("AMD", "AMD_dry", "AMD_wet")

reference_prefix <- Sys.getenv("MRMOSS_REFERENCE_PREFIX", unset = "")
if (!nzchar(reference_prefix) || toupper(reference_prefix) == "ONLINE") {
  reference_prefix <- NULL
}

plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}

pop <- Sys.getenv("MRMOSS_POP", unset = "EUR")
output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "smoking_amd_example"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Exposure:", exposure_name, "\n")
cat("Outcomes:", paste(outcome_names, collapse = ", "), "\n")
cat("Output dir:", output_dir, "\n")
if (is.null(reference_prefix)) {
  cat("Clumping mode: online OpenGWAS (set OPENGWAS_JWT)\n")
} else {
  cat("Clumping mode: local EUR reference:", reference_prefix, "\n")
}

prepared_input <- mrmoss_input(
  exposure = exposure_name,
  outcomes = outcome_names,
  formatted_data = formatted_data,
  iv_threshold = 5e-7,
  reference_prefix = reference_prefix,
  plink_bin = plink_bin,
  pop = pop,
  verbose = TRUE
)

fit_result <- mrmoss(
  exposure = exposure_name,
  outcomes = outcome_names,
  mrmoss_input = prepared_input,
  pvalue_output = "both",
  output_dir = output_dir,
  output_prefix = "smoking_to_amd",
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(fit_result$files)
cat("\nResult table:\n")
print(fit_result$result)
