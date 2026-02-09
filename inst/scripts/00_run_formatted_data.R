#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop("Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').", call. = FALSE)
}
suppressPackageStartupMessages(library(MRMOSS))

parse_list <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = "")
if (!nzchar(formatted_dir)) {
  stop("Set MRMOSS_FORMATTED_DIR to your formatted-data directory.", call. = FALSE)
}

exposure_name <- Sys.getenv("MRMOSS_EXPOSURE", unset = "")
outcome_names <- parse_list(Sys.getenv("MRMOSS_OUTCOMES", unset = ""))
if (!nzchar(exposure_name) || length(outcome_names) == 0) {
  stop("Set MRMOSS_EXPOSURE and MRMOSS_OUTCOMES (comma-separated).", call. = FALSE)
}

reference_prefix <- Sys.getenv("MRMOSS_REFERENCE_PREFIX", unset = "")
if (!nzchar(reference_prefix) || toupper(reference_prefix) == "ONLINE") {
  reference_prefix <- NULL
}

plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}

pop <- Sys.getenv("MRMOSS_POP", unset = "EUR")
output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "mrmoss_formatted_run"))
output_prefix <- Sys.getenv("MRMOSS_OUTPUT_PREFIX", unset = "mrmoss_result")
pvalue_output <- Sys.getenv("MRMOSS_PVALUE_OUTPUT", unset = "both")
iv_threshold <- suppressWarnings(as.numeric(Sys.getenv("MRMOSS_IV_THRESHOLD", unset = "5e-7")))
if (!is.finite(iv_threshold) || iv_threshold <= 0) {
  stop("MRMOSS_IV_THRESHOLD must be a positive numeric value.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Formatted dir:", formatted_dir, "\n")
cat("Exposure:", exposure_name, "\n")
cat("Outcomes:", paste(outcome_names, collapse = ", "), "\n")
cat("Output dir:", output_dir, "\n")

check_tab <- mrmoss_check_formatted_inputs(
  formatted_dir = formatted_dir,
  traits = unique(c(exposure_name, outcome_names))
)
print(check_tab)
if (any(check_tab$status != "ok")) {
  stop("Formatted input checks failed. Fix files before running MRMOSS.", call. = FALSE)
}

prepared_input <- mrmoss_input(
  exposure = exposure_name,
  outcomes = outcome_names,
  formatted_dir = formatted_dir,
  iv_threshold = iv_threshold,
  reference_prefix = reference_prefix,
  plink_bin = plink_bin,
  pop = pop,
  verbose = TRUE
)

fit_result <- mrmoss(
  exposure = exposure_name,
  outcomes = outcome_names,
  mrmoss_input = prepared_input,
  pvalue_output = pvalue_output,
  output_dir = output_dir,
  output_prefix = output_prefix,
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(fit_result$files)
cat("\nResult table:\n")
print(fit_result$result)
