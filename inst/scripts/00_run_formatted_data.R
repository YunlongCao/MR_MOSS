#!/usr/bin/env Rscript

if (!requireNamespace("MRMOSS", quietly = TRUE)) {
  stop("Package 'MRMOSS' is not installed. Install with remotes::install_github('YunlongCao/MR_MOSS').", call. = FALSE)
}
suppressPackageStartupMessages(library(MRMOSS))

parse_flag <- function(x, default = FALSE) {
  if (!nzchar(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_list <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

parse_numeric_list <- function(x, default) {
  if (!nzchar(x)) {
    return(default)
  }
  vals <- suppressWarnings(as.numeric(parse_list(x)))
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) {
    stop("MRMOSS_IV_THRESHOLDS must contain positive numeric values, e.g. '5e-7,5e-8'.", call. = FALSE)
  }
  vals
}

formatted_dir <- Sys.getenv("MRMOSS_FORMATTED_DIR", unset = "")
if (!nzchar(formatted_dir)) {
  stop("Set MRMOSS_FORMATTED_DIR to your formatted-data directory.", call. = FALSE)
}

exposures <- parse_list(Sys.getenv("MRMOSS_EXPOSURES", unset = ""))
outcomes <- parse_list(Sys.getenv("MRMOSS_OUTCOMES", unset = ""))
if (length(exposures) == 0 || length(outcomes) == 0) {
  stop("Set both MRMOSS_EXPOSURES and MRMOSS_OUTCOMES (comma-separated trait names).", call. = FALSE)
}

reference_prefix <- Sys.getenv("MRMOSS_REFERENCE_PREFIX", unset = "")
if (!nzchar(reference_prefix)) {
  reference_prefix <- NULL
}

output_dir <- Sys.getenv("MRMOSS_OUTPUT_DIR", unset = file.path(getwd(), "results", "mrmoss_formatted_run"))
output_prefix <- Sys.getenv("MRMOSS_OUTPUT_PREFIX", unset = "mrmoss_result")
plink_bin <- Sys.getenv("MRMOSS_PLINK_BIN", unset = "")
if (!nzchar(plink_bin)) {
  plink_bin <- NULL
}

iv_thresholds <- parse_numeric_list(
  Sys.getenv("MRMOSS_IV_THRESHOLDS", unset = ""),
  default = c(5e-7, 5e-8)
)
rd <- suppressWarnings(as.numeric(Sys.getenv("MRMOSS_RD", unset = "1.2")))
if (!is.finite(rd) || rd <= 0) {
  stop("MRMOSS_RD must be a positive numeric value.", call. = FALSE)
}
include_other_methods <- parse_flag(Sys.getenv("MRMOSS_INCLUDE_OTHER_METHODS", unset = ""), default = FALSE)
pop <- Sys.getenv("MRMOSS_POP", unset = "EUR")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Formatted dir:", formatted_dir, "\n")
cat("Exposures:", paste(exposures, collapse = ", "), "\n")
cat("Outcomes:", paste(outcomes, collapse = ", "), "\n")
cat("Output dir:", output_dir, "\n")
cat("Output prefix:", output_prefix, "\n")
if (is.null(reference_prefix)) {
  cat("Clumping mode: online OpenGWAS (ieugwasr)\n")
  if (!nzchar(Sys.getenv("OPENGWAS_JWT", unset = ""))) {
    cat("Note: OPENGWAS_JWT is not set; OpenGWAS may return 401.\n")
  }
} else {
  cat("Clumping mode: local PLINK reference:", reference_prefix, "\n")
}

check_tab <- mrmoss_check_formatted_inputs(
  formatted_dir = formatted_dir,
  traits = unique(c(exposures, outcomes))
)
print(check_tab)
if (any(check_tab$status != "ok")) {
  stop("Formatted input checks failed. Fix files before running MR-MOSS.", call. = FALSE)
}

res <- mrmoss_run_analysis(
  exposures = exposures,
  outcomes = outcomes,
  formatted_dir = formatted_dir,
  reference_prefix = reference_prefix,
  output_dir = output_dir,
  output_prefix = output_prefix,
  iv_thresholds = iv_thresholds,
  rd = rd,
  plink_bin = plink_bin,
  pop = pop,
  include_other_methods = include_other_methods,
  cache_dir = file.path(output_dir, "MRdat_cache"),
  verbose = TRUE
)

cat("\nDone. Output files:\n")
print(res$files)
