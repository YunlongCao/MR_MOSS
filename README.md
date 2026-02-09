# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with summary statistics.

This repository now focuses on a **formatted data -> MR-MOSS result** workflow.

## Install

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

## Formatted Input Requirements

`mrmoss_run_analysis()` expects one file per trait under `formatted_dir`.

- File path rule: `file.path(formatted_dir, trait_name)`
- File format: tab-delimited text
- Required columns: `SNP`, `A1`, `A2`, `Z`, `N`, `P`
- Optional extra columns are allowed (for example `chi2`)

Column expectations:

- `SNP`: rsID-style variant identifier
- `A1`, `A2`: alleles (`A/C/G/T`)
- `Z`: SNP effect Z-score
- `N`: sample size (> 0)
- `P`: p-value in `[0, 1]`

Quick check before analysis:

```r
check_tab <- mrmoss_check_formatted_inputs(
  formatted_dir = "/path/to/formatted_data",
  traits = c("MyExposure", "Outcome1", "Outcome2")
)
check_tab
```

## Run MR-MOSS on Your Own Formatted Data

### Option 1 (recommended): local PLINK clumping reference

```r
formatted_dir <- "/path/to/formatted_data"
exposures <- "MyExposure"
outcomes <- c("Outcome1", "Outcome2", "Outcome3")

res <- mrmoss_run_analysis(
  exposures = exposures,
  outcomes = outcomes,
  formatted_dir = formatted_dir,
  reference_prefix = "/path/to/EUR",  # expects EUR.bed/.bim/.fam
  output_dir = "results/my_run",
  output_prefix = "my_mrmoss_result",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)

res$files
read.delim(res$files$mrmoss)
```

### Option 2: online clumping via `ieugwasr` (no local EUR files)

```r
# Sys.setenv(OPENGWAS_JWT = "<your_token>")  # usually required by OpenGWAS API

res <- mrmoss_run_analysis(
  exposures = "MyExposure",
  outcomes = c("Outcome1", "Outcome2", "Outcome3"),
  formatted_dir = "/path/to/formatted_data",
  reference_prefix = NULL,  # online clumping mode
  pop = "EUR",
  output_dir = "results/my_run_online",
  output_prefix = "my_mrmoss_result_online",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)
```

If OpenGWAS returns `401 Invalid token`, set `OPENGWAS_JWT` or switch to local clumping.

## Included Real-Data Example (Smoking -> AMD outcomes)

The repository includes these formatted files in `inst/extdata/formatted_example`:

- exposure: `Smoking_initiation`
- outcomes: `AMD`, `AMD_dry`, `AMD_wet`

Run directly in R:

```r
example_dir <- system.file("extdata", "formatted_example", package = "MRMOSS")

res <- mrmoss_run_analysis(
  exposures = "Smoking_initiation",
  outcomes = c("AMD", "AMD_dry", "AMD_wet"),
  formatted_dir = example_dir,
  reference_prefix = NULL,  # set to "/path/to/EUR" for local clumping
  pop = "EUR",
  output_dir = "results/smoking_amd_example",
  output_prefix = "smoking_to_amd",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)

res$files
```

## Script Entry Points

Generic formatted-data runner:

```bash
MRMOSS_FORMATTED_DIR=/path/to/formatted \
MRMOSS_EXPOSURES=MyExposure \
MRMOSS_OUTCOMES=Outcome1,Outcome2,Outcome3 \
MRMOSS_OUTPUT_DIR=results/my_run \
MRMOSS_OUTPUT_PREFIX=my_mrmoss_result \
Rscript inst/scripts/00_run_formatted_data.R
```

Real-data example runner:

```bash
Rscript inst/scripts/10_run_smoking_amd_example.R
```

## Reference Files (`EUR.bed/.bim/.fam`)

The full EUR reference panel is too large to store directly in this repository.
Two practical options are supported:

- Use your own local EUR reference (`reference_prefix = "/path/to/EUR"`)
- Use online clumping (`reference_prefix = NULL`)

If your group hosts EUR files on OneDrive (or similar), download them first and then set `reference_prefix`.
Template:

```bash
wget -O EUR.bed "<DIRECT_URL_TO_EUR.bed>"
wget -O EUR.bim "<DIRECT_URL_TO_EUR.bim>"
wget -O EUR.fam "<DIRECT_URL_TO_EUR.fam>"
```

## Manuscript Reproduction Scripts

Scripts kept for manuscript-specific workflows:

```bash
Rscript inst/scripts/01_format_raw_data.R
Rscript inst/scripts/02_run_negative_control.R
Rscript inst/scripts/03_run_positive_control.R
Rscript inst/scripts/04_run_amd_application.R
Rscript inst/scripts/05_run_all_realdata.R
```

These are not required for the standard formatted-data workflow.
