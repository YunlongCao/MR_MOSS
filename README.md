# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with summary statistics.

This repository includes a real-data quickstart subset from the manuscript so users can reproduce one complete MR-MOSS analysis:

- exposure: `Ground_coffee_consumption`
- outcomes: `Emotional_neglect`, `Physical_abuse`, `Sexual_abuse`

## Install

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

## Quickstart (Real Data, Function-Based)

### 1. Prepare bundled quickstart assets

The package ships split zip parts (`quickstart_real_assets.zip.part-aa/ab`) containing:

- 4 quickstart GWAS files (1 exposure + 3 outcomes)
- matching PLINK reference files (`EUR.bed/.bim/.fam`)

Use R to reconstruct and unzip them:

```r
library(MRMOSS)

qs <- mrmoss_prepare_quickstart_real_data(
  output_dir = file.path(getwd(), "quickstart_real_assets"),
  overwrite = FALSE,
  verbose = TRUE
)

qs$quickstart_root
qs$manifest_csv
qs$reference_prefix
```

### 2. Format quickstart raw files

```r
formatted_dir <- file.path(getwd(), "results", "quickstart_real", "formatted")
dir.create(formatted_dir, recursive = TRUE, showWarnings = FALSE)

format_summary <- mrmoss_batch_format_manifest(
  manifest_csv = qs$manifest_csv,
  output_dir = formatted_dir,
  traits = c(qs$exposure, qs$outcomes),
  root_dir = qs$quickstart_root,
  overwrite = TRUE,
  verbose = TRUE
)

format_summary
```

### 3. Run MR-MOSS for 1 exposure -> 3 outcomes

```r
output_dir <- file.path(getwd(), "results", "quickstart_real")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

res <- mrmoss_run_analysis(
  exposures = qs$exposure,
  outcomes = qs$outcomes,
  formatted_dir = formatted_dir,
  reference_prefix = qs$reference_prefix,
  output_dir = output_dir,
  output_prefix = "quickstart_real_result",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)

res$files
read.delim(res$files$mrmoss)
```

Expected files:

- `results/quickstart_real/quickstart_real_result.txt`
- `results/quickstart_real/quickstart_real_result_R_matrix.tsv`

## Optional: Download Zip Parts via `wget`

```bash
wget -O quickstart_real_assets.zip.part-aa \
  https://raw.githubusercontent.com/YunlongCao/MR_MOSS/main/inst/extdata/quickstart_real_assets.zip.part-aa
wget -O quickstart_real_assets.zip.part-ab \
  https://raw.githubusercontent.com/YunlongCao/MR_MOSS/main/inst/extdata/quickstart_real_assets.zip.part-ab

cat quickstart_real_assets.zip.part-aa quickstart_real_assets.zip.part-ab > quickstart_real_assets.zip
unzip -o quickstart_real_assets.zip
```

Or run:

```bash
bash inst/scripts/download_quickstart_real_assets.sh
```

## Optional Script Runner

```bash
Rscript inst/scripts/00_quickstart_real.R
```

## Full Manuscript-Scale Pipelines (Optional)

```bash
Rscript inst/scripts/00_quickstart_raw_to_result.R
Rscript inst/scripts/01_format_raw_data.R
Rscript inst/scripts/02_run_negative_control.R
Rscript inst/scripts/03_run_positive_control.R
Rscript inst/scripts/04_run_amd_application.R
Rscript inst/scripts/05_run_all_realdata.R
```

These require external data paths (for example `MRMOSS_INTERNAL_DATA`, `MRMOSS_REFERENCE_PREFIX`).

## Troubleshooting

- If clumping fails, set `plink_bin` in `mrmoss_run_analysis()` or install `plinkbinr`.
- If no MR-MOSS row is produced, common causes are:
  - fewer than 3 shared IVs after clumping
  - fewer than 10 null SNPs for outcome-correlation estimation

## Core APIs

- `mrmoss_prepare_quickstart_real_data()`
- `mrmoss_batch_format_manifest()`
- `mrmoss_run_analysis()`
