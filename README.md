# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with summary statistics.

## 1. Install

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

## 2. Fastest Reproducible Run (Recommended)

This repository now includes a **self-contained demo quickstart** with:

- 4 small raw GWAS-like files (`Ground_coffee_consumption`, `Emotional_neglect`, `Physical_abuse`, `Sexual_abuse`)
- a tiny PLINK reference panel (`EUR.bed/.bim/.fam`)

So a new user can run end-to-end without preparing external data.

Run from a local clone:

```bash
Rscript inst/scripts/00_quickstart_demo.R
```

Default outputs:

- `results/quickstart_demo/quickstart_demo_result.txt`
- `results/quickstart_demo/quickstart_demo_result_R_matrix.tsv`
- `results/quickstart_demo/format_summary.tsv`

## 3. Download Demo Assets via `wget` (optional)

If you only want the quickstart assets (without cloning the full repository):

```bash
wget -O quickstart_demo_assets.zip \
  https://raw.githubusercontent.com/YunlongCao/MR_MOSS/main/inst/extdata/quickstart_demo_assets.zip

unzip -o quickstart_demo_assets.zip
```

Or use the helper script:

```bash
bash inst/scripts/download_quickstart_demo_assets.sh
```

## 4. Full Raw-Data Quickstart (manuscript-scale)

Script:

```bash
Rscript inst/scripts/00_quickstart_raw_to_result.R
```

This path uses the manuscript-style manifest and requires external raw data.
Set at least:

- `MRMOSS_INTERNAL_DATA` (directory containing required raw files)
- `MRMOSS_REFERENCE_PREFIX` (LD reference prefix, without `.bed/.bim/.fam`)

Optional:

- `MRMOSS_PLINK_BIN`
- `MRMOSS_FORMATTED_DIR`, `MRMOSS_OUTPUT_DIR`
- `MRMOSS_INCLUDE_OTHER_METHODS=true`

## 5. Full Profile Scripts (optional)

For full real-data profiles, also set data roots used in `inst/config/raw_data_manifest.csv`:

- `MVP_SUMMARYDATA`
- `UKB_SUMMARYDATA`
- `GWAS_CATALOG_SUMMARYDATA`
- `MRMOSS_INTERNAL_DATA`

Then run:

```bash
Rscript inst/scripts/01_format_raw_data.R
Rscript inst/scripts/02_run_negative_control.R
Rscript inst/scripts/03_run_positive_control.R
Rscript inst/scripts/04_run_amd_application.R
Rscript inst/scripts/05_run_all_realdata.R
```

Generic runner:

```bash
Rscript inst/scripts/run_profile.R amd_application
```

## 6. Troubleshooting

- If quickstart script says package is missing: install `MRMOSS` first.
- If clumping fails: set `MRMOSS_PLINK_BIN` to a valid PLINK executable path.
- If no MR-MOSS row is produced, common causes are:
  - too few shared IVs after clumping (needs at least 3)
  - too few null SNPs for outcome-correlation estimation (needs at least 10)
- In conda environments, if compiler wrapper is missing, set `R_MAKEVARS_USER` to use system `gcc/g++` and reinstall.

## 7. Core Files

- Main C++ implementation: `src/MRMOSS_PX.cpp`
- Full manifest: `inst/config/raw_data_manifest.csv`
- Demo manifest: `inst/extdata/quickstart_demo/manifest_quickstart_demo.csv`
- Demo quickstart script: `inst/scripts/00_quickstart_demo.R`
