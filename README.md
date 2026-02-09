# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with summary statistics.

This repository is organized so users can run:

- raw summary statistics -> formatted files
- formatted files -> final MR-MOSS result table

## 1. Prerequisites

- R >= 4.1
- A working C/C++ toolchain for compiling `src/`
- PLINK binary (or `plinkbinr`)

If `Rscript` is not in PATH, activate your R environment first (for example your conda env).

## 2. Install

### Option A (recommended): install from GitHub

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

### Option B: install from a local clone

```bash
git clone git@github.com:YunlongCao/MR_MOSS.git
cd MR_MOSS
R CMD INSTALL .
```

## 3. Quickstart (minimum end-to-end run)

This is the easiest path to verify the pipeline from raw input to final MR-MOSS output.

Quickstart uses:

- exposure: `Ground_coffee_consumption`
- outcomes: `Emotional_neglect`, `Physical_abuse`, `Sexual_abuse`
- IV threshold: `5e-7`

### Required environment variables for quickstart

Only two paths are required:

- `MRMOSS_INTERNAL_DATA`: directory containing quickstart raw files
- `MRMOSS_REFERENCE_PREFIX`: LD reference prefix (without `.bed/.bim/.fam`)

Optional:

- `MRMOSS_PLINK_BIN`: explicit path to PLINK executable
- `MRMOSS_FORMATTED_DIR`, `MRMOSS_OUTPUT_DIR`: output locations

### One-command quickstart (from repository root)

```bash
export MRMOSS_INTERNAL_DATA=/path/to/internal_data
export MRMOSS_REFERENCE_PREFIX=/path/to/reference/EUR
# optional: export MRMOSS_PLINK_BIN=/path/to/plink

Rscript inst/scripts/00_quickstart_raw_to_result.R
```

Expected outputs (default under `results/quickstart_minimal`):

- `quickstart_groundcoffee_negctrl.txt`
- `quickstart_groundcoffee_negctrl_R_matrix.tsv`
- `format_summary.tsv`

## 4. Full manuscript profiles (optional)

These require additional data roots referenced in `inst/config/raw_data_manifest.csv`:

- `MVP_SUMMARYDATA`
- `UKB_SUMMARYDATA`
- `GWAS_CATALOG_SUMMARYDATA`
- `MRMOSS_INTERNAL_DATA`

Run profile scripts:

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

By default, wrapper scripts set `include_other_methods = FALSE`.
Enable IVW/RAPS/Egger/MRMix with:

```bash
export MRMOSS_INCLUDE_OTHER_METHODS=true
```

## 5. Core files

- Main C++ implementation: `src/MRMOSS_PX.cpp`
- Data manifest: `inst/config/raw_data_manifest.csv`
- Quickstart script: `inst/scripts/00_quickstart_raw_to_result.R`
- Main analysis API: `mrmoss_run_analysis()`

## 6. Troubleshooting

If compilation fails in conda environments (for example `x86_64-conda-linux-gnu-c++: not found`), set `R_MAKEVARS_USER` to point to a Makevars file that uses system `gcc/g++`, then reinstall.
