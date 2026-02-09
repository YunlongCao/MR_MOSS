# MRMOSS

MRMOSS is packaged as an installable R project for multi-outcome Mendelian randomization.

This repository now supports a verified minimal end-to-end pipeline:

- start from raw summary statistics
- format traits to MRMOSS standard columns
- run MR-MOSS and write final result tables

## 1. Install

```r
# install.packages("devtools")
devtools::install_github("<your-org>/<your-repo>")
library(MRMOSS)
```

## 2. Configure Data Paths

Set data roots used by `inst/config/raw_data_manifest.csv`:

```r
Sys.setenv(
  MVP_SUMMARYDATA = "/path/to/MVP_summarydata",
  UKB_SUMMARYDATA = "/path/to/ukb_summarydata",
  GWAS_CATALOG_SUMMARYDATA = "/path/to/GWAS_catalog_summarydata",
  MRMOSS_INTERNAL_DATA = "/path/to/internal_data"
)
```

Set LD reference panel prefix (without `.bed/.bim/.fam`):

```r
reference_prefix <- "/path/to/reference/EUR"
```

If `plinkbinr` is unavailable, provide PLINK path with `MRMOSS_PLINK_BIN`.

## 3. Quickstart (Verified Raw -> Result)

The fastest reproducible path is a minimal real-data run:

- exposure: `Ground_coffee_consumption`
- outcomes: `Emotional_neglect`, `Physical_abuse`, `Sexual_abuse`
- threshold: `5e-7`

### Option A: one script

From a repository clone:

```bash
Rscript inst/scripts/00_quickstart_raw_to_result.R
```

Useful env vars:

- `MRMOSS_PROJECT_ROOT` (default `.`)
- `MRMOSS_MANIFEST` (default packaged manifest)
- `MRMOSS_FORMATTED_DIR`
- `MRMOSS_OUTPUT_DIR`
- `MRMOSS_REFERENCE_PREFIX`
- `MRMOSS_PLINK_BIN`
- `MRMOSS_OVERWRITE_FORMAT` (`true/false`)
- `MRMOSS_INCLUDE_OTHER_METHODS` (`true/false`, default `false`)

### Option B: direct function calls

```r
manifest <- system.file("config/raw_data_manifest.csv", package = "MRMOSS")
traits <- c("Ground_coffee_consumption", "Emotional_neglect", "Physical_abuse", "Sexual_abuse")

mrmoss_batch_format_manifest(
  manifest_csv = manifest,
  output_dir = "data/formatted_summary_data_quickstart",
  traits = traits,
  overwrite = FALSE,
  verbose = TRUE
)

mrmoss_run_analysis(
  exposures = "Ground_coffee_consumption",
  outcomes = c("Emotional_neglect", "Physical_abuse", "Sexual_abuse"),
  formatted_dir = "data/formatted_summary_data_quickstart",
  reference_prefix = "/path/to/reference/EUR",
  output_dir = "results/quickstart_minimal",
  output_prefix = "quickstart_groundcoffee_negctrl",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)
```

## 4. Full Manuscript Profiles (Optional)

Profile wrappers are still available:

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

By default, script wrappers set `include_other_methods = FALSE` to avoid extra package requirements.
Enable comparisons (IVW/RAPS/Egger/MRMix) with:

```bash
export MRMOSS_INCLUDE_OTHER_METHODS=true
```

## 5. Outputs

Each run writes:

- `<prefix>.txt`: MR-MOSS estimates and p-values
- `<prefix>_othermethods.txt`: optional comparison methods (only when enabled)
- `<prefix>_R_matrix.tsv`: estimated outcome correlation matrix

Quickstart also writes `format_summary.tsv` in the output directory.

## 6. Notes

- Core C++ implementation: `src/MRMOSS_PX.cpp`
- Manifest: `inst/config/raw_data_manifest.csv`
- Raw-to-formatted scripts: `inst/scripts/`

If installation fails in a conda environment due to C/C++ wrapper mismatch, set `R_MAKEVARS_USER` to use system `gcc/g++` and reinstall.
