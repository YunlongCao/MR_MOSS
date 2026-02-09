# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with summary statistics.

## Install

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

If `install_github()` hits GitHub API limit (`HTTP 403 rate limit exceeded`), use:

```bash
git clone https://github.com/YunlongCao/MR_MOSS.git
```

```r
remotes::install_local("MR_MOSS")
library(MRMOSS)
```

## MRAPSS-Like Quick Run (Built-In Data)

MRMOSS now provides built-in formatted example GWAS data:

- exposure: `Smoking_initiation`
- outcomes: `AMD`, `AMD_dry`, `AMD_wet`

Quick run in one command:

```r
MRres <- mrmoss_run_example()
MRres$files
MRres$mrmoss
```

Or explicit style (similar spirit to MRAPSS `data(...)` workflow):

```r
example_data <- mrmoss_example_formatted_data()

MRres <- mrmoss_run_analysis(
  exposures = "Smoking_initiation",
  outcomes = c("AMD", "AMD_dry", "AMD_wet"),
  formatted_data = example_data,
  reference_prefix = mrmoss_example_reference(),
  output_dir = "results/smoking_amd_example",
  output_prefix = "smoking_to_amd",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)

MRres$files
MRres$mrmoss
```

`mrmoss_example_reference()` returns a bundled lightweight local reference used for demonstration.

## Formatted Input Requirements

For your own data, each trait must contain columns:

- `SNP`, `A1`, `A2`, `Z`, `N`, `P`

You can provide inputs in two ways:

1. `formatted_dir`: one trait file per trait name (`file.path(formatted_dir, trait_name)`)
2. `formatted_data`: a named in-memory list of data.frames/data.tables

Quick file check:

```r
check_tab <- mrmoss_check_formatted_inputs(
  formatted_dir = "/path/to/formatted_data",
  traits = c("MyExposure", "Outcome1", "Outcome2")
)
check_tab
```

## Run MR-MOSS on Your Own Data

### Option 1: local PLINK clumping reference (recommended)

```r
res <- mrmoss_run_analysis(
  exposures = "MyExposure",
  outcomes = c("Outcome1", "Outcome2", "Outcome3"),
  formatted_dir = "/path/to/formatted_data",
  reference_prefix = "/path/to/EUR",  # expects EUR.bed/.bim/.fam
  output_dir = "results/my_run",
  output_prefix = "my_mrmoss_result",
  iv_thresholds = c(5e-7),
  rd = 1.2,
  include_other_methods = FALSE,
  verbose = TRUE
)
```

### Option 2: online clumping via OpenGWAS (`ieugwasr`)

```r
Sys.setenv(OPENGWAS_JWT = "<your_token>")

res <- mrmoss_run_analysis(
  exposures = "MyExposure",
  outcomes = c("Outcome1", "Outcome2", "Outcome3"),
  formatted_dir = "/path/to/formatted_data",
  reference_prefix = NULL,
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

Built-in smoking->AMD example runner:

```bash
Rscript inst/scripts/10_run_smoking_amd_example.R
```

`inst/scripts/10_run_smoking_amd_example.R` uses bundled local reference by default.
Set `MRMOSS_REFERENCE_PREFIX=ONLINE` to use online clumping.

## Manuscript Reproduction Scripts

These scripts are kept for manuscript-specific workflows:

```bash
Rscript inst/scripts/01_format_raw_data.R
Rscript inst/scripts/02_run_negative_control.R
Rscript inst/scripts/03_run_positive_control.R
Rscript inst/scripts/04_run_amd_application.R
Rscript inst/scripts/05_run_all_realdata.R
```
