# MRMOSS

MRMOSS is an R package for multi-outcome Mendelian randomization with GWAS summary statistics.

## Install

```r
# install.packages("remotes")
remotes::install_github("YunlongCao/MR_MOSS")
library(MRMOSS)
```

`MRMOSS` contains C++ code. If installation fails at compile step, install a C/C++ build toolchain first.

## Built-in formatted real-data example

The package includes formatted GWAS files for:

- Exposure: `Smoking_initiation`
- Outcomes: `AMD`, `AMD_dry`, `AMD_wet`

## Formatted data requirements

Each trait file (or each in-memory data.frame/data.table) must contain:

- `SNP` (rsID)
- `A1` (effect allele)
- `A2` (other allele)
- `Z` (z score)
- `N` (sample size)
- `P` (p-value)

You can provide formatted data to MRMOSS in two ways:

- `formatted_dir`: one file per trait, file name = trait name
- `formatted_data`: named list in memory

## Usage (1 exposure -> 3 outcomes)

### 1) Load formatted data

```r
library(MRMOSS)

exposure_name <- "Smoking_initiation"
outcome_names <- c("AMD", "AMD_dry", "AMD_wet")

formatted_data <- mrmoss_example_formatted_data()
```

### 2) Build MRMOSS input from formatted data

Use a real EUR LD reference for local clumping (recommended):

```r
prepared_input <- mrmoss_input(
  exposure = exposure_name,
  outcomes = outcome_names,
  formatted_data = formatted_data,
  iv_threshold = 5e-7,
  reference_prefix = "/path/to/EUR",   # EUR.bed / EUR.bim / EUR.fam
  pop = "EUR"
)
```

### 3) Run MRMOSS

`pvalue_output` supports:

- `"outcome"`: only outcome-specific p-values
- `"global"`: only global p-value
- `"both"`: both global and outcome-specific p-values

```r
fit_result <- mrmoss(
  exposure = exposure_name,
  outcomes = outcome_names,
  mrmoss_input = prepared_input,
  pvalue_output = "both",
  output_dir = "results/smoking_to_amd",
  output_prefix = "smoking_to_amd"
)

fit_result$result
fit_result$files
```

## Online clumping alternative (no local EUR files)

If local EUR reference files are unavailable, use OpenGWAS online clumping:

```r
Sys.setenv(OPENGWAS_JWT = "<your_token>")

prepared_input <- mrmoss_input(
  exposure = exposure_name,
  outcomes = outcome_names,
  formatted_data = formatted_data,
  iv_threshold = 5e-7,
  reference_prefix = NULL,
  pop = "EUR"
)
```

If online clumping returns `401`, check `OPENGWAS_JWT`.
