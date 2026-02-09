# Manifest Columns (`raw_data_manifest.csv`)

- `trait`: output trait name (formatted file name)
- `raw_path`: raw input path (supports `${ENV_VAR}` placeholders)
- `source_type`:
  - `beta_se`: derive `Z = beta / se`
  - `formatted`: file already has `SNP/A1/A2/Z/N/P`
- `snp_col`, `a1_col`, `a2_col`, `beta_col`, `se_col`, `z_col`, `p_col`, `n_col`: source column names
- `n_value`: fixed sample size when `n_col` is absent
- `note`: free-text annotation

Environment-variable placeholders are expanded by `mrmoss_expand_env_vars()`.
