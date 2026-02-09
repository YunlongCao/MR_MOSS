#' Read Manifest for Raw-to-Formatted Conversion
#'
#' @param manifest_csv Path to CSV manifest.
#' @return A data.table manifest.
mrmoss_read_manifest <- function(manifest_csv) {
  manifest <- data.table::fread(manifest_csv, na.strings = c("", "NA", "NaN"))

  required <- c("trait", "raw_path", "source_type")
  missing_cols <- setdiff(required, colnames(manifest))
  if (length(missing_cols) > 0) {
    stop(sprintf("Manifest is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  manifest
}

#' Format One Trait to MRMOSS Standard Columns
#'
#' @param raw_file Raw summary-statistics file.
#' @param source_type Either "beta_se" or "formatted".
#' @param snp_col SNP column name.
#' @param a1_col Effect allele column name.
#' @param a2_col Other allele column name.
#' @param beta_col Beta column name.
#' @param se_col Standard error column name.
#' @param z_col Z-score column name.
#' @param p_col P-value column name.
#' @param n_col Sample size column name.
#' @param n_value Constant sample size when `n_col` is absent.
#' @return A data.table with columns SNP/A1/A2/Z/N/chi2/P.
mrmoss_format_trait <- function(raw_file,
                                source_type,
                                snp_col = NULL,
                                a1_col = NULL,
                                a2_col = NULL,
                                beta_col = NULL,
                                se_col = NULL,
                                z_col = NULL,
                                p_col = NULL,
                                n_col = NULL,
                                n_value = NULL) {
  dt <- data.table::fread(raw_file)

  source_type <- tolower(as.character(source_type)[1])
  if (source_type %in% c("formatted", "formatted_ready")) {
    needed <- c("SNP", "A1", "A2", "Z", "N", "P")
    miss <- setdiff(needed, colnames(dt))
    if (length(miss) > 0) {
      stop(sprintf("Formatted source %s missing columns: %s", raw_file, paste(miss, collapse = ", ")), call. = FALSE)
    }

    out <- data.table::data.table(
      SNP = as.character(dt[["SNP"]]),
      A1 = toupper(as.character(dt[["A1"]])),
      A2 = toupper(as.character(dt[["A2"]])),
      Z = suppressWarnings(as.numeric(dt[["Z"]])),
      N = suppressWarnings(as.numeric(dt[["N"]])),
      P = suppressWarnings(as.numeric(dt[["P"]]))
    )
  } else if (source_type %in% c("beta_se", "betase")) {
    snp_col <- mrmoss_to_character(snp_col)
    a1_col <- mrmoss_to_character(a1_col)
    a2_col <- mrmoss_to_character(a2_col)
    p_col <- mrmoss_to_character(p_col)

    if (is.null(snp_col) || is.null(a1_col) || is.null(a2_col) || is.null(p_col)) {
      stop("For source_type='beta_se', snp_col/a1_col/a2_col/p_col must be provided.", call. = FALSE)
    }

    z_col <- mrmoss_to_character(z_col)
    beta_col <- mrmoss_to_character(beta_col)
    se_col <- mrmoss_to_character(se_col)

    if (!is.null(z_col)) {
      if (!z_col %in% colnames(dt)) {
        stop(sprintf("z_col '%s' not found in %s", z_col, raw_file), call. = FALSE)
      }
      z <- suppressWarnings(as.numeric(dt[[z_col]]))
    } else {
      if (is.null(beta_col) || is.null(se_col)) {
        stop("Provide either z_col or beta_col+se_col.", call. = FALSE)
      }
      if (!beta_col %in% colnames(dt) || !se_col %in% colnames(dt)) {
        stop(sprintf("beta_col/se_col not found in %s", raw_file), call. = FALSE)
      }
      beta <- suppressWarnings(as.numeric(dt[[beta_col]]))
      se <- suppressWarnings(as.numeric(dt[[se_col]]))
      z <- beta / se
    }

    n_col <- mrmoss_to_character(n_col)
    n_value <- mrmoss_to_numeric(n_value)
    if (!is.null(n_col)) {
      if (!n_col %in% colnames(dt)) {
        stop(sprintf("n_col '%s' not found in %s", n_col, raw_file), call. = FALSE)
      }
      n <- suppressWarnings(as.numeric(dt[[n_col]]))
    } else if (!is.na(n_value)) {
      n <- rep(n_value, nrow(dt))
    } else {
      stop("Either n_col or n_value must be provided.", call. = FALSE)
    }

    cols_check <- c(snp_col, a1_col, a2_col, p_col)
    cols_miss <- cols_check[!cols_check %in% colnames(dt)]
    if (length(cols_miss) > 0) {
      stop(sprintf("Columns not found in %s: %s", raw_file, paste(cols_miss, collapse = ", ")), call. = FALSE)
    }

    out <- data.table::data.table(
      SNP = as.character(dt[[snp_col]]),
      A1 = toupper(as.character(dt[[a1_col]])),
      A2 = toupper(as.character(dt[[a2_col]])),
      Z = z,
      N = n,
      P = suppressWarnings(as.numeric(dt[[p_col]]))
    )
  } else {
    stop(sprintf("Unknown source_type '%s'", source_type), call. = FALSE)
  }

  out <- out[
    !is.na(out$SNP) & !is.na(out$A1) & !is.na(out$A2) &
      !is.na(out$Z) & !is.na(out$N) & !is.na(out$P),
    ,
    drop = FALSE
  ]
  out <- out[is.finite(out$Z) & is.finite(out$N) & is.finite(out$P) & out$N > 0, , drop = FALSE]
  out <- out[out$A1 %in% c("A", "C", "G", "T") & out$A2 %in% c("A", "C", "G", "T"), , drop = FALSE]

  # Keep one row per SNP to match downstream MR assumptions.
  out <- out[!duplicated(out$SNP), , drop = FALSE]
  out$chi2 <- out$Z^2

  out <- data.table::as.data.table(out)
  data.table::setcolorder(out, c("SNP", "A1", "A2", "Z", "N", "chi2", "P"))
  out
}

#' Batch Format Traits from Manifest
#'
#' @param manifest_csv Path to manifest CSV.
#' @param output_dir Output directory for formatted files.
#' @param traits Optional subset of trait names.
#' @param root_dir Optional root used for relative paths in manifest.
#' @param overwrite Whether to overwrite existing formatted files.
#' @param verbose Whether to print progress.
#' @return Summary data.table.
mrmoss_batch_format_manifest <- function(manifest_csv,
                                         output_dir,
                                         traits = NULL,
                                         root_dir = NULL,
                                         overwrite = FALSE,
                                         verbose = TRUE) {
  manifest <- as.data.frame(mrmoss_read_manifest(manifest_csv), stringsAsFactors = FALSE)
  if (!is.null(traits)) {
    manifest <- manifest[manifest$trait %in% traits, , drop = FALSE]
  }

  if (nrow(manifest) == 0) {
    stop("No traits selected from manifest.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summary_rows <- vector("list", nrow(manifest))

  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, , drop = FALSE]
    trait <- row$trait[[1]]
    raw_path <- mrmoss_resolve_path(row$raw_path[[1]], root_dir = root_dir)
    out_file <- file.path(output_dir, trait)

    if (file.exists(out_file) && !overwrite) {
      mrmoss_message(verbose, sprintf("[skip] %s exists", out_file))
      summary_rows[[i]] <- data.table::data.table(
        trait = trait,
        raw_file = raw_path,
        output_file = out_file,
        n_snp = NA_integer_,
        status = "skipped"
      )
      next
    }

    if (!file.exists(raw_path)) {
      stop(sprintf("Raw file does not exist for trait '%s': %s", trait, raw_path), call. = FALSE)
    }

    mrmoss_message(verbose, sprintf("[format] %s <- %s", trait, raw_path))

    formatted <- mrmoss_format_trait(
      raw_file = raw_path,
      source_type = row$source_type[[1]],
      snp_col = row$snp_col[[1]],
      a1_col = row$a1_col[[1]],
      a2_col = row$a2_col[[1]],
      beta_col = row$beta_col[[1]],
      se_col = row$se_col[[1]],
      z_col = row$z_col[[1]],
      p_col = row$p_col[[1]],
      n_col = row$n_col[[1]],
      n_value = row$n_value[[1]]
    )

    data.table::fwrite(formatted, out_file, sep = "\t")

    summary_rows[[i]] <- data.table::data.table(
      trait = trait,
      raw_file = raw_path,
      output_file = out_file,
      n_snp = nrow(formatted),
      status = "ok"
    )
  }

  data.table::rbindlist(summary_rows, fill = TRUE)
}

#' Load Formatted Trait Files
#'
#' @param formatted_dir Directory containing formatted trait files.
#' @param traits Character vector of trait names.
#' @return Named list of data.tables.
mrmoss_load_formatted_traits <- function(formatted_dir, traits) {
  out <- vector("list", length(traits))
  names(out) <- traits

  for (trait in traits) {
    f <- file.path(formatted_dir, trait)
    if (!file.exists(f)) {
      stop(sprintf("Formatted trait file is missing: %s", f), call. = FALSE)
    }
    dt <- data.table::fread(f)
    needed <- c("SNP", "A1", "A2", "Z", "N", "P")
    miss <- setdiff(needed, colnames(dt))
    if (length(miss) > 0) {
      stop(sprintf("Formatted trait file %s missing columns: %s", f, paste(miss, collapse = ", ")), call. = FALSE)
    }
    out[[trait]] <- data.table::data.table(
      SNP = as.character(dt$SNP),
      A1 = toupper(as.character(dt$A1)),
      A2 = toupper(as.character(dt$A2)),
      Z = as.numeric(dt$Z),
      N = as.numeric(dt$N),
      P = as.numeric(dt$P)
    )
  }

  out
}
