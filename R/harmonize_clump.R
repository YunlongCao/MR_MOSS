mrmoss_complement_allele <- function(allele) {
  ifelse(
    allele == "A", "T",
    ifelse(
      allele == "T", "A",
      ifelse(allele == "G", "C", ifelse(allele == "C", "G", NA_character_))
    )
  )
}

#' Harmonize Exposure and Outcome Summary Statistics
#'
#' @param dat1 Exposure data with columns SNP/A1/A2/Z/N/P.
#' @param dat2 Outcome data with columns SNP/A1/A2/Z/N/P.
#' @return Harmonized data for MR estimation.
mrmoss_harmonize_data <- function(dat1, dat2) {
  dat <- merge(dat1, dat2, by = "SNP", suffixes = c(".exp", ".out"))

  a1e <- toupper(dat$A1.exp)
  a2e <- toupper(dat$A2.exp)
  a1o <- toupper(dat$A1.out)
  a2o <- toupper(dat$A2.out)

  same <- a1e == a1o & a2e == a2o
  flip <- a1e == a2o & a2e == a1o

  c1o <- mrmoss_complement_allele(a1o)
  c2o <- mrmoss_complement_allele(a2o)

  same_comp <- a1e == c1o & a2e == c2o
  flip_comp <- a1e == c2o & a2e == c1o

  keep <- same | flip | same_comp | flip_comp
  flip_sign <- flip | flip_comp

  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0) {
    return(data.table::data.table())
  }

  z_out <- dat$Z.out
  z_out[flip_sign[keep]] <- -z_out[flip_sign[keep]]

  out <- data.table::data.table(
    SNP = dat$SNP,
    A1 = dat$A1.exp,
    A2 = dat$A2.exp,
    b.exp = dat$Z.exp / sqrt(dat$N.exp),
    b.out = z_out / sqrt(dat$N.out),
    se.exp = 1 / sqrt(dat$N.exp),
    se.out = 1 / sqrt(dat$N.out),
    pval.exp = dat$P.exp,
    pval.out = dat$P.out
  )

  out <- out[
    is.finite(out$b.exp) & is.finite(out$b.out) &
      is.finite(out$se.exp) & is.finite(out$se.out),
    ,
    drop = FALSE
  ]
  out
}

#' LD Clumping for Instrument Selection
#'
#' @param dat Harmonized data from `mrmoss_harmonize_data`.
#' @param iv_threshold IV p-value threshold on exposure.
#' @param reference_prefix Prefix to PLINK reference panel (e.g. "/path/EUR").
#'   If `NULL`, use the online OpenGWAS clumping service via `ieugwasr`.
#' @param plink_bin Path to PLINK binary for local clumping.
#' @param pop Population code for ieugwasr (default "EUR").
#' @param clump_kb Clumping window in kb.
#' @param clump_r2 Clumping r2 threshold.
#' @param clump_p Secondary clumping p threshold.
#' @param f_min Minimum F statistic.
#' @return Clumped MR dataset.
mrmoss_clump_data <- function(dat,
                              iv_threshold,
                              reference_prefix = NULL,
                              plink_bin = NULL,
                              pop = "EUR",
                              clump_kb = 1000,
                              clump_r2 = 0.001,
                              clump_p = 0.999,
                              f_min = 10) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) {
    stop("Package 'ieugwasr' is required for LD clumping.", call. = FALSE)
  }

  if (nrow(dat) == 0) {
    return(dat)
  }

  dat <- data.table::as.data.table(dat)
  dat <- dat[dat$pval.exp <= iv_threshold, , drop = FALSE]
  if (nrow(dat) == 0) {
    return(dat)
  }

  dat$Fexp <- (dat$b.exp / dat$se.exp)^2

  clump_input <- data.frame(rsid = dat$SNP, pval = dat$pval.exp)
  if (!is.null(reference_prefix) && nzchar(reference_prefix)) {
    reference_prefix <- mrmoss_resolve_path(reference_prefix)
    if (!file.exists(paste0(reference_prefix, ".bed"))) {
      stop(
        sprintf("Reference bed file not found for local clumping: %s.bed", reference_prefix),
        call. = FALSE
      )
    }

    plink_bin <- mrmoss_get_plink_binary(plink_bin)
    clumped <- ieugwasr::ld_clump(
      clump_input,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = clump_p,
      pop = pop,
      bfile = reference_prefix,
      plink_bin = plink_bin
    )
  } else {
    clumped <- tryCatch(
      ieugwasr::ld_clump(
        clump_input,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = clump_p,
        pop = pop
      ),
      error = function(e) {
        stop(
          paste0(
            "Online clumping via OpenGWAS failed. ",
            "Set OPENGWAS_JWT for ieugwasr, or use local clumping with `reference_prefix`. ",
            "Original error: ", conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
  }

  out <- dat[dat$SNP %in% clumped$rsid, , drop = FALSE]
  out <- out[out$Fexp >= f_min, , drop = FALSE]
  out
}

#' Estimate Outcome Correlation Matrix from Null SNPs
#'
#' @param outcome_data_list Named list of formatted outcome data.tables.
#' @param p_null_threshold Null-SNP p-value threshold.
#' @return Positive-definite correlation matrix.
mrmoss_estimate_outcome_correlation <- function(outcome_data_list, p_null_threshold = 1e-5) {
  if (length(outcome_data_list) < 2) {
    return(matrix(1, nrow = length(outcome_data_list), ncol = length(outcome_data_list)))
  }

  common_snps <- Reduce(intersect, lapply(outcome_data_list, function(x) x$SNP))
  if (length(common_snps) == 0) {
    stop("No common SNPs across outcomes for correlation estimation.", call. = FALSE)
  }

  z_tables <- lapply(seq_along(outcome_data_list), function(i) {
    dt0 <- outcome_data_list[[i]]
    keep <- dt0$SNP %in% common_snps & dt0$P > p_null_threshold
    dt <- data.table::data.table(SNP = dt0$SNP[keep], Z = dt0$Z[keep])
    data.table::setnames(dt, "Z", paste0("Z", i))
    dt
  })

  merged <- Reduce(function(x, y) merge(x, y, by = "SNP", all = FALSE), z_tables)
  if (nrow(merged) < 10) {
    stop("Too few null SNPs remained for correlation estimation.", call. = FALSE)
  }

  z_matrix <- as.matrix(merged[, -1, drop = FALSE])
  R <- stats::cor(z_matrix, use = "pairwise.complete.obs")
  mrmoss_make_spd(R)
}
