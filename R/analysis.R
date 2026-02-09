mrmoss_other_methods <- function(b_exp, b_out, n1, n2) {
  p <- length(b_exp)
  out <- list(
    IVW_estimate = NA_real_, IVW_se = NA_real_, IVW_pvalue = NA_real_,
    RAPS_estimate = NA_real_, RAPS_se = NA_real_, RAPS_pvalue = NA_real_,
    MRMix_estimate = NA_real_, MRMix_se = NA_real_, MRMix_pvalue = NA_real_,
    Egger_estimate = NA_real_, Egger_se = NA_real_, Egger_pvalue = NA_real_
  )

  if (requireNamespace("TwoSampleMR", quietly = TRUE)) {
    ivw <- tryCatch(
      TwoSampleMR::mr_ivw(b_exp, b_out, rep(1 / sqrt(n1), p), rep(1 / sqrt(n2), p)),
      error = function(e) NULL
    )
    if (!is.null(ivw)) {
      out$IVW_estimate <- ivw$b
      out$IVW_se <- ivw$se
      out$IVW_pvalue <- ivw$pval
    }

    egger <- tryCatch(
      TwoSampleMR::mr_egger_regression(b_exp, b_out, rep(1 / sqrt(n1), p), rep(1 / sqrt(n2), p)),
      error = function(e) NULL
    )
    if (!is.null(egger)) {
      out$Egger_estimate <- egger$b
      out$Egger_se <- egger$se
      out$Egger_pvalue <- egger$pval
    }
  }

  if (requireNamespace("mr.raps", quietly = TRUE)) {
    raps_df <- data.frame(
      beta.exposure = b_exp,
      beta.outcome = b_out,
      se.exposure = rep(1 / sqrt(n1), p),
      se.outcome = rep(1 / sqrt(n2), p)
    )

    raps <- tryCatch(mr.raps::mr.raps(raps_df, diagnostics = FALSE), error = function(e) NULL)
    if (!is.null(raps)) {
      out$RAPS_estimate <- raps$beta.hat
      out$RAPS_se <- raps$beta.se
      out$RAPS_pvalue <- 2 * stats::pnorm(abs(raps$beta.hat / raps$beta.se), lower.tail = FALSE)
    }
  }

  if (requireNamespace("MRMix", quietly = TRUE)) {
    mrmix <- tryCatch(
      MRMix::MRMix(b_exp, b_out, rep(1 / sqrt(n1), p), rep(1 / sqrt(n2), p)),
      error = function(e) NULL
    )
    if (!is.null(mrmix)) {
      out$MRMix_estimate <- mrmix$theta
      out$MRMix_se <- mrmix$SE_theta
      out$MRMix_pvalue <- mrmix$pvalue_theta
    }
  }

  out
}

mrmoss_prepare_outcome_specific_mrdat <- function(exposure_dat,
                                                  outcome_dat,
                                                  iv_threshold,
                                                  reference_prefix,
                                                  plink_bin,
                                                  pop,
                                                  clump_kb,
                                                  clump_r2,
                                                  clump_p,
                                                  f_min) {
  harmonized <- mrmoss_harmonize_data(exposure_dat, outcome_dat)
  if (nrow(harmonized) == 0) {
    return(harmonized)
  }

  mrmoss_clump_data(
    dat = harmonized,
    iv_threshold = iv_threshold,
    reference_prefix = reference_prefix,
    plink_bin = plink_bin,
    pop = pop,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = clump_p,
    f_min = f_min
  )
}

mrmoss_standardize_formatted_trait <- function(dt, trait_name = NULL) {
  dt <- data.table::as.data.table(dt)
  needed <- c("SNP", "A1", "A2", "Z", "N", "P")
  miss <- setdiff(needed, colnames(dt))
  if (length(miss) > 0) {
    who <- if (!is.null(trait_name) && nzchar(trait_name)) sprintf(" '%s'", trait_name) else ""
    stop(sprintf("Formatted trait%s is missing columns: %s", who, paste(miss, collapse = ", ")), call. = FALSE)
  }

  out <- data.table::data.table(
    SNP = as.character(dt$SNP),
    A1 = toupper(as.character(dt$A1)),
    A2 = toupper(as.character(dt$A2)),
    Z = as.numeric(dt$Z),
    N = as.numeric(dt$N),
    P = as.numeric(dt$P)
  )

  keep <- !is.na(out$SNP) & !is.na(out$A1) & !is.na(out$A2) &
    is.finite(out$Z) & is.finite(out$N) & is.finite(out$P) &
    out$N > 0 &
    out$A1 %in% c("A", "C", "G", "T") &
    out$A2 %in% c("A", "C", "G", "T")
  out <- out[keep, ]
  out
}

#' Build MR-MOSS Input from Formatted Data
#'
#' Builds the clumped and harmonized MR-MOSS input for one exposure and
#' multiple outcomes.
#'
#' @param exposure Exposure trait name.
#' @param outcomes Outcome trait names.
#' @param formatted_dir Directory with formatted trait files. Ignored when
#'   `formatted_data` is provided.
#' @param formatted_data Optional named list of formatted trait data.frames or
#'   data.tables. Each element must contain `SNP`, `A1`, `A2`, `Z`, `N`, `P`.
#' @param iv_threshold Instrument p-value threshold for clumping.
#' @param reference_prefix Optional PLINK reference prefix for local clumping.
#'   If `NULL`, online OpenGWAS clumping is used via `ieugwasr`.
#' @param plink_bin Optional PLINK binary path for local clumping.
#' @param pop Population code for clumping.
#' @param clump_kb LD window size.
#' @param clump_r2 LD r2 threshold.
#' @param clump_p Secondary clumping p-value threshold.
#' @param f_min Minimum F statistic for instruments.
#' @param p_null_threshold Threshold for null SNPs in outcome-correlation estimation.
#' @param n2 Optional outcome sample size. Default is mean outcome N.
#' @param cache_dir Optional cache directory to save per-outcome clumped inputs.
#' @param verbose Whether to print progress.
#' @return A list containing `gamma_hat`, `Gamma_hat`, `R`, and metadata.
#' @export
mrmoss_input <- function(exposure,
                         outcomes,
                         formatted_dir = NULL,
                         formatted_data = NULL,
                         iv_threshold = 5e-7,
                         reference_prefix = NULL,
                         plink_bin = NULL,
                         pop = "EUR",
                         clump_kb = 1000,
                         clump_r2 = 0.001,
                         clump_p = 0.999,
                         f_min = 10,
                         p_null_threshold = 1e-5,
                         n2 = NULL,
                         cache_dir = NULL,
                         verbose = TRUE) {
  exposure <- as.character(exposure)[1]
  outcomes <- unique(as.character(outcomes))
  if (!nzchar(exposure) || length(outcomes) == 0) {
    stop("`exposure` must be one trait name and `outcomes` must be non-empty.", call. = FALSE)
  }

  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(reference_prefix) && nzchar(reference_prefix)) {
    reference_prefix <- mrmoss_resolve_path(reference_prefix)
    if (!file.exists(paste0(reference_prefix, ".bed"))) {
      stop(sprintf("Reference bed file not found: %s.bed", reference_prefix), call. = FALSE)
    }
    plink_bin <- mrmoss_get_plink_binary(plink_bin)
    mrmoss_message(verbose, sprintf("[clump] local reference: %s", reference_prefix))
  } else {
    reference_prefix <- NULL
    plink_bin <- NULL
    mrmoss_message(verbose, sprintf("[clump] online OpenGWAS service (pop = %s)", pop))
    if (!nzchar(Sys.getenv("OPENGWAS_JWT", unset = ""))) {
      mrmoss_message(
        verbose,
        "[clump] OPENGWAS_JWT is not set. If online clumping fails with 401, set OPENGWAS_JWT or use local reference_prefix."
      )
    }
  }

  if (!is.null(formatted_data)) {
    if (!is.list(formatted_data) || is.null(names(formatted_data))) {
      stop("`formatted_data` must be a named list.", call. = FALSE)
    }

    need_traits <- unique(c(exposure, outcomes))
    missing_traits <- setdiff(need_traits, names(formatted_data))
    if (length(missing_traits) > 0) {
      stop(
        sprintf("`formatted_data` is missing required traits: %s", paste(missing_traits, collapse = ", ")),
        call. = FALSE
      )
    }

    trait_data <- lapply(need_traits, function(tr) {
      mrmoss_standardize_formatted_trait(formatted_data[[tr]], trait_name = tr)
    })
    names(trait_data) <- need_traits

    exp_dat <- trait_data[[exposure]]
    outcome_data <- trait_data[outcomes]
    mrmoss_message(verbose, sprintf("[input] using in-memory formatted_data for %d traits", length(need_traits)))
  } else {
    if (is.null(formatted_dir) || !nzchar(as.character(formatted_dir)[1])) {
      stop("Provide either `formatted_dir` or `formatted_data`.", call. = FALSE)
    }
    exp_dat <- mrmoss_load_formatted_traits(formatted_dir, exposure)[[1]]
    outcome_data <- mrmoss_load_formatted_traits(formatted_dir, outcomes)
  }

  R <- as.matrix(mrmoss_estimate_outcome_correlation(
    outcome_data,
    p_null_threshold = p_null_threshold
  ))
  storage.mode(R) <- "double"
  if (nrow(R) != length(outcomes) || ncol(R) != length(outcomes)) {
    stop(
      sprintf(
        "Outcome correlation matrix dimension mismatch: expected %dx%d, got %dx%d.",
        length(outcomes), length(outcomes), nrow(R), ncol(R)
      ),
      call. = FALSE
    )
  }
  rownames(R) <- outcomes
  colnames(R) <- outcomes

  if (is.null(n2)) {
    n2 <- floor(mean(vapply(outcome_data, function(x) x$N[1], numeric(1)), na.rm = TRUE))
  }

  n1 <- as.integer(round(exp_dat$N[1]))
  mrdat_list <- vector("list", length(outcomes))
  names(mrdat_list) <- outcomes

  for (outcome in outcomes) {
    out_dat <- outcome_data[[outcome]]
    mrdat <- mrmoss_prepare_outcome_specific_mrdat(
      exposure_dat = exp_dat,
      outcome_dat = out_dat,
      iv_threshold = iv_threshold,
      reference_prefix = reference_prefix,
      plink_bin = plink_bin,
      pop = pop,
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p = clump_p,
      f_min = f_min
    )
    mrdat_list[[outcome]] <- mrdat

    if (!is.null(cache_dir)) {
      cache_file <- file.path(
        cache_dir,
        paste0("MRdat_", exposure, "_", outcome, "_", format(iv_threshold, scientific = TRUE), ".rds")
      )
      saveRDS(mrdat, cache_file)
    }
  }

  zero_outcomes <- names(mrdat_list)[vapply(mrdat_list, nrow, integer(1)) == 0]
  if (length(zero_outcomes) > 0) {
    stop(
      sprintf(
        "At least one outcome has zero instruments after clumping: %s",
        paste(zero_outcomes, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  common_snps <- Reduce(intersect, lapply(mrdat_list, function(x) x$SNP))
  if (length(common_snps) < 3) {
    stop(sprintf("Too few shared instruments after clumping: %d", length(common_snps)), call. = FALSE)
  }
  common_snps <- as.character(common_snps)

  exp_idx <- match(common_snps, as.character(mrdat_list[[1]]$SNP))
  if (anyNA(exp_idx)) {
    stop("Internal error: exposure SNP alignment failed after intersection.", call. = FALSE)
  }
  b_exp <- as.numeric(mrdat_list[[1]]$b.exp[exp_idx])

  b_out <- vapply(
    mrdat_list,
    FUN = function(dt) {
      out_idx <- match(common_snps, as.character(dt$SNP))
      if (anyNA(out_idx)) {
        stop("Internal error: outcome SNP alignment failed after intersection.", call. = FALSE)
      }
      as.numeric(dt$b.out[out_idx])
    },
    FUN.VALUE = numeric(length(common_snps))
  )
  b_out <- as.matrix(b_out)
  storage.mode(b_out) <- "double"

  out <- list(
    exposure = exposure,
    outcomes = outcomes,
    iv_threshold = iv_threshold,
    NO_of_IVs = length(common_snps),
    common_snps = common_snps,
    n1 = as.integer(round(n1)),
    n2 = as.integer(round(n2)),
    gamma_hat = b_exp,
    Gamma_hat = b_out,
    R = R
  )
  class(out) <- "mrmoss_input"
  out
}

#' Run Core MR-MOSS Model on Prepared Input
#'
#' @param exposure Exposure trait name.
#' @param outcomes Outcome trait names.
#' @param mrmoss_input Prepared input object from `mrmoss_input()`.
#' @param pvalue_output One of `"both"`, `"outcome"`, `"global"`.
#' @param maxiter Maximum PX-EM iterations.
#' @param output_dir Optional output directory.
#' @param output_prefix Output basename when writing files.
#' @param verbose Whether to print progress.
#' @return A list with result table, raw fit object, and optional output files.
#' @export
mrmoss <- function(exposure,
                   outcomes,
                   mrmoss_input,
                   pvalue_output = c("both", "outcome", "global"),
                   maxiter = 1000000,
                   output_dir = NULL,
                   output_prefix = "mrmoss_result",
                   verbose = TRUE) {
  pvalue_output <- match.arg(pvalue_output)
  exposure <- as.character(exposure)[1]
  outcomes <- as.character(outcomes)

  if (!is.list(mrmoss_input)) {
    stop("`mrmoss_input` must be a list created by `mrmoss_input()`.", call. = FALSE)
  }
  required <- c("exposure", "outcomes", "gamma_hat", "Gamma_hat", "R", "n1", "n2", "NO_of_IVs")
  miss <- setdiff(required, names(mrmoss_input))
  if (length(miss) > 0) {
    stop(sprintf("`mrmoss_input` is missing fields: %s", paste(miss, collapse = ", ")), call. = FALSE)
  }
  if (!identical(exposure, as.character(mrmoss_input$exposure))) {
    stop(sprintf("Exposure mismatch: '%s' vs input '%s'.", exposure, as.character(mrmoss_input$exposure)), call. = FALSE)
  }

  input_outcomes <- as.character(mrmoss_input$outcomes)
  out_idx <- match(outcomes, input_outcomes)
  if (length(outcomes) == 0 || anyNA(out_idx)) {
    stop("`outcomes` must be a non-empty subset in the same naming space as `mrmoss_input$outcomes`.", call. = FALSE)
  }

  gamma_hat <- as.numeric(mrmoss_input$gamma_hat)
  Gamma_hat <- as.matrix(mrmoss_input$Gamma_hat)[, out_idx, drop = FALSE]
  R <- as.matrix(mrmoss_input$R)[out_idx, out_idx, drop = FALSE]
  n1 <- as.integer(round(mrmoss_input$n1))
  n2 <- as.integer(round(mrmoss_input$n2))

  m <- length(outcomes)
  theta <- c(rep(0, m), 0.001, rep(0.001, m), 1, rep(0.8, m))
  test <- as.integer(seq_len(m))

  fit <- MRMOSS_PX_cpp(
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    R = R,
    n1 = n1,
    n2 = n2,
    theta0 = theta,
    test = test,
    maxiter = as.integer(maxiter)
  )

  row <- data.table::data.table(
    exposure = exposure,
    IV_Threshold = as.numeric(mrmoss_input$iv_threshold),
    NO_of_IVs = as.integer(mrmoss_input$NO_of_IVs),
    iteration = as.integer(fit$iteration)
  )
  if (pvalue_output %in% c("both", "global")) {
    row$Overall_pvalue <- as.numeric(fit$pvalue_overall)
  }

  for (k in seq_along(outcomes)) {
    row[[paste0("outcome.", k)]] <- outcomes[k]
    row[[paste0("MOSS_estimates.", k)]] <- as.numeric(fit$beta[k])
    if (pvalue_output %in% c("both", "outcome")) {
      row[[paste0("pvalue.", k)]] <- as.numeric(fit$pvalue[k])
    }
  }

  files <- list(mrmoss = NULL, correlation = NULL)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(output_dir, paste0(output_prefix, ".txt"))
    r_file <- file.path(output_dir, paste0(output_prefix, "_R_matrix.tsv"))
    data.table::fwrite(row, out_file, sep = "\t")
    data.table::fwrite(data.table::as.data.table(R, keep.rownames = "outcome"), r_file, sep = "\t")
    files$mrmoss <- out_file
    files$correlation <- r_file
    mrmoss_message(verbose, sprintf("[output] wrote %s", out_file))
  }

  list(
    result = row,
    fit = fit,
    pvalue_output = pvalue_output,
    files = files
  )
}

#' Run MR-MOSS Real-Data Analysis
#'
#' @param exposures Exposure trait names.
#' @param outcomes Outcome trait names.
#' @param formatted_dir Directory with formatted trait files.
#'   Ignored when `formatted_data` is provided.
#' @param formatted_data Optional named list of formatted trait data.frames or
#'   data.tables. Each element must contain columns `SNP`, `A1`, `A2`, `Z`,
#'   `N`, `P`.
#' @param reference_prefix Optional PLINK reference prefix (without
#'   `.bed/.bim/.fam`) for local clumping.
#'   If `NULL`, LD clumping uses the online OpenGWAS service via `ieugwasr`.
#' @param output_dir Output directory.
#' @param output_prefix Output basename (without extension).
#' @param iv_thresholds Instrument p-value thresholds.
#' @param n2 Optional outcome sample size. Default is mean outcome N.
#' @param plink_bin Optional PLINK binary path for local clumping.
#' @param pop Population code.
#' @param clump_kb LD window size.
#' @param clump_r2 LD r2 threshold.
#' @param clump_p Secondary p threshold for clumping.
#' @param f_min Minimum F statistic for instruments.
#' @param p_null_threshold Threshold for null SNPs when estimating outcome correlation.
#' @param maxiter Maximum PX-EM iterations.
#' @param include_other_methods Whether to compute IVW/RAPS/Egger/MRMix.
#' @param cache_dir Optional directory for saving MRdat intermediate files.
#' @param verbose Print progress.
#' @return List with MR-MOSS table, other-method table, and output paths.
mrmoss_run_analysis <- function(exposures,
                                outcomes,
                                formatted_dir = NULL,
                                formatted_data = NULL,
                                reference_prefix = NULL,
                                output_dir,
                                output_prefix,
                                iv_thresholds = c(5e-07, 5e-08),
                                n2 = NULL,
                                plink_bin = NULL,
                                pop = "EUR",
                                clump_kb = 1000,
                                clump_r2 = 0.001,
                                clump_p = 0.999,
                                f_min = 10,
                                p_null_threshold = 1e-5,
                                maxiter = 1000000,
                                include_other_methods = TRUE,
                                cache_dir = NULL,
                                verbose = TRUE) {
  exposures <- unique(as.character(exposures))
  outcomes <- unique(as.character(outcomes))
  if (length(exposures) == 0 || length(outcomes) == 0) {
    stop("`exposures` and `outcomes` must each contain at least one trait name.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!is.null(reference_prefix) && nzchar(reference_prefix)) {
    reference_prefix <- mrmoss_resolve_path(reference_prefix)
    if (!file.exists(paste0(reference_prefix, ".bed"))) {
      stop(sprintf("Reference bed file not found: %s.bed", reference_prefix), call. = FALSE)
    }
    plink_bin <- mrmoss_get_plink_binary(plink_bin)
    mrmoss_message(verbose, sprintf("[clump] local reference: %s", reference_prefix))
  } else {
    reference_prefix <- NULL
    plink_bin <- NULL
    mrmoss_message(verbose, sprintf("[clump] online OpenGWAS service (pop = %s)", pop))
    if (!nzchar(Sys.getenv("OPENGWAS_JWT", unset = ""))) {
      mrmoss_message(
        verbose,
        "[clump] OPENGWAS_JWT is not set. If online clumping fails with 401, set OPENGWAS_JWT or use local reference_prefix."
      )
    }
  }

  if (!is.null(formatted_data)) {
    if (!is.list(formatted_data) || is.null(names(formatted_data))) {
      stop("`formatted_data` must be a named list.", call. = FALSE)
    }

    need_traits <- unique(c(exposures, outcomes))
    missing_traits <- setdiff(need_traits, names(formatted_data))
    if (length(missing_traits) > 0) {
      stop(
        sprintf(
          "`formatted_data` is missing required traits: %s",
          paste(missing_traits, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    trait_data <- lapply(need_traits, function(tr) {
      mrmoss_standardize_formatted_trait(formatted_data[[tr]], trait_name = tr)
    })
    names(trait_data) <- need_traits

    exposure_data <- trait_data[exposures]
    outcome_data <- trait_data[outcomes]
    mrmoss_message(verbose, sprintf("[input] using in-memory formatted_data for %d traits", length(need_traits)))
  } else {
    if (is.null(formatted_dir) || !nzchar(as.character(formatted_dir)[1])) {
      stop("Provide either `formatted_dir` or `formatted_data`.", call. = FALSE)
    }
    exposure_data <- mrmoss_load_formatted_traits(formatted_dir, exposures)
    outcome_data <- mrmoss_load_formatted_traits(formatted_dir, outcomes)
  }

  R <- as.matrix(mrmoss_estimate_outcome_correlation(
    outcome_data,
    p_null_threshold = p_null_threshold
  ))
  storage.mode(R) <- "double"
  if (!is.matrix(R)) {
    stop("Outcome correlation object `R` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(R) != length(outcomes) || ncol(R) != length(outcomes)) {
    stop(
      sprintf(
        "Outcome correlation matrix dimension mismatch: expected %dx%d, got %dx%d.",
        length(outcomes), length(outcomes), nrow(R), ncol(R)
      ),
      call. = FALSE
    )
  }
  colnames(R) <- outcomes
  rownames(R) <- outcomes

  if (is.null(n2)) {
    n2 <- floor(mean(vapply(outcome_data, function(x) x$N[1], numeric(1)), na.rm = TRUE))
  }

  m <- length(outcomes)
  theta <- c(rep(0, m), 0.001, rep(0.001, m), 1, rep(0.8, m))
  test <- as.integer(seq_len(m))

  moss_rows <- list()
  other_rows <- list()

  for (exposure in exposures) {
    mrmoss_message(verbose, sprintf("[exposure] %s", exposure))
    exp_dat <- exposure_data[[exposure]]
    n1 <- as.integer(round(exp_dat$N[1]))

    for (thr in iv_thresholds) {
      mrmoss_message(verbose, sprintf("  [threshold] %.1e", thr))

      mrdat_list <- vector("list", length(outcomes))
      names(mrdat_list) <- outcomes

      for (outcome in outcomes) {
        out_dat <- outcome_data[[outcome]]

        mrdat <- mrmoss_prepare_outcome_specific_mrdat(
          exposure_dat = exp_dat,
          outcome_dat = out_dat,
          iv_threshold = thr,
          reference_prefix = reference_prefix,
          plink_bin = plink_bin,
          pop = pop,
          clump_kb = clump_kb,
          clump_r2 = clump_r2,
          clump_p = clump_p,
          f_min = f_min
        )

        mrdat_list[[outcome]] <- mrdat

        if (!is.null(cache_dir)) {
          cache_file <- file.path(cache_dir, paste0("MRdat_", exposure, "_", outcome, "_", format(thr, scientific = TRUE), ".rds"))
          saveRDS(mrdat, cache_file)
        }
      }

      if (any(vapply(mrdat_list, nrow, integer(1)) == 0)) {
        mrmoss_message(verbose, "    [skip] at least one outcome has zero instruments after clumping")
        next
      }

      common_snps <- Reduce(intersect, lapply(mrdat_list, function(x) x$SNP))
      if (length(common_snps) < 3) {
        mrmoss_message(verbose, sprintf("    [skip] too few shared instruments (%d)", length(common_snps)))
        next
      }

      common_snps <- as.character(common_snps)

      exp_idx <- match(common_snps, as.character(mrdat_list[[1]]$SNP))
      if (anyNA(exp_idx)) {
        stop("Internal error: exposure SNP alignment failed after intersection.", call. = FALSE)
      }
      b_exp <- as.numeric(mrdat_list[[1]]$b.exp[exp_idx])

      b_out <- vapply(
        mrdat_list,
        FUN = function(dt) {
          out_idx <- match(common_snps, as.character(dt$SNP))
          if (anyNA(out_idx)) {
            stop("Internal error: outcome SNP alignment failed after intersection.", call. = FALSE)
          }
          as.numeric(dt$b.out[out_idx])
        },
        FUN.VALUE = numeric(length(common_snps))
      )
      b_out <- as.matrix(b_out)
      storage.mode(b_out) <- "double"
      if (!is.matrix(b_out)) {
        stop("Internal error: Gamma_hat is not a matrix after alignment.", call. = FALSE)
      }
      if (nrow(b_out) != length(b_exp)) {
        stop(
          sprintf(
            "Internal error: gamma_hat/Gamma_hat row mismatch (%d vs %d).",
            length(b_exp), nrow(b_out)
          ),
          call. = FALSE
        )
      }
      if (ncol(b_out) != length(outcomes)) {
        stop(
          sprintf(
            "Internal error: Gamma_hat column mismatch (expected %d outcomes, got %d).",
            length(outcomes), ncol(b_out)
          ),
          call. = FALSE
        )
      }

      fit <- MRMOSS_PX_cpp(
        gamma_hat = b_exp,
        Gamma_hat = b_out,
        R = R,
        n1 = n1,
        n2 = as.integer(round(n2)),
        theta0 = theta,
        test = test,
        maxiter = as.integer(maxiter)
      )

      moss_row <- data.table::data.table(
        exposure = exposure,
        IV_Threshold = thr,
        NO_of_IVs = length(common_snps),
        Overall_pvalue = fit$pvalue_overall,
        iteration = fit$iteration
      )

      for (k in seq_along(outcomes)) {
        moss_row[[paste0("outcome.", k)]] <- outcomes[k]
        moss_row[[paste0("MOSS_estimates.", k)]] <- as.numeric(fit$beta[k])
        moss_row[[paste0("pvalue.", k)]] <- as.numeric(fit$pvalue[k])
      }

      moss_rows[[length(moss_rows) + 1]] <- moss_row

      if (isTRUE(include_other_methods)) {
        for (k in seq_along(outcomes)) {
          other <- mrmoss_other_methods(
            b_exp = b_exp,
            b_out = b_out[, k],
            n1 = n1,
            n2 = n2
          )

          other_rows[[length(other_rows) + 1]] <- data.table::data.table(
            exposure = exposure,
            outcome = outcomes[k],
            IV_Threshold = thr,
            NO_of_IVs = length(common_snps),
            IVW_estimate = other$IVW_estimate,
            IVW_se = other$IVW_se,
            IVW_pvalue = other$IVW_pvalue,
            RAPS_estimate = other$RAPS_estimate,
            RAPS_se = other$RAPS_se,
            RAPS_pvalue = other$RAPS_pvalue,
            MRMix_estimate = other$MRMix_estimate,
            MRMix_se = other$MRMix_se,
            MRMix_pvalue = other$MRMix_pvalue,
            Egger_estimate = other$Egger_estimate,
            Egger_se = other$Egger_se,
            Egger_pvalue = other$Egger_pvalue
          )
        }
      }
    }
  }

  moss_tab <- if (length(moss_rows) > 0) data.table::rbindlist(moss_rows, fill = TRUE) else data.table::data.table()
  other_tab <- if (length(other_rows) > 0) data.table::rbindlist(other_rows, fill = TRUE) else data.table::data.table()

  moss_file <- NULL
  other_file <- NULL
  r_file <- file.path(output_dir, paste0(output_prefix, "_R_matrix.tsv"))

  if (nrow(moss_tab) > 0) {
    moss_file <- file.path(output_dir, paste0(output_prefix, ".txt"))
    data.table::fwrite(moss_tab, moss_file, sep = "\t")
  } else {
    mrmoss_message(verbose, "No MR-MOSS rows were produced; returning mrmoss output path as NULL.")
  }
  if (nrow(other_tab) > 0) {
    other_file <- file.path(output_dir, paste0(output_prefix, "_othermethods.txt"))
    data.table::fwrite(other_tab, other_file, sep = "\t")
  } else if (isTRUE(include_other_methods)) {
    mrmoss_message(verbose, "No comparison-method rows were produced; returning other_methods path as NULL.")
  }
  data.table::fwrite(data.table::as.data.table(R, keep.rownames = "outcome"), r_file, sep = "\t")

  list(
    mrmoss = moss_tab,
    other_methods = other_tab,
    correlation = R,
    files = list(mrmoss = moss_file, other_methods = other_file, correlation = r_file)
  )
}

#' Run One Built-in Manuscript Profile
#'
#' @param profile Profile name in `mrmoss_profiles()`.
#' @param formatted_dir Directory with formatted trait files.
#' @param reference_prefix Optional PLINK reference prefix for local clumping.
#' @param output_dir Output directory for this profile.
#' @param manifest_csv Optional manifest to auto-format missing traits before running.
#' @param root_dir Optional root directory for relative paths in manifest.
#' @param overwrite_format Whether to overwrite existing formatted files.
#' @param ... Extra arguments passed to `mrmoss_run_analysis`.
#' @return Result from `mrmoss_run_analysis`.
mrmoss_run_profile <- function(profile,
                               formatted_dir,
                               reference_prefix = NULL,
                               output_dir,
                               manifest_csv = NULL,
                               root_dir = NULL,
                               overwrite_format = FALSE,
                               ...) {
  profiles <- mrmoss_profiles()
  if (!profile %in% names(profiles)) {
    stop(sprintf("Unknown profile '%s'. Available: %s", profile, paste(names(profiles), collapse = ", ")), call. = FALSE)
  }

  cfg <- profiles[[profile]]

  if (!is.null(manifest_csv)) {
    need_traits <- unique(c(cfg$exposures, cfg$outcomes))
    mrmoss_batch_format_manifest(
      manifest_csv = manifest_csv,
      output_dir = formatted_dir,
      traits = need_traits,
      root_dir = root_dir,
      overwrite = overwrite_format,
      verbose = TRUE
    )
  }

  mrmoss_run_analysis(
    exposures = cfg$exposures,
    outcomes = cfg$outcomes,
    formatted_dir = formatted_dir,
    reference_prefix = reference_prefix,
    output_dir = output_dir,
    output_prefix = cfg$output_prefix,
    iv_thresholds = cfg$iv_thresholds,
    ...
  )
}
