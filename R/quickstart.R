#' Prepare Real-Data Quickstart Assets from Packaged Zip Parts
#'
#' Reconstructs `quickstart_real_assets.zip` from split parts shipped in
#' `inst/extdata`, unzips the bundle, and returns paths needed by the core
#' MR-MOSS APIs.
#'
#' @param output_dir Directory where the quickstart bundle will be extracted.
#' @param overwrite Whether to re-create zip and re-extract when files already exist.
#' @param verbose Whether to print progress messages.
#' @return A named list with quickstart paths and trait names.
#' @export
mrmoss_prepare_quickstart_real_data <- function(output_dir = file.path(getwd(), "quickstart_real_assets"),
                                                overwrite = FALSE,
                                                verbose = TRUE) {
  part_aa <- system.file("extdata", "quickstart_real_assets.zip.part-aa", package = "MRMOSS")
  part_ab <- system.file("extdata", "quickstart_real_assets.zip.part-ab", package = "MRMOSS")

  if (!nzchar(part_aa) || !nzchar(part_ab)) {
    stop(
      paste(
        "Bundled quickstart zip parts are missing from this installation.",
        "Reinstall from GitHub and retry."
      ),
      call. = FALSE
    )
  }

  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  quickstart_root <- file.path(output_dir, "quickstart_real")
  manifest_csv <- file.path(quickstart_root, "manifest_quickstart_real.csv")
  reference_prefix <- file.path(quickstart_root, "reference", "EUR")

  if (file.exists(manifest_csv) && file.exists(paste0(reference_prefix, ".bed")) && !isTRUE(overwrite)) {
    mrmoss_message(verbose, sprintf("[quickstart] using existing assets in %s", quickstart_root))
  } else {
    zip_file <- file.path(output_dir, "quickstart_real_assets.zip")

    mrmoss_message(verbose, sprintf("[quickstart] rebuilding zip in %s", zip_file))
    mrmoss_concat_binary_files(c(part_aa, part_ab), zip_file)

    mrmoss_message(verbose, sprintf("[quickstart] extracting zip to %s", output_dir))
    utils::unzip(zip_file, exdir = output_dir, overwrite = TRUE)
  }

  if (!file.exists(manifest_csv)) {
    stop(sprintf("Quickstart manifest not found after extraction: %s", manifest_csv), call. = FALSE)
  }
  if (!file.exists(paste0(reference_prefix, ".bed"))) {
    stop(sprintf("Quickstart reference panel missing: %s.bed", reference_prefix), call. = FALSE)
  }

  list(
    quickstart_root = quickstart_root,
    manifest_csv = manifest_csv,
    reference_prefix = reference_prefix,
    exposure = "Ground_coffee_consumption",
    outcomes = c("Emotional_neglect", "Physical_abuse", "Sexual_abuse")
  )
}

mrmoss_concat_binary_files <- function(input_files, output_file, chunk_size = 1024L * 1024L) {
  out_con <- file(output_file, open = "wb")
  on.exit(close(out_con), add = TRUE)

  for (in_file in input_files) {
    in_con <- file(in_file, open = "rb")
    tryCatch(
      {
        repeat {
          buf <- readBin(in_con, what = "raw", n = chunk_size)
          if (length(buf) == 0) {
            break
          }
          writeBin(buf, out_con)
        }
      },
      finally = close(in_con)
    )
  }
}
