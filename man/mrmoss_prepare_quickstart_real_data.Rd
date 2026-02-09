% Generated manually for package documentation
\name{mrmoss_prepare_quickstart_real_data}
\alias{mrmoss_prepare_quickstart_real_data}
\title{Prepare Real-Data Quickstart Assets from Packaged Zip Parts}
\usage{
mrmoss_prepare_quickstart_real_data(
  output_dir = file.path(getwd(), "quickstart_real_assets"),
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{output_dir}{Directory where the quickstart bundle will be extracted.}

\item{overwrite}{Whether to re-create zip and re-extract when files already exist.}

\item{verbose}{Whether to print progress messages.}
}
\value{
A named list with quickstart paths and trait names.
}
\description{
Reconstructs \code{quickstart_real_assets.zip} from split parts shipped in
\code{inst/extdata}, unzips the bundle, and returns paths needed by the core
MR-MOSS APIs.
}
