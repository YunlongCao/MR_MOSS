#' Low-Level C++ Backend Interfaces
#'
#' Advanced interfaces to the internal C++ routines used by MR-MOSS.
#'
#' @param gamma_hat Numeric vector of exposure effect estimates.
#' @param Gamma_hat Numeric matrix of outcome effect estimates.
#' @param R Outcome correlation matrix.
#' @param n1 Exposure sample size.
#' @param n2 Outcome sample size.
#' @param theta Numeric parameter vector for log-likelihood evaluation.
#' @param theta0 Numeric initial parameter vector for PX-EM fitting.
#' @param test Integer vector of outcome indices to test.
#' @param maxiter Maximum number of PX-EM iterations.
#' @param rd Multiplicative update factor for residual standard deviations.
#'
#' @return
#' `loglikelihood_cpp()` returns a numeric scalar log-likelihood value.
#'
#' `MRMOSS_PX_cpp()` returns a list with model estimates, p-values, and
#' optimization diagnostics.
#'
#' @name mrmoss_cpp_backend
#' @aliases loglikelihood_cpp MRMOSS_PX_cpp
NULL
