#' Multilinear Principal Component Analysis
#'
#' Validate tensor inputs and delegate multilinear PCA estimation to
#' `rTensor::mpca()`. The final tensor mode is treated as the uncompressed
#' measurement mode, following the `rTensor` convention.
#'
#' @param X A numeric array or an `rTensor::Tensor` object with at least two
#'   modes.
#' @param ranks Positive integer vector giving the retained rank for each mode
#'   except the final measurement mode.
#' @param max_iter Integer greater than or equal to 2 giving the maximum number
#'   of ALS iterations.
#' @param tol Positive convergence tolerance.
#' @return The result returned by `rTensor::mpca()`, containing the extended
#'   core tensor `Z_ext`, factor matrices `U`, convergence flag `conv`,
#'   reconstruction `est`, explained norm percentage, and residual
#'   diagnostics.
#' @export
mpca_decomposition <- function(
  X,
  ranks,
  max_iter = 25L,
  tol = 1e-5
) {
  X <- as_decomposition_tensor(X)

  if (X@num_modes < 2L) {
    stop("X must have at least two modes for MPCA.")
  }

  validate_mpca_ranks(ranks, X@modes)

  if (
    !is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      !is.finite(max_iter) ||
      max_iter < 2 ||
      max_iter %% 1 != 0
  ) {
    stop("max_iter must be an integer greater than or equal to 2.")
  }

  if (
    !is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0
  ) {
    stop("tol must be a positive finite number.")
  }

  rTensor::mpca(
    X,
    ranks = as.integer(ranks),
    max_iter = as.integer(max_iter),
    tol = tol
  )
}


validate_mpca_ranks <- function(ranks, modes) {
  compressed_modes <- modes[-length(modes)]

  if (
    !is.numeric(ranks) ||
      length(ranks) != length(compressed_modes) ||
      anyNA(ranks) ||
      any(!is.finite(ranks)) ||
      any(ranks <= 0) ||
      any(ranks %% 1 != 0)
  ) {
    stop(
      "ranks must contain one positive integer for each compressed tensor mode."
    )
  }

  if (any(ranks > compressed_modes)) {
    stop("ranks cannot exceed the corresponding compressed dimensions.")
  }

  invisible(TRUE)
}
