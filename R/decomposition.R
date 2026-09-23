#' Higher-Order Singular Value Decomposition
#'
#' Validate a tensor and delegate HOSVD estimation to `rTensor::hosvd()`.
#'
#' @param X A numeric array or an `rTensor::Tensor` object.
#' @param ranks Optional positive integer vector giving the retained rank for
#'   each tensor mode. When `NULL`, no truncation is requested.
#' @return The decomposition returned by `rTensor::hosvd()`, containing the
#'   core tensor `Z`, factor matrices `U`, reconstructed tensor `est`,
#'   and Frobenius residual `fnorm_resid`.
#' @export
hosvd_decomposition <- function(X, ranks = NULL) {
  X <- as_decomposition_tensor(X)

  if (!is.null(ranks)) {
    validate_decomposition_ranks(ranks, X@modes)
  }

  rTensor::hosvd(X, ranks = ranks)
}


#' Tucker Decomposition
#'
#' Validate a tensor and delegate Tucker decomposition to `rTensor::tucker()`.
#'
#' @param X A numeric array or an `rTensor::Tensor` object.
#' @param ranks Positive integer vector giving the retained rank for each tensor
#'   mode.
#' @param max_iter Integer greater than or equal to 2 giving the maximum number
#'   of ALS iterations.
#' @param tol Positive convergence tolerance.
#' @return The decomposition returned by `rTensor::tucker()`, including the
#'   core tensor `Z`, factor matrices `U`, convergence flag `conv`,
#'   reconstructed tensor `est`, explained norm percentage, and residual
#'   diagnostics.
#' @export
tucker_decomposition <- function(
  X,
  ranks,
  max_iter = 25L,
  tol = 1e-5
) {
  X <- as_decomposition_tensor(X)
  validate_decomposition_ranks(ranks, X@modes)

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

  rTensor::tucker(
    X,
    ranks = as.integer(ranks),
    max_iter = as.integer(max_iter),
    tol = tol
  )
}


as_decomposition_tensor <- function(X) {
  if (inherits(X, "Tensor")) {
    return(X)
  }

  if (!is.array(X) || !is.numeric(X)) {
    stop("X must be a numeric array or an rTensor Tensor.")
  }

  if (length(dim(X)) < 2L) {
    stop("X must have at least two modes.")
  }

  rTensor::as.tensor(X)
}


validate_decomposition_ranks <- function(ranks, modes) {
  if (
    !is.numeric(ranks) ||
      length(ranks) != length(modes) ||
      anyNA(ranks) ||
      any(!is.finite(ranks)) ||
      any(ranks <= 0) ||
      any(ranks %% 1 != 0)
  ) {
    stop("ranks must contain one positive integer for each tensor mode.")
  }

  if (any(ranks > modes)) {
    stop("ranks cannot exceed the corresponding tensor dimensions.")
  }

  invisible(TRUE)
}
