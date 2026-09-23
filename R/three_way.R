#' Perform Matrix Operations
#'
#' This function performs specific matrix operations including multiplication of a matrix and
#' the transpose of another and calculates differences based on tensor multiplication with a diagonal matrix.
#'
#' @param a A numeric matrix.
#' @param b Another numeric matrix.
#' @return A list containing:
#'   - `result_ab_t`: the product of `a` and the transpose of `b`.
#'   - `result_diff1`: the difference between the vectorized product of `a` and `b^T` and the product of `b %x% diag(nrow(a))` and `a`.
#'   - `result_diff2`: the difference between the product of `b %x% diag(nrow(a))` and `a` and the tensor product of `b` and `a`.
#' @examples
#' a <- matrix(1:10, ncol = 1)
#' b <- matrix(1:5, ncol = 1)
#' perform_matrix_operations(a, b)
#' @export
perform_matrix_operations <- function(a, b) {
  if (!is.matrix(a) || !is.matrix(b)) {
    stop("Both arguments must be matrices.")
  }

  result_ab_t <- a %*% t(b)
  result_diff1 <- c(result_ab_t) - (b %x% diag(nrow(a))) %*% a
  result_diff2 <- (b %x% diag(nrow(a))) %*% a - b %x% a
  
  list(
    result_ab_t = result_ab_t,
    result_diff1 = result_diff1,
    result_diff2 = result_diff2
  )
}


#' Compute Various Tensor Products of Two Vectors
#'
#' This function computes all possible tensor products of two input vectors i1 and i2,
#' considering combinations of tensoring them up to three times in all possible orders.
#'
#' @param i1 A numeric vector or matrix to be used in tensor products.
#' @param i2 A numeric vector or matrix to be used in tensor products.
#' @return A list containing various tensor products:
#'   - `i1xi1xi1`: Tensor product of i1 with itself three times.
#'   - `i1xi1xi2`: Tensor product of i1 with itself twice and then with i2.
#'   - `i1xi2xi1`: Tensor product of i1 with i2 and then with i1.
#'   - `i1xi2xi2`: Tensor product of i1 with i2 twice.
#'   - `i2xi1xi1`: Tensor product of i2 with i1 twice.
#'   - `i2xi1xi2`: Tensor product of i2 with i1 and then with i2.
#'   - `i2xi2xi1`: Tensor product of i2 with itself and then with i1.
#'   - `i2xi2xi2`: Tensor product of i2 with itself three times.
#' @examples
#' i1 <- matrix(c(1, 0), ncol = 1)
#' i2 <- matrix(c(0, 1), ncol = 1)
#' tensor_products(i1, i2)
#' @export
tensor_products <- function(i1, i2) {
  if (!is.numeric(i1) || !is.numeric(i2)) {
    stop("Both i1 and i2 must be numeric vectors or matrices.")
  }

  list(
    i1xi1xi1 = i1 %x% i1 %x% i1,
    i1xi1xi2 = i1 %x% i1 %x% i2,
    i1xi2xi1 = i1 %x% i2 %x% i1,
    i1xi2xi2 = i1 %x% i2 %x% i2,
    i2xi1xi1 = i2 %x% i1 %x% i1,
    i2xi1xi2 = i2 %x% i1 %x% i2,
    i2xi2xi1 = i2 %x% i2 %x% i1,
    i2xi2xi2 = i2 %x% i2 %x% i2
  )
}


#' Perform Complex Matrix Operations Involving Tensor Products
#'
#' This function computes complex matrix operations involving multiple matrix products,
#' tensor products, and different kinds of matrix transformations. Specifically, it calculates
#' a product involving matrices A, B, C, and G, then computes differences based on subsetting
#' and matrix transformations to explore specific patterns in the results.
#'
#' @param A Matrix of dimensions I x P.
#' @param B Matrix of dimensions J x Q.
#' @param C Matrix of dimensions K x R.
#' @param G Matrix of dimensions P x (Q * R).
#' @param I Number of rows in matrix A.
#' @param J Number of columns in matrix B relevant for subsetting.
#' @param K Number of steps in sequence generation.
#' @param P Number of columns in matrix A and rows in matrix G.
#' @param Q Number of columns in matrix B and part of the dimensions for G.
#' @param R Number of columns in matrix C and part of the dimensions for G.
#' @return A list containing:
#'   - `X`: Result of the matrix product A %*% G %*% t(C %x% B).
#'   - `row_col_diff`: Difference computation for the second row across the first J columns.
#'   - `seq_diff`: Difference computation along a sequence derived from matrix dimensions.
#' @examples
#' I <- 10; J <- 4; K <- 3
#' P <- 4; Q <- 3; R <- 2
#' A <- matrix(rnorm(I * P), I, P)
#' B <- matrix(rnorm(J * Q), J, Q)
#' C <- matrix(rnorm(K * R), K, R)
#' G <- matrix(rnorm(P * Q * R), P, Q * R)
#' results <- complex_matrix_operations(A, B, C, G, I, J, K, P, Q, R)
#' @export
complex_matrix_operations <- function(A, B, C, G, I, J, K, P, Q, R) {
    if (!all(dim(A) == c(I, P), dim(B) == c(J, Q), dim(C) == c(K, R), dim(G) == c(P, Q * R))) {
        stop("Input matrices do not match specified dimensions.")
    }

    # Compute the product of matrix A, G, and the transpose of the tensor product of C and B
    X <- A %*% G %*% t(C %x% B)

    # Difference calculation for the second row across the first J columns
    row_col_diff <- X[2, 1:J] - B %*% t(C[1, 1] * G[, 1:Q] + C[1, 2] * G[, Q + (1:Q)]) %*% A[2, ]

    # Difference calculation along a sequence derived from matrix dimensions
    seq_diff <- X[1, seq(1, J * K, J)] - C %*% t(B[1, 1] * G[, seq(1, Q * R, Q)] +
                                                     B[1, 2] * G[, seq(2, Q * R, Q)] +
                                                     B[1, 3] * G[, seq(3, Q * R, Q)]) %*% A[1, ]

    list(
        X = X,
        row_col_diff = row_col_diff,
        seq_diff = seq_diff
    )
}

#' Perform CANDECOMP/PARAFAC (CP) Decomposition
#'
#' Adapt a three-mode tensor to the CP decomposition implemented by
#' `rTensor::cp()` and return the three factor matrices with their component
#' weights.
#'
#' @param X A numeric three-dimensional array or an `rTensor::Tensor` object.
#' @param R Positive integer giving the number of CP components.
#' @param max_iter Integer greater than or equal to 2 giving the maximum number
#'   of ALS iterations.
#' @param conv_eps Positive convergence tolerance passed to `rTensor::cp()`.
#' @return A list containing:
#'   - `A`: factor matrix for the first mode.
#'   - `B`: factor matrix for the second mode.
#'   - `C`: factor matrix for the third mode.
#'   - `lambda`: component weights returned by `rTensor::cp()`.
#' @examples
#' set.seed(123)
#' X <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
#' result <- CPfunc(X, R = 2, max_iter = 25, conv_eps = 1e-5)
#' @export
CPfunc <- function(X, R, max_iter = 25L, conv_eps = 1e-5) {
  if (!inherits(X, "Tensor")) {
    if (!is.array(X) || !is.numeric(X)) {
      stop("X must be a numeric array or an rTensor Tensor.")
    }
    X <- rTensor::as.tensor(X)
  }

  if (X@num_modes != 3L) {
    stop("CPfunc currently supports only three-mode tensors.")
  }
  if (!is.numeric(X@data)) {
    stop("X must contain numeric data.")
  }

  if (
    !is.numeric(R) ||
      length(R) != 1L ||
      is.na(R) ||
      !is.finite(R) ||
      R <= 0 ||
      R %% 1 != 0
  ) {
    stop("R must be a positive integer.")
  }

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
    !is.numeric(conv_eps) ||
      length(conv_eps) != 1L ||
      is.na(conv_eps) ||
      !is.finite(conv_eps) ||
      conv_eps <= 0
  ) {
    stop("conv_eps must be a positive finite number.")
  }

  cp_result <- rTensor::cp(
    X,
    num_components = as.integer(R),
    max_iter = as.integer(max_iter),
    tol = conv_eps
  )

  if (
    !is.list(cp_result) ||
      !all(c("U", "lambdas") %in% names(cp_result)) ||
      length(cp_result$U) != 3L
  ) {
    stop("Unexpected result returned by rTensor::cp().")
  }

  list(
    A = cp_result$U[[1L]],
    B = cp_result$U[[2L]],
    C = cp_result$U[[3L]],
    lambda = cp_result$lambdas
  )
}


#' Wrapper for CP Decomposition
#'
#' Convenience wrapper around `CPfunc()` retaining the historical `dims`
#' argument name for the CP rank.
#'
#' @param X A numeric three-dimensional array or an `rTensor::Tensor` object.
#' @param dims Positive scalar giving the CP rank.
#' @param max_iter Integer greater than or equal to 2 giving the maximum number
#'   of ALS iterations.
#' @param conv_eps Positive convergence tolerance.
#' @return The same four-element list returned by `CPfunc()`.
#' @examples
#' set.seed(123)
#' X <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
#' result <- cp_decomposition_wrapper(
#'   X,
#'   dims = 2,
#'   max_iter = 25,
#'   conv_eps = 1e-5
#' )
#' @export
cp_decomposition_wrapper <- function(
  X,
  dims,
  max_iter = 25L,
  conv_eps = 1e-5
) {
  CPfunc(
    X = X,
    R = dims,
    max_iter = max_iter,
    conv_eps = conv_eps
  )
}
