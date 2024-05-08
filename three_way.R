library(rTensor)

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



complex_matrix_operations <- function(A, B, C, G, I, J, K, P, Q, R) {
    X <- A %*% G %*% t(C %x% B)
    row_col_diff <- X[2, 1:J] - B %*% t(C[1, 1] * G[, 1:Q] + C[1, 2] * G[, Q + (1:Q)]) %*% A[2, ]
    seq_diff <- X[1, seq(1, J * K, J)] - C %*% t(B[1, 1] * G[, seq(1, Q * R, Q)] + B[1, 2] * G[, seq(2, Q * R, Q)] + B[1, 3] * G[, seq(3, Q * R, Q)]) %*% A[1, ]

    list(
        X = X,
        row_col_diff = row_col_diff,
        seq_diff = seq_diff
    )
}

# Adjusted CPfunc to handle outputs as a list
CPfunc <- function(X, R, max_iter, conv_eps) {
    # Ensure X is a tensor
    if (!inherits(X, "Tensor")) {
        X <- as.tensor(X)
    }

    # Perform CP decomposition using rTensor's cp function
    cp_result <- cp(X, num_components = R, max_iter = max_iter, tol = conv_eps)

    # Check if cp_result is a list and handle accordingly
    if (is.list(cp_result)) {
        # Access the factor matrices and lambda values assuming they are stored in a list
        A = cp_result$U[[1]]
        B = cp_result$U[[2]]
        C = cp_result$U[[3]]
        lambda = cp_result$lambda
    } else {
        # Fallback if cp_result is an S4 object and the original method works
        A = cp_result@U[[1]]
        B = cp_result@U[[2]]
        C = cp_result@U[[3]]
        lambda = cp_result@lambda
    }

    # Return the formatted result
    result <- list(
        A = A,
        B = B,
        C = C,
        lambda = lambda
    )

    return(result)
}

# Wrapper function to handle the decomposition
cp_decomposition_wrapper <- function(X, dims, max_iter, conv_eps) {
    op <- CPfunc(X, dims[[1]], max_iter, conv_eps)
    list(
        A = op$A,
        B = op$B,
        C = op$C,
        lambda = op$lambda
    )
}
# Test perform_matrix_operations function
a <- matrix(runif(10), nrow = 10, ncol = 1)
b <- matrix(runif(5), nrow = 5, ncol = 1)

# Run the function
result <- perform_matrix_operations(a, b)

# Print results to manually check or use assertive checks
print(result$result_ab_t) # Should be a 10x5 matrix
print(result$result_diff1) # Check specific properties that should be mathematically true
print(result$result_diff2)


# Test tensor_products function
i1 <- matrix(c(1, 0), nrow = 2, ncol = 1)
i2 <- matrix(c(0, 1), nrow = 2, ncol = 1)

# Run the function
tensors <- tensor_products(i1, i2)

# Print and inspect the outputs
print(tensors$i1xi1xi1) # Expect a 2x2x2 tensor of specific format


I <- 10
J <- 4
K <- 3
P <- 4
Q <- 3
R <- 2
A <- matrix(rnorm(I * P), nrow = I, ncol = P)
B <- matrix(rnorm(J * Q), nrow = J, ncol = Q)
C <- matrix(rnorm(K * R), nrow = K, ncol = R)
G <- matrix(rnorm(P * Q * R), nrow = P, ncol = Q * R)

# Run the function
complex_results <- complex_matrix_operations(A, B, C, G, I, J, K, P, Q, R)

# Check each component
print(complex_results$X)
print(complex_results$row_col_diff)
print(complex_results$seq_diff)


# Example usage
set.seed(123)
X <- array(rnorm(2*3*4), dim = c(2, 3, 4))
dims = c(2, 2, 2)
conv_eps = 1e-6  # Convergence threshold
max_iter = 1000  # Maximum number of iterations

# Use the wrapper function with the new CPfunc
cp_results <- cp_decomposition_wrapper(X, dims, max_iter, conv_eps)

# Print results to inspect
print(cp_results$A)
print(cp_results$B)
print(cp_results$C)
print(cp_results$lambda)



# Assuming 'X' is a tensor you want to decompose
X <- as.tensor(array(runif(24), dim = c(3, 4, 2)))

# Try-Catch to handle any potential errors gracefully
tryCatch({
    cp_result <- cp(X, num_components = 2)
    # Access components safely
    if ("lambda" %in% names(cp_result)) {
        lambdas <- cp_result$lambda
        print(lambdas)
    } else {
        print("Lambda values are not available in the result.")
    }
}, error = function(e) {
    print(paste("An error occurred:", e$message))
})