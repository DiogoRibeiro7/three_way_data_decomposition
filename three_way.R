library(rTensor)

perform_matrix_operations <- function(a, b) {
    result_ab_t <- a %*% t(b)
    result_diff1 <- c(result_ab_t) - (b %x% diag(nrow(a))) %*% a
    result_diff2 <- (b %x% diag(nrow(a))) %*% a - b %x% a

    list(
        result_ab_t = result_ab_t,
        result_diff1 = result_diff1,
        result_diff2 = result_diff2
    )
}

tensor_products <- function(i1, i2) {
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

# Simplified custom CP decomposition function
CPfunc <- function(X, R, max_iter, conv_eps) {
  # Ensure X is a tensor
  if (!inherits(X, "Tensor")) {
    X <- as.tensor(X)
  }

  # Perform CP decomposition using rTensor's cp function
  cp_result <- cp(X, num_components = R, max_iter = max_iter, tol = conv_eps)

  # Format the result
  result <- list(
    A = cp_result@U[[1]],  # First mode factor matrix
    B = cp_result@U[[2]],  # Second mode factor matrix
    C = cp_result@U[[3]],  # Third mode factor matrix
    lambda = cp_result@lambda  # Extracting lambda values
  )

  return(result)
}

# Adjusted wrapper function to match the structure of CPfunc
cp_decomposition_wrapper <- function(X, dims, max_iter, conv_eps) {
  op <- CPfunc(X, dims[[1]], conv_eps, max_iter)
  list(
    A = op$A,
    B = op$B,
    C = op$C,
    lambda = op$lambda  # Ensuring lambda is correctly passed through
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
dims <- c(2, 2, 2)
conv_eps <- 1e-6  # Convergence threshold
max_iter <- 1000  # Maximum number of iterations

# Use the wrapper function with the new CPfunc
cp_results <- cp_decomposition_wrapper(X, dims, max_iter, conv_eps)

# Print results to inspect
print(cp_results$A)
print(cp_results$B)
print(cp_results$C)
print(cp_results$lambda)  # Printing lambda values

