# Manual examples for three_way.R
#
# These examples are intentionally kept outside the reusable implementation so
# sourcing three_way.R does not execute simulations, print output, or mutate the
# global workspace.

source("three_way.R")

set.seed(123)

# Matrix operations
a <- matrix(runif(10), nrow = 10, ncol = 1)
b <- matrix(runif(5), nrow = 5, ncol = 1)

result <- perform_matrix_operations(a, b)
print(result$result_ab_t)
print(result$result_diff1)
print(result$result_diff2)

# Tensor products
i1 <- matrix(c(1, 0), nrow = 2, ncol = 1)
i2 <- matrix(c(0, 1), nrow = 2, ncol = 1)

tensors <- tensor_products(i1, i2)
print(tensors$i1xi1xi1)

# Structured matrix/tensor operations
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

complex_results <- complex_matrix_operations(A, B, C, G, I, J, K, P, Q, R)
print(complex_results$X)
print(complex_results$row_col_diff)
print(complex_results$seq_diff)

# CP decomposition
#
# Keep this example explicit about the required Tensor input and scalar rank.
X <- rTensor::as.tensor(array(rnorm(2 * 3 * 4), dim = c(2, 3, 4)))
cp_results <- cp_decomposition_wrapper(
  X,
  dims = 2,
  max_iter = 1000,
  conv_eps = 1e-6
)

print(cp_results$A)
print(cp_results$B)
print(cp_results$C)
print(cp_results$lambda)
