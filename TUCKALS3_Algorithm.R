# Load necessary library
if (!requireNamespace("methods", quietly = TRUE)) {
    install.packages("methods")
}
library(methods)

# Define dimensions and initialize tensor X
I <- 5  # First dimension size
J <- 4  # Second dimension size
K <- 3  # Third dimension size
X <- array(runif(I * J * K), dim = c(I, J, K))  # Initialize X with random values

# Initialize test matrices for each mode
A_test <- matrix(runif(I * 3), nrow = I, ncol = 3)  # Test matrix for mode 1
B_test <- matrix(runif(J * 3), nrow = J, ncol = 3)  # Test matrix for mode 2, used directly in multiplication
C_test <- matrix(runif(K * 3), nrow = K, ncol = 3)  # Test matrix for mode 3, used directly in multiplication

# Compute Kronecker product correctly for modes B and C
test_matrix_A <- kronecker(B_test, C_test)  # Assume this gives us the correct dimensions

# Function to perform n-mode product
nmode_product <- function(tensor, matrix, mode) {
  dimensions <- dim(tensor)
  mode_dim <- dimensions[mode]
  other_dims <- dimensions[-mode]
  
  # Unfold tensor along the specified mode
  tensor_unfolded <- aperm(tensor, c(mode, setdiff(1:length(dimensions), mode)))
  tensor_unfolded <- matrix(tensor_unfolded, nrow = mode_dim, ncol = prod(other_dims))
  
  # Check dimensions for matrix multiplication
  if (ncol(matrix) != nrow(tensor_unfolded)) {
    stop("Non-conformable dimensions for matrix multiplication: matrix columns = ", ncol(matrix), ", tensor unfolded rows = ", nrow(tensor_unfolded))
  }
  
  # Perform multiplication
  result_matrix <- matrix %*% tensor_unfolded
  
  # Fold the result back into a tensor of the appropriate dimensions
  new_dimensions <- c(nrow(matrix), other_dims)
  result_tensor <- array(result_matrix, dim = new_dimensions)
  perm_order <- c(mode, setdiff(1:length(new_dimensions), mode))
  result_tensor <- aperm(result_tensor, order = match(1:length(dimensions), perm_order))
  
  return(result_tensor)
}

# Execute n-mode product function
test_result_A <- nmode_product(X, test_matrix_A, 1)
print("Dimensions of result after mode-1 multiplication:", dim(test_result_A))

