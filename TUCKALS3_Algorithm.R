# Load necessary library
if (!requireNamespace("rTensor", quietly = TRUE)) {
    install.packages("rTensor")
}
library(rTensor)

# Define dimensions and initialize tensor X
I <- 5  # First dimension size
J <- 4  # Second dimension size
K <- 3  # Third dimension size
X <- array(runif(I * J * K), dim = c(I, J, K))  # Initialize X with random values

# Number of components in each mode
p <- 3
q <- 3
r <- 3

# Initialize factor matrices A, B, C with random values
set.seed(123)
A <- matrix(runif(I * p), nrow = I, ncol = p)
B <- matrix(runif(J * q), nrow = J, ncol = q)
C <- matrix(runif(K * r), nrow = K, ncol = r)

# Function to perform n-mode product
nmode_product <- function(tensor, matrix, mode) {
  # Calculate the dimension order for unfolding
  dimensions <- dim(tensor)
  mode_dim <- dimensions[mode]
  other_dims <- dimensions[-mode]
  
  # Unfold the tensor along the specified mode
  unfolded_tensor <- aperm(tensor, c(mode, setdiff(1:length(dimensions), mode)))
  unfolded_tensor <- matrix(unfolded_tensor, nrow = mode_dim, ncol = prod(other_dims))
  
  # Matrix multiplication
  if (ncol(matrix) != nrow(unfolded_tensor)) {
    stop("Non-conformable dimensions for matrix multiplication")
  }
  result_matrix <- matrix %*% unfolded_tensor
  
  # Fold the result back into a tensor
  new_dims <- c(nrow(matrix), other_dims)
  folded_tensor <- array(result_matrix, dim = new_dims)
  back_perm_order <- c(mode, setdiff(1:length(new_dims), mode))
  result_tensor <- aperm(folded_tensor, order = match(1:length(dimensions), back_perm_order))
  return(result_tensor)
}

# Tucker decomposition iterations
max_iter <- 100
epsilon <- 1e-6
for (iter in 1:max_iter) {
  # Update factor matrices A, B, C
  G <- nmode_product(X, kronecker(kronecker(C, B), A), 1)
  A_new <- qr.solve(t(G) %*% G, t(G) %*% matrix(A, ncol = p))
  G <- nmode_product(X, kronecker(kronecker(C, A_new), B), 2)
  B_new <- qr.solve(t(G) %*% G, t(G) %*% matrix(B, ncol = q))
  G <- nmode_product(X, kronecker(kronecker(B_new, A_new), C), 3)
  C_new <- qr.solve(t(G) %*% G, t(G) %*% matrix(C, ncol = r))

  # Check for convergence
  if (max(abs(A - A_new), abs(B - B_new), abs(C - C_new)) < epsilon) {
    A <- A_new
    B <- B_new
    C <- C_new
    break
  }
  A <- A_new
  B <- B_new
  C <- C_new
}

# Compute the core tensor
G_core <- nmode_product(X, kronecker(kronecker(C, B), A), 1)
G_core <- nmode_product(G_core, kronecker(kronecker(C, A), B), 2)
G_core <- nmode_product(G_core, kronecker(kronecker(B, A), C), 3)

# Output results
list(A = A, B = B, C = C, Core = G_core)
