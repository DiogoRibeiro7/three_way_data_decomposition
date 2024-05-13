#' Generate a random sample from a mixture of Gaussian distributions
#'
#' This function generates a random sample from a specified mixture of Gaussian
#' distributions. Each component of the mixture is defined by its mean vector,
#' covariance matrix, and a prior probability. The function also calculates the
#' posterior probabilities of each sample belonging to each Gaussian component.
#'
#' @param n Integer, the number of observations to generate.
#' @param p Numeric vector, prior probabilities for each Gaussian component.
#'           Must sum to 1.
#' @param M Numeric matrix, where each row represents the mean vector of a
#'           Gaussian component.
#' @param Sig List of covariance matrices, each corresponding to a Gaussian
#'           component. The list should have the same length as the number of
#'           rows in M.
#'
#' @return A list containing:
#'   - `X`: a matrix of generated data where each row corresponds to an observation.
#'   - `U`: a matrix of posterior probabilities with rows corresponding to observations
#'     and columns to Gaussian components.
#'   - `z`: an integer vector indicating the most likely component each observation
#'     was generated from.
#'
#' @examples
#' # Define parameters
#' n <- 100
#' p <- c(0.5, 0.5)
#' M <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)
#' Sig <- list(matrix(c(1, 0, 0, 1), nrow = 2), matrix(c(1, 0.5, 0.5, 2), nrow = 2))
#'
#' # Generate data
#' result <- genmixhet(n, p, M, Sig)
#' print(result$X)  # Generated data
#' print(result$U)  # Posterior probabilities
#' print(result$z)  # Component assignments
#'
#' @export
genmixhet <- function(n, p, M, Sig) {
  # n: number of observations to generate
  # p: vector of prior probabilities for each Gaussian component
  # M: matrix of means where each row corresponds to a Gaussian component
  # Sig: list of covariance matrices for each Gaussian component

  G <- nrow(M)    # Number of Gaussian components
  J <- ncol(M)    # Dimensionality of the data
  
  U <- matrix(0, n, G)
  ones_n <- rep(1, n)
  X <- matrix(0, n, J)
  
  # Generate the number of observations for each class
  cp <- cumsum(p)
  x <- runif(n)
  z <- findInterval(x, cp) + 1
  
  # Generate the data
  for (g in 1:G) {
    ind <- which(z == g)
    ng <- length(ind)
    X[ind, ] <- matrix(M[g, ], ng, J, byrow = TRUE) + matrix(rnorm(ng * J), ng, J) %*% chol(Sig[[g]])
  }
  
  # Compute the matrix U for posterior probabilities
  for (g in 1:G) {
    P <- svd(Sig[[g]])$u
    L <- svd(Sig[[g]])$d
    Q <- svd(Sig[[g]])$v
    lam <- diag(L)
    U[, g] <- -0.5 * sum(log(lam)) - 0.5 * rowSums(((X - matrix(rep(M[g, ], each = n), n, J)) %*% Q %*% diag(1/sqrt(lam)))^2)
  }
  
  # Adjust for very low probabilities
  ind <- U < -700
  U[ind] <- exp(-700)
  U[!ind] <- exp(U[!ind])
  U <- U %*% diag(p)
  U <- t(t(U) / rowSums(U))
  
  return(list(X = X, U = U, z = z))
}
