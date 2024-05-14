#' t2mixt: Joint Tucker3 Model & Gaussian Mixture (ML Solution)
#'
#' This function implements a joint Tucker3 model and Gaussian mixture model using maximum likelihood (ML) estimation.
#'
#' @param X A matrix of observations, where each column is an observation.
#' @param U A matrix of initial cluster memberships.
#' @param TB A matrix of initial Tucker basis.
#' @param eps A numeric value specifying the convergence threshold.
#' @param dis A binary value indicating whether to display the output (1) or not (0).
#' @return A list containing:
#' \describe{
#'   \item{U}{Updated matrix of cluster memberships.}
#'   \item{TB}{Updated matrix of Tucker basis.}
#'   \item{SV}{Estimated covariance matrix.}
#'   \item{Y}{Updated matrix of factor scores.}
#'   \item{like}{Log-likelihood of the fitted model.}
#'   \item{bic}{Bayesian Information Criterion (BIC) of the fitted model.}
#' }
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' U <- matrix(runif(100 * 3), 100, 3)
#' TB <- matrix(runif(10 * 3), 10, 3)
#' eps <- 1e-6
#' dis <- 1
#' result <- t2mixt(X, U, TB, eps, dis)
#' str(result)
#' }
#' @export
t2mixt <- function(X, U, TB, eps, dis) {
  # joint Tucker3 model & Gaussian mixture (ML solution)
  # X = [X_1 ...X_K]
  # X  = U*G*kron(C,B)'
  # Xy = B*Gy*kron(C,U)'
  # Xz = C*Gz*kron(B,U)'
  # 21/09/2016
  
  d1 <- nrow(X)
  n1 <- ncol(U)
  d2 <- nrow(TB)
  n2 <- ncol(TB)
  SV <- diag(d2)
  su <- colSums(U)
  p <- su / d1
  likeold <- -Inf
  dif <- 1
  it <- 0
  Xbar <- diag(1 / su) %*% t(U) %*% X
  Y <- Xbar %*% TB
  M <- Y %*% t(TB)
  onesd1 <- rep(1, d1)
  od2 <- rep(1, d2)
  
  while (dif > eps) {
    it <- it + 1
    
    # update SV
    SV <- matrix(0, d2, d2)
    for (g in 1:n1) {
      Xd <- X - onesd1 %*% t(M[g,])
      SV <- SV + t(Xd) %*% diag(U[, g]) %*% Xd
    }
    SV <- SV / d1
    
    # update TB
    WB <- t(Xbar) %*% diag(su) %*% Xbar
    svd_SV <- svd(SV)
    l <- svd_SV$d
    SVR <- svd_SV$u %*% diag(sqrt(l)) %*% t(svd_SV$u)
    iSVR <- svd_SV$u %*% diag(1 / sqrt(l)) %*% t(svd_SV$u)
    svd_WB <- svd(iSVR %*% WB %*% iSVR)
    TB <- SVR %*% svd_WB$u[, 1:n2]
    
    # update Y
    STm1 <- solve(SV)
    Y <- Xbar %*% STm1 %*% TB
    M <- Y %*% t(TB)
    
    # update U
    lfig <- matrix(0, d1, n1)
    for (g in 1:n1) {
      lfig[, g] <- -0.5 * rowSums((X - onesd1 %*% t(M[g,])) %*% svd_SV$u %*% diag(1 / sqrt(l)) ^ 2)
    }
    ind <- (lfig < -7e2)
    U <- exp(-7e2 * ind + lfig * (1 - ind)) %*% diag(p)
    U <- diag(1 / rowSums(U)) %*% U
    su <- colSums(U)
    Xbar <- diag(1 / su) %*% t(U) %*% X
    
    # update p
    p <- su / d1
    
    # stopping rule
    suppressWarnings({
      UlU <- U * log(U)
      UlU[is.nan(UlU)] <- 0
      sulu <- sum(UlU)
    })
    like <- -0.5 * d1 * log(det(SV)) + d1 * sum(p * log(p)) + sum(U * lfig) - sulu
    dif <- like - likeold
    likeold <- like
  }
  
  if (dif < -0.000001) {
    warning('-------------------------------------- error t2mixt --------------------------------------------')
    print(dif)
    print(svd(SV)$d)
  }
  
  np <- n1 - 1 + d2 + (n1 - 1) * n2 + (d2 - n2) * n2 + (d2 * d2 + d2) / 2 - 1
  bic <- 2 * like - log(d1) * np
  if (dis == 1) {
    message(sprintf('T2mix: dif=%g, iter=%g, like=%g, np=%g, BIC=%g', dif, it, like, np, bic))
  }
  
  list(U = U, TB = TB, SV = SV, Y = Y, like = like, bic = bic)
}

# Example usage:
# X <- matrix(rnorm(100 * 10), 100, 10)
# U <- matrix(runif(100 * 3), 100, 3)
# TB <- matrix(runif(10 * 3), 10, 3)
# eps <- 1e-6
# dis <- 1
# result <- t2mixt(X, U, TB, eps, dis)
# str(result)
