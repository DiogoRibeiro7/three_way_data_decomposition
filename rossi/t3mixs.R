#' t3mixs: Joint Tucker3 Model & Gaussian Mixture (ML Solution)
#'
#' This function implements a joint Tucker3 model and Gaussian mixture model using maximum likelihood (ML) estimation.
#'
#' @param X A matrix of observations, where each column is an observation.
#' @param U A matrix of initial cluster memberships.
#' @param TB A matrix of initial Tucker basis B.
#' @param TC A matrix of initial Tucker basis C.
#' @param SV An initial covariance matrix for B.
#' @param SO An initial covariance matrix for C.
#' @param eps A numeric value specifying the convergence threshold.
#' @param dis A binary value indicating whether to display the output (1) or not (0).
#' @return A list containing:
#' \describe{
#'   \item{U}{Updated matrix of cluster memberships.}
#'   \item{TB}{Updated matrix of Tucker basis B.}
#'   \item{TC}{Updated matrix of Tucker basis C.}
#'   \item{SO}{Updated covariance matrix for C.}
#'   \item{SV}{Updated covariance matrix for B.}
#'   \item{Y}{Updated matrix of factor scores.}
#'   \item{like}{Log-likelihood of the fitted model.}
#'   \item{bic}{Bayesian Information Criterion (BIC) of the fitted model.}
#' }
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' U <- matrix(runif(100 * 3), 100, 3)
#' TB <- matrix(runif(10 * 3), 10, 3)
#' TC <- matrix(runif(10 * 3), 10, 3)
#' SV <- diag(10)
#' SO <- diag(10)
#' eps <- 1e-6
#' dis <- 1
#' result <- t3mixs(X, U, TB, TC, SV, SO, eps, dis)
#' str(result)
#' }
#' @export
t3mixs <- function(X, U, TB, TC, SV, SO, eps, dis) {
  d1 <- nrow(X)
  n1 <- ncol(U)
  d2 <- nrow(TB)
  n2 <- ncol(TB)
  d3 <- nrow(TC)
  n3 <- ncol(TC)
  su <- colSums(U)
  p <- su / d1
  likeold <- -Inf
  dif <- 1
  it <- 0
  Xbar <- diag(1 / su) %*% t(U) %*% X
  M <- Xbar
  Xm <- matrix(0, d2, d3)
  Mm <- matrix(0, d2, d3)
  lfig <- matrix(0, d1, n1)
  onesd1 <- rep(1, d1)
  
  while (dif > eps) {
    it <- it + 1
    
    # update SO
    A <- 0
    SVm1 <- solve(SV)
    for (i in 1:d1) {
      Xm[] <- X[i,]
      for (g in 1:n1) {
        Mm[] <- M[g,]
        A <- A + U[i, g] * t(Xm - Mm) %*% SVm1 %*% (Xm - Mm)
      }
    }
    SO <- A / (d1 * d2)
    
    # update SV
    A <- 0
    SOm1 <- solve(SO)
    for (i in 1:d1) {
      Xm[] <- X[i,]
      for (g in 1:n1) {
        Mm[] <- M[g,]
        A <- A + U[i, g] * (Xm - Mm) %*% SOm1 %*% t(Xm - Mm)
      }
    }
    SV <- A / (d1 * d3)
    
    # update TB
    WB <- 0
    SOm1 <- solve(SO)
    TCTCp <- SOm1 %*% TC %*% t(TC) %*% SOm1
    for (g in 1:n1) {
      Xm[] <- Xbar[g,]
      WB <- WB + su[g] * Xm %*% TCTCp %*% t(Xm)
    }
    svd_SV <- svd(SV)
    l <- svd_SV$d
    SVR <- svd_SV$u %*% diag(sqrt(l)) %*% t(svd_SV$u)
    iSVR <- svd_SV$u %*% diag(1 / sqrt(l)) %*% t(svd_SV$u)
    svd_WB <- svd(iSVR %*% WB %*% iSVR)
    TB <- SVR %*% svd_WB$u[, 1:n2]
    
    # update TC
    WC <- 0
    SVm1 <- solve(SV)
    TBTBp <- SVm1 %*% TB %*% t(TB) %*% SVm1
    for (g in 1:n1) {
      Xm[] <- Xbar[g,]
      WC <- WC + su[g] * t(Xm) %*% TBTBp %*% Xm
    }
    SOR <- sqrtm(SO)
    iSOR <- solve(SOR)
    svd_WC <- svd(iSOR %*% WC %*% iSOR)
    TC <- SOR %*% svd_WC$u[, 1:n3]
    
    # update Y (eta)
    iSO <- solve(SO)
    iSV <- solve(SV)
    Y <- Xbar %*% kronecker(iSO %*% TC, iSV %*% TB)
    M <- Y %*% kronecker(t(TC), t(TB))
    
    # update U
    iSORiSVR <- kronecker(iSOR, iSVR)
    for (g in 1:n1) {
      lfig[, g] <- -0.5 * rowSums((X - onesd1 %*% t(M[g,])) %*% iSORiSVR ^ 2)
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
    like <- -0.5 * d1 * d2 * log(det(SO)) - 0.5 * d1 * d3 * log(det(SV)) + d1 * sum(p * log(p)) + sum(U * lfig) - sulu
    dif <- like - likeold
    likeold <- like
  }
  
  if (dif < -0.000001) {
    warning('------------------------------------------- error -------------------------------------------------')
    print(dif)
  }
  
  np <- n1 - 1 + d2 * d3 + (n1 - 1) * n2 * n3 + (d2 - n2) * n2 + (d3 - n3) * n3 + (d2 * d2 + d2) / 2 + (d3 * d3 + d3) / 2 - 1
  bic <- 2 * like - log(d1) * np
  if (dis == 1) {
    message(sprintf('T3mix: dif=%g, iter=%g, like=%g, np=%g, BIC=%g', dif, it, like, np, bic))
  }
  
  list(U = U, TB = TB, TC = TC, SO = SO, SV = SV, Y = Y, like = like, bic = bic)
}

# Example usage:
# X <- matrix(rnorm(100 * 10), 100, 10)
# U <- matrix(runif(100 * 3), 100, 3)
# TB <- matrix(runif(10 * 3), 10, 3)
# TC <- matrix(runif(10 * 3), 10, 3)
# SV <- diag(10)
# SO <- diag(10)
# eps <- 1e-6
# dis <- 1
# result <- t3mixs(X, U, TB, TC, SV, SO, eps, dis)
# str(result)
