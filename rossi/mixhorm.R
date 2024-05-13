
mixhom <- function(X, U, eps, dis) {
  # Initialize variables
  n <- nrow(X)
  J <- ncol(X)
  G <- ncol(U)
  su <- colSums(U)
  ones_n <- matrix(1, n, 1)
  oJ <- rep(1, J)
  dif <- 1
  likeold <- -Inf
  it <- 0
  like <- 0
  
  # Main loop
  while (dif > eps) {
    it <- it + 1
    
    # Update Mmu
    Mmu <- (diag(1/su) %*% t(U)) %*% X
    
    # Update Sig
    Sig <- matrix(0, J, J)
    for (g in 1:G) {
      Xd <- X - matrix(Mmu[g, ], n, J, byrow = TRUE)
      Sig <- Sig + t(Xd) %*% diag(U[, g]) %*% Xd
    }
    Sig <- Sig / n
    svd_Sig <- svd(Sig)
    Q <- svd_Sig$u
    L <- svd_Sig$d
    lam <- diag(L)
    
    # Update U
    lfig <- matrix(0, n, G)
    for (g in 1:G) {
      lfig[, g] <- -0.5 * sum(log(lam)) - 0.5 * rowSums(((X - matrix(Mmu[g, ], n, J, byrow = TRUE)) %*% Q %*% diag(1/sqrt(lam)))^2)
    }
    
    ind <- lfig < -700
    U[ind] <- exp(-700)
    U[!ind] <- exp(lfig[!ind])
    U <- U %*% diag(su / n)
    U <- sweep(U, 1, rowSums(U), "/")
    su <- colSums(U)
    
    # Update p
    p <- su / n
    
    # Stopping rule
    ulu <- U * log(U)
    ulu[is.nan(ulu)] <- 0
    sulu <- sum(ulu)
    like <- n * sum(p * log(p)) + sum(colSums(U * lfig)) - sulu
    dif <- like - likeold
    likeold <- like
  }
  
  np <- G - 1 + G * J + (J * J + J) * 0.5
  bic <- 2 * like - log(n) * np
  aic <- 2 * like - 2 * np
  
  # Display results if dis == 1
  if (dis == 1) {
    cat(sprintf("mixhom(%g): dif=%g, iter=%g, like=%g, BIC=%g, AIC=%g\n", G, dif, it, like, bic, aic))
  }
  
  if (dif < -0.000001) {
    cat("Error in convergence\n")
  }
  
  return(list(U = U, Mmu = Mmu, Sig = Sig, dif = dif, like = like, bic = bic, it = it))
}
