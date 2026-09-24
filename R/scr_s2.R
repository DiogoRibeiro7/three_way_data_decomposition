#' Fit the Two-Way SCR Baseline
#'
#' Fit model S2 from Rocci, Vichi, and Ranalli: simultaneous clustering and
#' dimensionality reduction after vectorising the three-way observations.
#'
#' The implementation follows the authors' MATLAB reference algorithm. It
#' alternates updates of the common covariance, discriminating basis, latent
#' group scores, posterior memberships, and mixing proportions.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Numeric initial posterior-membership matrix with
#'   observations in rows and mixture components in columns.
#' @param basis Numeric initial reduction basis. It must have `ncol(X)` rows;
#'   its number of columns is the retained S2 dimension.
#' @param tolerance Positive convergence threshold for the reference
#'   log-likelihood increment.
#' @param max_iter Positive integer maximum number of iterations.
#' @param display Logical; print a compact convergence summary when `TRUE`.
#' @return A list containing posterior memberships `U`, reduction basis
#'   `TB`, common covariance `SV`, latent group scores `Y`, fitted group
#'   means `M`, final likelihood increment `dif`, reference log-likelihood
#'   `like`, reference BIC score `bic`, iteration count `it`, mixing
#'   proportions `probabilities`, convergence flag `converged`, and the
#'   likelihood trajectory `log_likelihood_trace`.
#' @export
fit_scr_s2 <- function(
  X,
  membership,
  basis,
  tolerance = 1e-6,
  max_iter = 1000L,
  display = FALSE
) {
  validate_scr_s2_inputs(
    X = X,
    membership = membership,
    basis = basis,
    tolerance = tolerance,
    max_iter = max_iter,
    display = display
  )

  n <- nrow(X)
  dimension <- ncol(X)
  groups <- ncol(membership)
  rank <- ncol(basis)

  U <- membership
  TB <- basis
  SV <- diag(dimension)

  group_mass <- colSums(U)
  if (any(group_mass <= 0)) {
    stop("each component must have positive initial membership mass.")
  }

  probabilities <- group_mass / n
  Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")
  Y <- Xbar %*% TB
  M <- Y %*% t(TB)

  previous_likelihood <- -Inf
  difference <- Inf
  iteration <- 0L
  likelihood <- -Inf
  likelihood_trace <- numeric(0)
  converged <- FALSE

  while (iteration < max_iter && difference > tolerance) {
    iteration <- iteration + 1L

    # Common covariance update.
    SV <- matrix(0, nrow = dimension, ncol = dimension)
    for (g in seq_len(groups)) {
      centered <- sweep(X, 2L, M[g, ], "-")
      weighted <- centered * sqrt(U[, g])
      SV <- SV + crossprod(weighted)
    }
    SV <- SV / n
    SV <- (SV + t(SV)) / 2

    covariance_svd <- svd(SV)
    eigenvalues <- covariance_svd$d

    if (
      any(!is.finite(eigenvalues)) ||
        any(eigenvalues <= sqrt(.Machine$double.eps))
    ) {
      stop("S2 covariance became singular during fitting.")
    }

    P <- covariance_svd$u
    covariance_root <- P %*%
      diag(sqrt(eigenvalues), nrow = dimension) %*%
      t(P)
    inverse_covariance_root <- P %*%
      diag(1 / sqrt(eigenvalues), nrow = dimension) %*%
      t(P)

    # Reduced discriminating subspace update.
    weighted_centroids <- Xbar * sqrt(group_mass)
    WB <- crossprod(weighted_centroids)
    reduced_svd <- svd(
      inverse_covariance_root %*% WB %*% inverse_covariance_root
    )
    TB <- covariance_root %*%
      reduced_svd$u[, seq_len(rank), drop = FALSE]

    # Latent group scores and reconstructed group means.
    Y <- Xbar %*% solve(SV, TB)
    M <- Y %*% t(TB)

    # Posterior membership update. The transformed residual must be squared
    # after both matrix products, matching the MATLAB reference expression.
    log_kernel <- matrix(0, nrow = n, ncol = groups)
    whitening <- P %*%
      diag(1 / sqrt(eigenvalues), nrow = dimension)

    for (g in seq_len(groups)) {
      centered <- sweep(X, 2L, M[g, ], "-")
      transformed <- centered %*% whitening
      log_kernel[, g] <- -0.5 * rowSums(transformed^2)
    }

    kernel <- exp(pmax(log_kernel, -700))
    weighted_kernel <- sweep(
      kernel,
      MARGIN = 2L,
      STATS = probabilities,
      FUN = "*"
    )
    normalizer <- rowSums(weighted_kernel)

    if (any(!is.finite(normalizer)) || any(normalizer <= 0)) {
      stop("S2 posterior normalization failed during fitting.")
    }

    U <- weighted_kernel / normalizer
    group_mass <- colSums(U)

    if (any(group_mass <= 0)) {
      stop("an S2 component acquired zero membership mass during fitting.")
    }

    Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")
    probabilities <- group_mass / n

    entropy_terms <- U * log(U)
    entropy_terms[is.nan(entropy_terms)] <- 0
    entropy <- sum(entropy_terms)

    log_determinant <- as.numeric(
      determinant(SV, logarithm = TRUE)$modulus
    )

    likelihood <- -0.5 * n * log_determinant +
      n * sum(probabilities * log(probabilities)) +
      sum(U * log_kernel) -
      entropy

    difference <- likelihood - previous_likelihood
    likelihood_trace <- c(likelihood_trace, likelihood)
    previous_likelihood <- likelihood

    if (is.finite(difference) && difference <= tolerance) {
      converged <- difference >= -tolerance
      break
    }
  }

  n_parameters <- groups - 1L +
    dimension +
    (groups - 1L) * rank +
    (dimension - rank) * rank +
    (dimension * dimension + dimension) / 2 -
    1L

  bic <- 2 * likelihood - log(n) * n_parameters

  if (difference < -1e-6) {
    warning(
      "S2 reference likelihood decreased before termination.",
      call. = FALSE
    )
  }

  if (display) {
    message(
      sprintf(
        "S2: dif=%g, iter=%d, like=%g, np=%g, BIC=%g",
        difference,
        iteration,
        likelihood,
        n_parameters,
        bic
      )
    )
  }

  list(
    U = U,
    TB = TB,
    SV = SV,
    Y = Y,
    M = M,
    dif = difference,
    like = likelihood,
    bic = bic,
    it = iteration,
    probabilities = probabilities,
    converged = converged,
    log_likelihood_trace = likelihood_trace
  )
}


#' Legacy SCR Name for the Two-Way Model
#'
#' Compatibility wrapper for the original MATLAB/R `t2mixt()` interface.
#'
#' @param X See `fit_scr_s2()`.
#' @param U Initial posterior-membership matrix.
#' @param TB Initial reduction basis.
#' @param eps Convergence threshold.
#' @param dis Numeric/logical display flag retained for MATLAB compatibility.
#' @param max_iter Maximum number of iterations.
#' @return A list containing the historical six fields `U`, `TB`, `SV`,
#'   `Y`, `like`, and `bic`.
#' @export
t2mixt <- function(
  X,
  U,
  TB,
  eps = 1e-6,
  dis = 0,
  max_iter = 1000L
) {
  fit <- fit_scr_s2(
    X = X,
    membership = U,
    basis = TB,
    tolerance = eps,
    max_iter = max_iter,
    display = isTRUE(as.logical(dis))
  )

  fit[c("U", "TB", "SV", "Y", "like", "bic")]
}


validate_scr_s2_inputs <- function(
  X,
  membership,
  basis,
  tolerance,
  max_iter,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 2L ||
      ncol(X) < 2L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop(
      "X must be a finite numeric matrix with at least two observations and two variables."
    )
  }

  if (
    !is.matrix(membership) ||
      !is.numeric(membership) ||
      nrow(membership) != nrow(X) ||
      ncol(membership) < 2L ||
      anyNA(membership) ||
      any(!is.finite(membership)) ||
      any(membership < 0)
  ) {
    stop(
      "membership must be a finite non-negative numeric matrix with one row per observation and at least two components."
    )
  }

  row_totals <- rowSums(membership)
  if (any(row_totals <= 0)) {
    stop("each membership row must have positive mass.")
  }
  if (max(abs(row_totals - 1)) > 1e-10) {
    stop("each membership row must sum to one.")
  }

  if (
    !is.matrix(basis) ||
      !is.numeric(basis) ||
      nrow(basis) != ncol(X) ||
      ncol(basis) < 1L ||
      ncol(basis) > ncol(X) ||
      anyNA(basis) ||
      any(!is.finite(basis))
  ) {
    stop(
      "basis must be a finite numeric matrix with ncol(X) rows and between 1 and ncol(X) columns."
    )
  }

  if (qr(basis)$rank < ncol(basis)) {
    stop("basis must have full column rank.")
  }

  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      !is.finite(tolerance) ||
      tolerance <= 0
  ) {
    stop("tolerance must be a positive finite number.")
  }

  if (
    !is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      !is.finite(max_iter) ||
      max_iter <= 0 ||
      max_iter %% 1 != 0
  ) {
    stop("max_iter must be a positive integer.")
  }

  if (!is.logical(display) || length(display) != 1L || is.na(display)) {
    stop("display must be TRUE or FALSE.")
  }

  invisible(TRUE)
}
