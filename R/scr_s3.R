#' Fit the Three-Way SCR Baseline
#'
#' Fit model S3 from Rocci, Vichi, and Ranalli: simultaneous clustering and
#' dimensionality reduction for three-way data with a Tucker2 mean structure
#' and separable common covariance.
#'
#' The observations are supplied in vectorised form, with `nrow(variable_basis)`
#' variables and `nrow(occasion_basis)` occasions. The model uses
#' `kronecker(occasion_basis, variable_basis)` for the discriminating mean
#' subspace and `kronecker(occasion_covariance, variable_covariance)` for the
#' common covariance structure.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Numeric initial posterior-membership matrix.
#' @param variable_basis Numeric initial variable-mode basis `B`.
#' @param occasion_basis Numeric initial occasion-mode basis `C`.
#' @param variable_covariance Positive-definite initial covariance `Sigma_V`.
#' @param occasion_covariance Positive-definite initial covariance `Sigma_O`.
#' @param tolerance Positive convergence threshold for the reference
#'   log-likelihood increment.
#' @param max_iter Positive integer maximum number of iterations.
#' @param display Logical; print a compact convergence summary when `TRUE`.
#' @return A list containing posterior memberships `U`, variable basis `TB`,
#'   occasion basis `TC`, occasion covariance `SO`, variable covariance
#'   `SV`, latent group scores `Y`, reconstructed group means `M`, final
#'   likelihood increment `dif`, reference log-likelihood `like`, reference
#'   BIC score `bic`, iteration count `it`, mixing proportions
#'   `probabilities`, convergence flag `converged`, and likelihood trajectory
#'   `log_likelihood_trace`.
#' @export
fit_scr_s3 <- function(
  X,
  membership,
  variable_basis,
  occasion_basis,
  variable_covariance,
  occasion_covariance,
  tolerance = 1e-6,
  max_iter = 1000L,
  display = FALSE
) {
  validate_scr_s3_inputs(
    X = X,
    membership = membership,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = tolerance,
    max_iter = max_iter,
    display = display
  )

  n <- nrow(X)
  groups <- ncol(membership)
  n_variables <- nrow(variable_basis)
  variable_rank <- ncol(variable_basis)
  n_occasions <- nrow(occasion_basis)
  occasion_rank <- ncol(occasion_basis)

  U <- membership
  TB <- variable_basis
  TC <- occasion_basis
  SV <- variable_covariance
  SO <- occasion_covariance

  group_mass <- colSums(U)
  if (any(group_mass <= 0)) {
    stop("each component must have positive initial membership mass.")
  }

  probabilities <- group_mass / n
  Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")
  M <- Xbar

  previous_likelihood <- -Inf
  difference <- Inf
  iteration <- 0L
  likelihood <- -Inf
  likelihood_trace <- numeric(0)
  converged <- FALSE

  while (iteration < max_iter && difference > tolerance) {
    iteration <- iteration + 1L

    # Occasion covariance update.
    inverse_SV <- solve(SV)
    SO_accumulator <- matrix(
      0,
      nrow = n_occasions,
      ncol = n_occasions
    )

    for (i in seq_len(n)) {
      Xi <- matrix(
        X[i, ],
        nrow = n_variables,
        ncol = n_occasions
      )

      for (g in seq_len(groups)) {
        Mg <- matrix(
          M[g, ],
          nrow = n_variables,
          ncol = n_occasions
        )
        residual <- Xi - Mg
        SO_accumulator <- SO_accumulator +
          U[i, g] * crossprod(residual, inverse_SV %*% residual)
      }
    }

    SO <- SO_accumulator / (n * n_variables)
    SO <- (SO + t(SO)) / 2
    validate_positive_definite_matrix(SO, "occasion covariance")

    # Variable covariance update.
    inverse_SO <- solve(SO)
    SV_accumulator <- matrix(
      0,
      nrow = n_variables,
      ncol = n_variables
    )

    for (i in seq_len(n)) {
      Xi <- matrix(
        X[i, ],
        nrow = n_variables,
        ncol = n_occasions
      )

      for (g in seq_len(groups)) {
        Mg <- matrix(
          M[g, ],
          nrow = n_variables,
          ncol = n_occasions
        )
        residual <- Xi - Mg
        SV_accumulator <- SV_accumulator +
          U[i, g] * residual %*% inverse_SO %*% t(residual)
      }
    }

    SV <- SV_accumulator / (n * n_occasions)
    SV <- (SV + t(SV)) / 2
    validate_positive_definite_matrix(SV, "variable covariance")

    # Variable-mode basis update.
    inverse_SO <- solve(SO)
    projected_TC <- inverse_SO %*% TC
    TC_metric <- projected_TC %*% t(projected_TC)
    WB <- matrix(0, nrow = n_variables, ncol = n_variables)

    for (g in seq_len(groups)) {
      Xbar_g <- matrix(
        Xbar[g, ],
        nrow = n_variables,
        ncol = n_occasions
      )
      WB <- WB +
        group_mass[g] * Xbar_g %*% TC_metric %*% t(Xbar_g)
    }

    SV_parts <- symmetric_matrix_roots(SV, "variable covariance")
    reduced_variable <- svd(
      SV_parts$inverse_root %*% WB %*% SV_parts$inverse_root
    )
    TB <- SV_parts$root %*%
      reduced_variable$u[, seq_len(variable_rank), drop = FALSE]

    # Occasion-mode basis update.
    inverse_SV <- solve(SV)
    projected_TB <- inverse_SV %*% TB
    TB_metric <- projected_TB %*% t(projected_TB)
    WC <- matrix(0, nrow = n_occasions, ncol = n_occasions)

    for (g in seq_len(groups)) {
      Xbar_g <- matrix(
        Xbar[g, ],
        nrow = n_variables,
        ncol = n_occasions
      )
      WC <- WC +
        group_mass[g] * t(Xbar_g) %*% TB_metric %*% Xbar_g
    }

    SO_parts <- symmetric_matrix_roots(SO, "occasion covariance")
    reduced_occasion <- svd(
      SO_parts$inverse_root %*% WC %*% SO_parts$inverse_root
    )
    TC <- SO_parts$root %*%
      reduced_occasion$u[, seq_len(occasion_rank), drop = FALSE]

    # Latent group scores and reconstructed Tucker2 means.
    inverse_SO <- solve(SO)
    inverse_SV <- solve(SV)
    reduction <- kronecker(
      inverse_SO %*% TC,
      inverse_SV %*% TB
    )
    Y <- Xbar %*% reduction
    M <- Y %*% kronecker(t(TC), t(TB))

    # Posterior update. In MATLAB, the complete transformed residual is squared
    # elementwise: ((X - M_g) * kron(iSOR, iSVR)).^2.
    SV_parts <- symmetric_matrix_roots(SV, "variable covariance")
    SO_parts <- symmetric_matrix_roots(SO, "occasion covariance")
    whitening <- kronecker(
      SO_parts$inverse_root,
      SV_parts$inverse_root
    )

    log_kernel <- matrix(0, nrow = n, ncol = groups)
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
      stop("S3 posterior normalization failed during fitting.")
    }

    U <- weighted_kernel / normalizer
    group_mass <- colSums(U)

    if (any(group_mass <= 0)) {
      stop("an S3 component acquired zero membership mass during fitting.")
    }

    Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")
    probabilities <- group_mass / n

    entropy_terms <- U * log(U)
    entropy_terms[is.nan(entropy_terms)] <- 0
    entropy <- sum(entropy_terms)

    logdet_SO <- as.numeric(
      determinant(SO, logarithm = TRUE)$modulus
    )
    logdet_SV <- as.numeric(
      determinant(SV, logarithm = TRUE)$modulus
    )

    likelihood <- -0.5 * n * n_variables * logdet_SO -
      0.5 * n * n_occasions * logdet_SV +
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
    n_variables * n_occasions +
    (groups - 1L) * variable_rank * occasion_rank +
    (n_variables - variable_rank) * variable_rank +
    (n_occasions - occasion_rank) * occasion_rank +
    (n_variables * n_variables + n_variables) / 2 +
    (n_occasions * n_occasions + n_occasions) / 2 -
    1L

  bic <- 2 * likelihood - log(n) * n_parameters

  if (difference < -1e-6) {
    warning(
      "S3 reference likelihood decreased before termination.",
      call. = FALSE
    )
  }

  if (display) {
    message(
      sprintf(
        "S3: dif=%g, iter=%d, like=%g, np=%g, BIC=%g",
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
    TC = TC,
    SO = SO,
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


#' Legacy SCR Name for the Three-Way Model
#'
#' Compatibility wrapper for the original MATLAB/R `t3mixs()` interface.
#'
#' @param X See `fit_scr_s3()`.
#' @param U Initial posterior-membership matrix.
#' @param TB Initial variable-mode basis.
#' @param TC Initial occasion-mode basis.
#' @param SV Initial variable covariance.
#' @param SO Initial occasion covariance.
#' @param eps Convergence threshold.
#' @param dis Numeric/logical display flag retained for MATLAB compatibility.
#' @param max_iter Maximum number of iterations.
#' @return A list containing the historical eight fields `U`, `TB`, `TC`,
#'   `SO`, `SV`, `Y`, `like`, and `bic`.
#' @export
t3mixs <- function(
  X,
  U,
  TB,
  TC,
  SV,
  SO,
  eps = 1e-6,
  dis = 0,
  max_iter = 1000L
) {
  fit <- fit_scr_s3(
    X = X,
    membership = U,
    variable_basis = TB,
    occasion_basis = TC,
    variable_covariance = SV,
    occasion_covariance = SO,
    tolerance = eps,
    max_iter = max_iter,
    display = isTRUE(as.logical(dis))
  )

  fit[c("U", "TB", "TC", "SO", "SV", "Y", "like", "bic")]
}


symmetric_matrix_roots <- function(matrix, label) {
  eig <- eigen(matrix, symmetric = TRUE)
  values <- eig$values

  if (
    any(!is.finite(values)) ||
      any(values <= sqrt(.Machine$double.eps))
  ) {
    stop(paste(label, "must be positive definite."))
  }

  root <- eig$vectors %*%
    diag(sqrt(values), nrow = length(values)) %*%
    t(eig$vectors)
  inverse_root <- eig$vectors %*%
    diag(1 / sqrt(values), nrow = length(values)) %*%
    t(eig$vectors)

  list(root = root, inverse_root = inverse_root)
}


validate_positive_definite_matrix <- function(matrix, label) {
  if (
    !is.matrix(matrix) ||
      !is.numeric(matrix) ||
      nrow(matrix) != ncol(matrix) ||
      anyNA(matrix) ||
      any(!is.finite(matrix))
  ) {
    stop(paste(label, "must be a finite numeric square matrix."))
  }

  if (
    max(abs(matrix - t(matrix))) >
      sqrt(.Machine$double.eps)
  ) {
    stop(paste(label, "must be symmetric."))
  }

  tryCatch(
    chol(matrix),
    error = function(e) {
      stop(
        paste(label, "must be positive definite."),
        call. = FALSE
      )
    }
  )

  invisible(TRUE)
}


validate_scr_s3_inputs <- function(
  X,
  membership,
  variable_basis,
  occasion_basis,
  variable_covariance,
  occasion_covariance,
  tolerance,
  max_iter,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 2L ||
      ncol(X) < 4L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop(
      "X must be a finite numeric matrix with at least two observations."
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
    !is.matrix(variable_basis) ||
      !is.numeric(variable_basis) ||
      nrow(variable_basis) < 2L ||
      ncol(variable_basis) < 1L ||
      ncol(variable_basis) > nrow(variable_basis) ||
      anyNA(variable_basis) ||
      any(!is.finite(variable_basis))
  ) {
    stop("variable_basis must be a finite full-rank numeric matrix.")
  }

  if (qr(variable_basis)$rank < ncol(variable_basis)) {
    stop("variable_basis must have full column rank.")
  }

  if (
    !is.matrix(occasion_basis) ||
      !is.numeric(occasion_basis) ||
      nrow(occasion_basis) < 2L ||
      ncol(occasion_basis) < 1L ||
      ncol(occasion_basis) > nrow(occasion_basis) ||
      anyNA(occasion_basis) ||
      any(!is.finite(occasion_basis))
  ) {
    stop("occasion_basis must be a finite full-rank numeric matrix.")
  }

  if (qr(occasion_basis)$rank < ncol(occasion_basis)) {
    stop("occasion_basis must have full column rank.")
  }

  n_variables <- nrow(variable_basis)
  n_occasions <- nrow(occasion_basis)

  if (ncol(X) != n_variables * n_occasions) {
    stop(
      "ncol(X) must equal nrow(variable_basis) * nrow(occasion_basis)."
    )
  }

  if (
    !identical(
      dim(variable_covariance),
      c(n_variables, n_variables)
    )
  ) {
    stop(
      "variable_covariance dimensions must match the variable mode."
    )
  }

  if (
    !identical(
      dim(occasion_covariance),
      c(n_occasions, n_occasions)
    )
  ) {
    stop(
      "occasion_covariance dimensions must match the occasion mode."
    )
  }

  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

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
