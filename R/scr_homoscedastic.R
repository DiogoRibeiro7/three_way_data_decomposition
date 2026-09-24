#' Fit the Homoscedastic Gaussian Mixture Baseline
#'
#' Fit model H from the SCR comparison of Rocci, Vichi, and Ranalli: a Gaussian
#' mixture with unrestricted component means and one covariance matrix shared by
#' all components.
#'
#' The likelihood and BIC score follow the authors' MATLAB reference
#' implementation. In particular, the returned `bic` is
#' `2 * log_likelihood - log(n) * n_parameters`, so larger values are preferred.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Numeric initial posterior-membership matrix with
#'   observations in rows and mixture components in columns.
#' @param tolerance Positive convergence threshold for the increase in the
#'   reference log-likelihood objective.
#' @param max_iter Positive integer maximum number of iterations.
#' @param display Logical; print a compact convergence summary when `TRUE`.
#' @return A list containing posterior memberships `U`, component means
#'   `Mmu`, shared covariance `Sig`, final likelihood increment `dif`,
#'   reference log-likelihood `like`, reference BIC score `bic`, iteration
#'   count `it`, mixing proportions `probabilities`, AIC score `aic`,
#'   convergence flag `converged`, and the likelihood trajectory
#'   `log_likelihood_trace`.
#' @export
fit_homoscedastic_gaussian_mixture <- function(
  X,
  membership,
  tolerance = 1e-6,
  max_iter = 1000L,
  display = FALSE
) {
  validate_homogeneous_mixture_inputs(
    X = X,
    membership = membership,
    tolerance = tolerance,
    max_iter = max_iter,
    display = display
  )

  n <- nrow(X)
  dimension <- ncol(X)
  groups <- ncol(membership)

  U <- membership
  group_mass <- colSums(U)

  if (any(group_mass <= 0)) {
    stop("each component must have positive initial membership mass.")
  }

  probabilities <- group_mass / n
  previous_likelihood <- -Inf
  difference <- Inf
  iteration <- 0L
  likelihood_trace <- numeric(0)
  converged <- FALSE

  means <- matrix(NA_real_, nrow = groups, ncol = dimension)
  covariance <- matrix(NA_real_, nrow = dimension, ncol = dimension)
  log_kernel <- matrix(NA_real_, nrow = n, ncol = groups)
  likelihood <- -Inf

  while (iteration < max_iter && difference > tolerance) {
    iteration <- iteration + 1L

    group_mass <- colSums(U)
    if (any(group_mass <= 0)) {
      stop("a component acquired zero membership mass during fitting.")
    }

    means <- sweep(t(U) %*% X, 1L, group_mass, "/")

    covariance <- matrix(0, nrow = dimension, ncol = dimension)
    for (g in seq_len(groups)) {
      centered <- sweep(X, 2L, means[g, ], "-")
      weighted <- centered * sqrt(U[, g])
      covariance <- covariance + crossprod(weighted)
    }
    covariance <- covariance / n
    covariance <- (covariance + t(covariance)) / 2

    decomposition <- svd(covariance)
    singular_values <- decomposition$d

    if (
      any(!is.finite(singular_values)) ||
        any(singular_values <= sqrt(.Machine$double.eps))
    ) {
      stop("shared covariance became singular during fitting.")
    }

    whitening <- decomposition$v %*%
      diag(1 / sqrt(singular_values), nrow = dimension)

    for (g in seq_len(groups)) {
      centered <- sweep(X, 2L, means[g, ], "-")
      transformed <- centered %*% whitening

      log_kernel[, g] <- -0.5 * sum(log(singular_values)) -
        0.5 * rowSums(transformed^2)
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
      stop("posterior normalization failed during fitting.")
    }

    U <- weighted_kernel / normalizer
    group_mass <- colSums(U)
    probabilities <- group_mass / n

    entropy_terms <- U * log(U)
    entropy_terms[is.nan(entropy_terms)] <- 0
    entropy <- sum(entropy_terms)

    likelihood <- n * sum(probabilities * log(probabilities)) +
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
    groups * dimension +
    (dimension * dimension + dimension) / 2

  bic <- 2 * likelihood - log(n) * n_parameters
  aic <- 2 * likelihood - 2 * n_parameters

  if (difference < -1e-6) {
    warning(
      "reference likelihood decreased before termination.",
      call. = FALSE
    )
  }

  if (display) {
    message(
      sprintf(
        "H(%d): dif=%g, iter=%d, like=%g, BIC=%g, AIC=%g",
        groups,
        difference,
        iteration,
        likelihood,
        bic,
        aic
      )
    )
  }

  list(
    U = U,
    Mmu = means,
    Sig = covariance,
    dif = difference,
    like = likelihood,
    bic = bic,
    it = iteration,
    probabilities = probabilities,
    aic = aic,
    converged = converged,
    log_likelihood_trace = likelihood_trace
  )
}


#' Legacy SCR Name for the Homoscedastic Gaussian Mixture
#'
#' Compatibility wrapper for the original MATLAB/R `mixhom()` interface.
#'
#' @param X See `fit_homoscedastic_gaussian_mixture()`.
#' @param U Initial posterior-membership matrix.
#' @param eps Convergence threshold.
#' @param dis Numeric/logical display flag retained for MATLAB compatibility.
#' @param max_iter Maximum number of iterations.
#' @return A list containing the seven fields returned by the historical
#'   implementation: `U`, `Mmu`, `Sig`, `dif`, `like`, `bic`, and
#'   `it`.
#' @export
mixhom <- function(
  X,
  U,
  eps = 1e-6,
  dis = 0,
  max_iter = 1000L
) {
  fit <- fit_homoscedastic_gaussian_mixture(
    X = X,
    membership = U,
    tolerance = eps,
    max_iter = max_iter,
    display = isTRUE(as.logical(dis))
  )

  fit[c("U", "Mmu", "Sig", "dif", "like", "bic", "it")]
}


validate_homogeneous_mixture_inputs <- function(
  X,
  membership,
  tolerance,
  max_iter,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 2L ||
      ncol(X) < 1L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop("X must be a finite numeric matrix with at least two observations.")
  }

  if (
    !is.matrix(membership) ||
      !is.numeric(membership) ||
      nrow(membership) != nrow(X) ||
      ncol(membership) < 1L ||
      anyNA(membership) ||
      any(!is.finite(membership)) ||
      any(membership < 0)
  ) {
    stop(
      "membership must be a finite non-negative numeric matrix with one row per observation."
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
