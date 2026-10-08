#' Evaluate Tucker3 SCR Log-Likelihood Under a Nugget Covariance
#'
#' Evaluate the observed-data mixture log-likelihood for fixed Tucker3
#' component means, probabilities, and Kronecker covariance factors under
#'
#' `Sigma(tau) = Sigma_O %x% Sigma_V + tau I`.
#'
#' The Gaussian constant common to all compared models is omitted, matching the
#' likelihood convention used by the SCR fitters in this package.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param means Numeric matrix with one component mean per row.
#' @param probabilities Strictly positive mixing probabilities summing to one.
#' @param variable_covariance Positive-definite J-by-J covariance factor.
#' @param occasion_covariance Positive-definite K-by-K covariance factor.
#' @param nugget Non-negative scalar tau.
#' @param return_membership Logical; include posterior memberships.
#' @return A list containing log-likelihood, nugget, log determinant, covariance
#'   eigenvalues, and optionally posterior memberships.
#' @export
scr_tucker3_nugget_loglik <- function(
  X,
  means,
  probabilities,
  variable_covariance,
  occasion_covariance,
  nugget,
  return_membership = TRUE
) {
  validate_nugget_loglik_inputs(
    X = X,
    means = means,
    probabilities = probabilities,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    nugget = nugget,
    return_membership = return_membership
  )

  spectral <- nugget_spectral_decomposition(
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    nugget = nugget
  )

  n <- nrow(X)
  groups <- nrow(means)
  log_kernel <- matrix(0, nrow = n, ncol = groups)

  for (g in seq_len(groups)) {
    centered <- sweep(X, 2L, means[g, ], "-")
    transformed <- centered %*% spectral$vectors
    quadratic <- rowSums(
      sweep(
        transformed^2,
        MARGIN = 2L,
        STATS = spectral$eigenvalues,
        FUN = "/"
      )
    )

    log_kernel[, g] <- -0.5 * spectral$log_determinant -
      0.5 * quadratic +
      log(probabilities[g])
  }

  row_loglik <- row_log_sum_exp(log_kernel)
  log_likelihood <- sum(row_loglik)

  result <- list(
    log_likelihood = as.numeric(log_likelihood),
    nugget = nugget,
    log_determinant = spectral$log_determinant,
    eigenvalues = spectral$eigenvalues
  )

  if (return_membership) {
    centered_log <- sweep(
      log_kernel,
      MARGIN = 1L,
      STATS = row_loglik,
      FUN = "-"
    )
    membership <- exp(centered_log)
    membership <- membership / rowSums(membership)
    result$membership <- membership
  }

  result
}


#' Profile the Tucker3 SCR Nugget Parameter
#'
#' Optimize the non-negative nugget parameter for a fitted Tucker3 SCR model,
#' holding component means, mixing probabilities, and Kronecker covariance
#' factors fixed.
#'
#' @param X Numeric observation matrix used to fit the model.
#' @param fit Fitted Tucker3 SCR object returned by `fit_scr_s3_tucker3()`.
#' @param upper Optional finite positive upper search bound. When NULL, a
#'   scale-aware bound is constructed and expanded deterministically when the
#'   optimum lies near the current boundary.
#' @param tolerance Optimization tolerance.
#' @param boundary_tolerance Non-negative tolerance used to snap the optimum
#'   back to the separable boundary tau = 0.
#' @param expansion_factor Multiplicative factor for automatic upper-bound
#'   expansion.
#' @param max_expansions Maximum number of automatic bound expansions.
#' @return An object of class `scr_tucker3_nugget_profile` containing the
#'   profiled nugget, likelihood/BIC/ICL comparisons, posterior memberships,
#'   selected covariance diagnostics, and search metadata.
#' @export
profile_scr_tucker3_nugget <- function(
  X,
  fit,
  upper = NULL,
  tolerance = 1e-6,
  boundary_tolerance = 1e-8,
  expansion_factor = 4,
  max_expansions = 6L
) {
  validate_nugget_profile_inputs(
    X = X,
    fit = fit,
    upper = upper,
    tolerance = tolerance,
    boundary_tolerance = boundary_tolerance,
    expansion_factor = expansion_factor,
    max_expansions = max_expansions
  )

  variables <- nrow(fit$SV)
  occasions <- nrow(fit$SO)
  groups <- nrow(fit$M)

  tau0 <- scr_tucker3_nugget_loglik(
    X = X,
    means = fit$M,
    probabilities = fit$probabilities,
    variable_covariance = fit$SV,
    occasion_covariance = fit$SO,
    nugget = 0,
    return_membership = TRUE
  )

  auto_upper <- is.null(upper)
  if (auto_upper) {
    base_eigenvalues <- nugget_spectral_decomposition(
      variable_covariance = fit$SV,
      occasion_covariance = fit$SO,
      nugget = 0
    )$eigenvalues

    upper <- max(
      stats::median(base_eigenvalues),
      mean(base_eigenvalues),
      sqrt(.Machine$double.eps)
    )
  }

  objective <- function(tau) {
    -scr_tucker3_nugget_loglik(
      X = X,
      means = fit$M,
      probabilities = fit$probabilities,
      variable_covariance = fit$SV,
      occasion_covariance = fit$SO,
      nugget = tau,
      return_membership = FALSE
    )$log_likelihood
  }

  expansions <- 0L
  repeat {
    optimum <- stats::optimize(
      f = objective,
      interval = c(0, upper),
      tol = tolerance
    )

    near_upper <- optimum$minimum >= 0.95 * upper

    if (
      !auto_upper ||
        !near_upper ||
        expansions >= max_expansions
    ) {
      break
    }

    upper <- upper * expansion_factor
    expansions <- expansions + 1L
  }

  profiled_tau <- optimum$minimum
  profiled <- scr_tucker3_nugget_loglik(
    X = X,
    means = fit$M,
    probabilities = fit$probabilities,
    variable_covariance = fit$SV,
    occasion_covariance = fit$SO,
    nugget = profiled_tau,
    return_membership = TRUE
  )

  if (
    profiled_tau <= boundary_tolerance ||
      profiled$log_likelihood <=
        tau0$log_likelihood + boundary_tolerance
  ) {
    profiled_tau <- 0
    profiled <- tau0
  }

  ranks <- fit$ranks
  if (
    is.null(ranks) ||
      !all(c("P", "Q", "R") %in% names(ranks))
  ) {
    stop("fit must contain named Tucker3 ranks P, Q, and R.")
  }

  separable_parameters <- scr_tucker3_parameter_count(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = ranks[["P"]],
    variable_rank = ranks[["Q"]],
    occasion_rank = ranks[["R"]]
  )
  nugget_parameters <- separable_parameters + 1L

  n <- nrow(X)
  bic_separable <- 2 * tau0$log_likelihood -
    log(n) * separable_parameters
  bic_nugget <- 2 * profiled$log_likelihood -
    log(n) * nugget_parameters

  entropy <- scr_classification_entropy(profiled$membership)
  icl_nugget <- bic_nugget - 2 * entropy

  separable_entropy <- scr_classification_entropy(tau0$membership)
  icl_separable <- bic_separable - 2 * separable_entropy

  covariance_details <- scr_kronecker_nugget_covariance(
    variable_covariance = fit$SV,
    occasion_covariance = fit$SO,
    nugget = profiled_tau,
    details = TRUE
  )

  structure(
    list(
      nugget = profiled_tau,
      log_likelihood = profiled$log_likelihood,
      log_likelihood_separable = tau0$log_likelihood,
      likelihood_improvement =
        profiled$log_likelihood - tau0$log_likelihood,
      bic = bic_nugget,
      bic_separable = bic_separable,
      bic_improvement = bic_nugget - bic_separable,
      icl = icl_nugget,
      icl_separable = icl_separable,
      entropy = entropy,
      membership = profiled$membership,
      parameters = c(
        separable = as.integer(separable_parameters),
        nugget = as.integer(nugget_parameters)
      ),
      search = list(
        upper = upper,
        automatic_upper = auto_upper,
        expansions = as.integer(expansions),
        expansion_factor = expansion_factor,
        max_expansions = as.integer(max_expansions),
        tolerance = tolerance
      ),
      covariance = covariance_details
    ),
    class = "scr_tucker3_nugget_profile"
  )
}


nugget_spectral_decomposition <- function(
  variable_covariance,
  occasion_covariance,
  nugget
) {
  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

  if (
    !is.numeric(nugget) ||
      length(nugget) != 1L ||
      is.na(nugget) ||
      !is.finite(nugget) ||
      nugget < 0
  ) {
    stop("nugget must be a non-negative finite number.")
  }

  variable <- eigen(variable_covariance, symmetric = TRUE)
  occasion <- eigen(occasion_covariance, symmetric = TRUE)

  eigenvalues <- as.vector(
    outer(
      variable$values,
      occasion$values,
      FUN = "*"
    )
  ) + nugget

  vectors <- kronecker(
    occasion$vectors,
    variable$vectors
  )

  list(
    eigenvalues = eigenvalues,
    vectors = vectors,
    log_determinant = sum(log(eigenvalues))
  )
}


row_log_sum_exp <- function(matrix) {
  maxima <- apply(matrix, 1L, max)
  maxima + log(
    rowSums(
      exp(
        sweep(
          matrix,
          MARGIN = 1L,
          STATS = maxima,
          FUN = "-"
        )
      )
    )
  )
}


validate_nugget_loglik_inputs <- function(
  X,
  means,
  probabilities,
  variable_covariance,
  occasion_covariance,
  nugget,
  return_membership
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 1L ||
      ncol(X) < 1L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop("X must be a non-empty finite numeric matrix.")
  }

  if (
    !is.matrix(means) ||
      !is.numeric(means) ||
      nrow(means) < 2L ||
      ncol(means) != ncol(X) ||
      anyNA(means) ||
      any(!is.finite(means))
  ) {
    stop("means must be a finite G-by-ncol(X) numeric matrix.")
  }

  validate_mixture_probabilities(probabilities, nrow(means))
  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

  if (
    ncol(X) !=
      nrow(variable_covariance) * nrow(occasion_covariance)
  ) {
    stop(
      "X columns must equal variables times occasions implied by the covariance factors."
    )
  }

  if (
    !is.numeric(nugget) ||
      length(nugget) != 1L ||
      is.na(nugget) ||
      !is.finite(nugget) ||
      nugget < 0
  ) {
    stop("nugget must be a non-negative finite number.")
  }

  if (
    !is.logical(return_membership) ||
      length(return_membership) != 1L ||
      is.na(return_membership)
  ) {
    stop("return_membership must be TRUE or FALSE.")
  }

  invisible(TRUE)
}


validate_nugget_profile_inputs <- function(
  X,
  fit,
  upper,
  tolerance,
  boundary_tolerance,
  expansion_factor,
  max_expansions
) {
  required <- c("M", "SV", "SO", "probabilities", "ranks")
  if (
    !is.list(fit) ||
      !all(required %in% names(fit))
  ) {
    stop(
      "fit must be a Tucker3 fit containing M, SV, SO, probabilities, and ranks."
    )
  }

  validate_nugget_loglik_inputs(
    X = X,
    means = fit$M,
    probabilities = fit$probabilities,
    variable_covariance = fit$SV,
    occasion_covariance = fit$SO,
    nugget = 0,
    return_membership = TRUE
  )

  if (!is.null(upper)) {
    if (
      !is.numeric(upper) ||
        length(upper) != 1L ||
        is.na(upper) ||
        !is.finite(upper) ||
        upper <= 0
    ) {
      stop("upper must be NULL or a positive finite number.")
    }
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
    !is.numeric(boundary_tolerance) ||
      length(boundary_tolerance) != 1L ||
      is.na(boundary_tolerance) ||
      !is.finite(boundary_tolerance) ||
      boundary_tolerance < 0
  ) {
    stop("boundary_tolerance must be a non-negative finite number.")
  }

  if (
    !is.numeric(expansion_factor) ||
      length(expansion_factor) != 1L ||
      is.na(expansion_factor) ||
      !is.finite(expansion_factor) ||
      expansion_factor <= 1
  ) {
    stop("expansion_factor must be greater than one.")
  }

  validate_model_selection_dimension(
    max_expansions,
    "max_expansions",
    minimum = 0L
  )

  invisible(TRUE)
}
