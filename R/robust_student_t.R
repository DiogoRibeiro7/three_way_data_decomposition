#' Evaluate Tucker3 SCR Student-t Mixture Likelihood
#'
#' Evaluate a multivariate Student-t mixture likelihood for fixed Tucker3
#' component means, mixing probabilities, and separable common scale factors.
#'
#' The returned log-likelihood uses the same convention as the package Gaussian
#' likelihood: the Gaussian constant d/2 * log(2*pi) is omitted. With this
#' convention, the Student-t log kernel converges to the package Gaussian log
#' kernel as degrees_freedom tends to infinity.
#'
#' For squared Mahalanobis distance delta_ig, the latent robustness weight is
#'
#' `w_ig = (nu + d) / (nu + delta_ig)`.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param means Numeric G-by-d component-mean matrix.
#' @param probabilities Strictly positive mixing probabilities summing to one.
#' @param variable_scale Positive-definite J-by-J Student-t scale factor.
#' @param occasion_scale Positive-definite K-by-K Student-t scale factor.
#' @param degrees_freedom Student-t degrees of freedom nu, required to exceed 2.
#' @param return_membership Logical; include posterior memberships.
#' @param return_distances Logical; include squared Mahalanobis distances.
#' @param return_weights Logical; include component and observation weights.
#' @return A list with log-likelihood and requested posterior diagnostics.
#' @export
scr_tucker3_student_loglik <- function(
  X,
  means,
  probabilities,
  variable_scale,
  occasion_scale,
  degrees_freedom,
  return_membership = TRUE,
  return_distances = TRUE,
  return_weights = TRUE
) {
  validate_student_tucker3_inputs(
    X = X,
    means = means,
    probabilities = probabilities,
    variable_scale = variable_scale,
    occasion_scale = occasion_scale,
    degrees_freedom = degrees_freedom,
    return_membership = return_membership,
    return_distances = return_distances,
    return_weights = return_weights
  )

  n <- nrow(X)
  groups <- nrow(means)
  dimension <- ncol(X)

  variable <- eigen(variable_scale, symmetric = TRUE)
  occasion <- eigen(occasion_scale, symmetric = TRUE)

  eigenvalues <- as.vector(
    outer(
      variable$values,
      occasion$values,
      FUN = "*"
    )
  )
  eigenvectors <- kronecker(
    occasion$vectors,
    variable$vectors
  )
  log_determinant <- sum(log(eigenvalues))

  log_kernel <- matrix(0, nrow = n, ncol = groups)
  distances <- matrix(0, nrow = n, ncol = groups)

  nu <- degrees_freedom
  normalization <- lgamma((nu + dimension) / 2) -
    lgamma(nu / 2) -
    0.5 * dimension * log(nu) -
    0.5 * dimension * log(pi) +
    0.5 * dimension * log(2 * pi)

  for (g in seq_len(groups)) {
    centered <- sweep(X, 2L, means[g, ], "-")
    transformed <- centered %*% eigenvectors

    delta <- rowSums(
      sweep(
        transformed^2,
        MARGIN = 2L,
        STATS = eigenvalues,
        FUN = "/"
      )
    )
    distances[, g] <- delta

    log_kernel[, g] <- log(probabilities[g]) +
      normalization -
      0.5 * log_determinant -
      0.5 * (nu + dimension) * log1p(delta / nu)
  }

  row_loglik <- row_log_sum_exp(log_kernel)
  result <- list(
    log_likelihood = as.numeric(sum(row_loglik)),
    degrees_freedom = nu,
    dimension = dimension,
    log_determinant = log_determinant
  )

  membership <- NULL
  if (return_membership || return_weights) {
    centered_log <- sweep(
      log_kernel,
      MARGIN = 1L,
      STATS = row_loglik,
      FUN = "-"
    )
    membership <- exp(centered_log)
    membership <- membership / rowSums(membership)
  }

  if (return_membership) {
    result$membership <- membership
  }

  if (return_distances) {
    result$mahalanobis_squared <- distances
  }

  if (return_weights) {
    component_weights <- (nu + dimension) /
      (nu + distances)
    observation_weights <- rowSums(
      membership * component_weights
    )

    result$component_weights <- component_weights
    result$observation_weights <- observation_weights
  }

  result
}


#' Parameter Count for Student-t Tucker3 SCR
#'
#' Extend the Gaussian Tucker3 parameter count for a Student-t model.
#'
#' @param groups Number of mixture components G.
#' @param variables Number of variables J.
#' @param occasions Number of occasions K.
#' @param centroid_rank Centroid-mode rank P.
#' @param variable_rank Variable-mode rank Q.
#' @param occasion_rank Occasion-mode rank R.
#' @param estimate_degrees_freedom Logical; add one parameter when nu is
#'   estimated rather than fixed.
#' @return Integer parameter count.
#' @export
scr_tucker3_student_parameter_count <- function(
  groups,
  variables,
  occasions,
  centroid_rank,
  variable_rank,
  occasion_rank,
  estimate_degrees_freedom = FALSE
) {
  if (
    !is.logical(estimate_degrees_freedom) ||
      length(estimate_degrees_freedom) != 1L ||
      is.na(estimate_degrees_freedom)
  ) {
    stop("estimate_degrees_freedom must be TRUE or FALSE.")
  }

  base <- scr_tucker3_parameter_count(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  as.integer(
    base + if (estimate_degrees_freedom) 1L else 0L
  )
}


validate_student_tucker3_inputs <- function(
  X,
  means,
  probabilities,
  variable_scale,
  occasion_scale,
  degrees_freedom,
  return_membership,
  return_distances,
  return_weights
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

  validate_mixture_probabilities(
    probabilities,
    nrow(means)
  )
  validate_positive_definite_matrix(
    variable_scale,
    "variable scale"
  )
  validate_positive_definite_matrix(
    occasion_scale,
    "occasion scale"
  )

  if (
    ncol(X) != nrow(variable_scale) *
      nrow(occasion_scale)
  ) {
    stop(
      "X columns must equal variables times occasions implied by the scale factors."
    )
  }

  if (
    !is.numeric(degrees_freedom) ||
      length(degrees_freedom) != 1L ||
      is.na(degrees_freedom) ||
      !is.finite(degrees_freedom) ||
      degrees_freedom <= 2
  ) {
    stop("degrees_freedom must be a finite number greater than two.")
  }

  logical_controls <- list(
    return_membership = return_membership,
    return_distances = return_distances,
    return_weights = return_weights
  )

  for (name in names(logical_controls)) {
    value <- logical_controls[[name]]
    if (
      !is.logical(value) ||
        length(value) != 1L ||
        is.na(value)
    ) {
      stop(paste(name, "must be TRUE or FALSE."))
    }
  }

  invisible(TRUE)
}
