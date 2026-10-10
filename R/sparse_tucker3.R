#' Sparse Covariance-Metric Loading Projection
#'
#' Apply row-group soft thresholding to a loading matrix and re-orthonormalize
#' its columns in a covariance metric.
#'
#' For row j with Euclidean norm r_j, the thresholded row is
#'
#' `(1 - penalty / r_j)_+ * row_j`.
#'
#' The surviving matrix is then right-multiplied by the inverse square root of
#' its covariance-metric Gram matrix so that
#'
#' `t(B) %*% solve(Sigma, B) = I`.
#'
#' Right multiplication preserves rows that were thresholded exactly to zero.
#'
#' @param basis Numeric loading matrix.
#' @param covariance Positive-definite covariance defining the metric.
#' @param penalty Non-negative row-group penalty.
#' @param tolerance Numerical rank and identity tolerance.
#' @return A list containing the sparse metric-orthonormal basis, active support,
#'   original row norms, threshold multipliers, and metric Gram matrix.
#' @examples
#' basis <- diag(3)[, 1:2, drop = FALSE]
#' sparse <- scr_sparse_metric_basis(
#'   basis,
#'   covariance = diag(3),
#'   penalty = 0.1
#' )
#' sparse$active
#' #' @export
scr_sparse_metric_basis <- function(
  basis,
  covariance,
  penalty = 0,
  tolerance = 1e-10
) {
  validate_sparse_basis_inputs(
    basis = basis,
    covariance = covariance,
    penalty = penalty,
    tolerance = tolerance
  )

  row_norms <- sqrt(rowSums(basis^2))
  multipliers <- numeric(length(row_norms))

  positive <- row_norms > 0
  multipliers[positive] <- pmax(
    0,
    1 - penalty / row_norms[positive]
  )

  thresholded <- sweep(
    basis,
    MARGIN = 1L,
    STATS = multipliers,
    FUN = "*"
  )

  active <- which(rowSums(thresholded^2) > tolerance^2)

  if (length(active) < ncol(basis)) {
    stop(
      "sparse thresholding leaves fewer active rows than the requested rank.",
      call. = FALSE
    )
  }

  inverse_covariance_basis <- solve(covariance, thresholded)
  gram <- crossprod(thresholded, inverse_covariance_basis)
  gram <- (gram + t(gram)) / 2

  decomposition <- eigen(gram, symmetric = TRUE)
  if (
    any(!is.finite(decomposition$values)) ||
      min(decomposition$values) <= tolerance
  ) {
    stop(
      "sparse thresholding makes the loading basis rank deficient in the covariance metric.",
      call. = FALSE
    )
  }

  if (
    penalty == 0 &&
      max(abs(gram - diag(ncol(basis)))) <= tolerance
  ) {
    sparse_basis <- basis
  } else {
    inverse_root <- decomposition$vectors %*%
      diag(
        1 / sqrt(decomposition$values),
        nrow = length(decomposition$values)
      ) %*%
      t(decomposition$vectors)

    sparse_basis <- thresholded %*% inverse_root
  }

  metric_gram <- crossprod(
    sparse_basis,
    solve(covariance, sparse_basis)
  )

  if (
    max(
      abs(
        metric_gram - diag(ncol(sparse_basis))
      )
    ) > sqrt(tolerance)
  ) {
    stop(
      "failed to restore covariance-metric orthonormality.",
      call. = FALSE
    )
  }

  zero_rows <- setdiff(seq_len(nrow(basis)), active)
  if (
    length(zero_rows) &&
      max(abs(sparse_basis[zero_rows, , drop = FALSE])) >
        sqrt(tolerance)
  ) {
    stop("metric re-normalization changed a thresholded zero row.")
  }

  list(
    basis = sparse_basis,
    active = as.integer(active),
    inactive = as.integer(zero_rows),
    row_norms = row_norms,
    multipliers = multipliers,
    metric_gram = metric_gram,
    penalty = penalty
  )
}


#' Sparse Tucker3 Mean Projection
#'
#' Sparsify the variable and occasion loading matrices of a Tucker3 mean model,
#' keep the centroid-mode basis fixed, and re-estimate the Tucker core by
#' orthogonal projection in the weighted whitened centroid tensor.
#'
#' This is a post-fit sparse projection. It does not optimize a penalized
#' mixture likelihood.
#'
#' @param group_centroids Numeric G-by-(J*K) posterior group centroids.
#' @param group_mass Positive posterior group masses.
#' @param centroid_basis Numeric G-by-P centroid loading matrix.
#' @param variable_basis Numeric J-by-Q variable loading matrix.
#' @param occasion_basis Numeric K-by-R occasion loading matrix.
#' @param variable_covariance Positive-definite J-by-J covariance.
#' @param occasion_covariance Positive-definite K-by-K covariance.
#' @param variable_penalty Non-negative row-group penalty for variables.
#' @param occasion_penalty Non-negative row-group penalty for occasions.
#' @param tolerance Numerical tolerance.
#' @return A list containing sparse loading matrices, re-estimated Tucker core,
#'   reconstructed component means, active supports, and weighted residual.
#' @export
scr_sparse_tucker3_projection <- function(
  group_centroids,
  group_mass,
  centroid_basis,
  variable_basis,
  occasion_basis,
  variable_covariance,
  occasion_covariance,
  variable_penalty = 0,
  occasion_penalty = 0,
  tolerance = 1e-10
) {
  validate_sparse_tucker3_projection_inputs(
    group_centroids = group_centroids,
    group_mass = group_mass,
    centroid_basis = centroid_basis,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    variable_penalty = variable_penalty,
    occasion_penalty = occasion_penalty,
    tolerance = tolerance
  )

  groups <- nrow(group_centroids)
  variables <- nrow(variable_basis)
  occasions <- nrow(occasion_basis)
  probabilities <- group_mass / sum(group_mass)

  grand_vector <- as.numeric(
    crossprod(probabilities, group_centroids)
  )
  grand_mean <- matrix(
    grand_vector,
    nrow = variables,
    ncol = occasions
  )

  sparse_variable <- scr_sparse_metric_basis(
    basis = variable_basis,
    covariance = variable_covariance,
    penalty = variable_penalty,
    tolerance = tolerance
  )
  sparse_occasion <- scr_sparse_metric_basis(
    basis = occasion_basis,
    covariance = occasion_covariance,
    penalty = occasion_penalty,
    tolerance = tolerance
  )

  variable_roots <- symmetric_matrix_roots(
    variable_covariance,
    "variable covariance"
  )
  occasion_roots <- symmetric_matrix_roots(
    occasion_covariance,
    "occasion covariance"
  )

  whitened_variable <- variable_roots$inverse_root %*%
    sparse_variable$basis
  whitened_occasion <- occasion_roots$inverse_root %*%
    sparse_occasion$basis

  weighted_group_basis <- sweep(
    centroid_basis,
    MARGIN = 1L,
    STATS = sqrt(group_mass),
    FUN = "*"
  )

  group_gram <- crossprod(weighted_group_basis)
  if (
    max(
      abs(
        group_gram -
          diag(ncol(weighted_group_basis))
      )
    ) > sqrt(tolerance)
  ) {
    stop(
      "centroid_basis must be orthonormal in the group-mass metric.",
      call. = FALSE
    )
  }

  weighted_tensor <- array(
    0,
    dim = c(groups, variables, occasions)
  )

  for (g in seq_len(groups)) {
    centroid <- matrix(
      group_centroids[g, ],
      nrow = variables,
      ncol = occasions
    )
    centered <- centroid - grand_mean
    whitened <- variable_roots$inverse_root %*%
      centered %*%
      occasion_roots$inverse_root

    weighted_tensor[g, , ] <-
      sqrt(group_mass[g]) * whitened
  }

  core <- weighted_tensor
  core <- array_mode_product(
    core,
    t(weighted_group_basis),
    1L
  )
  core <- array_mode_product(
    core,
    t(whitened_variable),
    2L
  )
  core <- array_mode_product(
    core,
    t(whitened_occasion),
    3L
  )

  reconstructed_weighted <- core
  reconstructed_weighted <- array_mode_product(
    reconstructed_weighted,
    weighted_group_basis,
    1L
  )
  reconstructed_weighted <- array_mode_product(
    reconstructed_weighted,
    whitened_variable,
    2L
  )
  reconstructed_weighted <- array_mode_product(
    reconstructed_weighted,
    whitened_occasion,
    3L
  )

  weighted_residual <- sqrt(
    sum(
      (weighted_tensor - reconstructed_weighted)^2
    )
  )

  means <- scr_tucker3_means(
    grand_mean = grand_mean,
    centroid_basis = centroid_basis,
    variable_basis = sparse_variable$basis,
    occasion_basis = sparse_occasion$basis,
    core = core,
    probabilities = probabilities,
    tolerance = max(tolerance, 1e-8)
  )

  list(
    grand_mean = grand_mean,
    centroid_basis = centroid_basis,
    variable_basis = sparse_variable$basis,
    occasion_basis = sparse_occasion$basis,
    core = core,
    means = means,
    weighted_residual = weighted_residual,
    active_variables = sparse_variable$active,
    inactive_variables = sparse_variable$inactive,
    active_occasions = sparse_occasion$active,
    inactive_occasions = sparse_occasion$inactive,
    variable_penalty = variable_penalty,
    occasion_penalty = occasion_penalty,
    variable_details = sparse_variable,
    occasion_details = sparse_occasion
  )
}


#' Evaluate a Sparse Tucker3 Support Path
#'
#' Apply sparse Tucker3 projection over a deterministic grid of variable and
#' occasion penalties. Infeasible penalty pairs are retained explicitly rather
#' than terminating the path.
#'
#' @inheritParams scr_sparse_tucker3_projection
#' @param variable_penalties Numeric vector of non-negative variable penalties.
#' @param occasion_penalties Numeric vector of non-negative occasion penalties.
#' @return An object of class `scr_sparse_tucker3_path` containing a tidy
#'   support/error table and the per-grid projection objects.
#' @export
scr_sparse_tucker3_path <- function(
  group_centroids,
  group_mass,
  centroid_basis,
  variable_basis,
  occasion_basis,
  variable_covariance,
  occasion_covariance,
  variable_penalties,
  occasion_penalties,
  tolerance = 1e-10
) {
  validate_penalty_vector(
    variable_penalties,
    "variable_penalties"
  )
  validate_penalty_vector(
    occasion_penalties,
    "occasion_penalties"
  )

  variable_penalties <- sort(unique(variable_penalties))
  occasion_penalties <- sort(unique(occasion_penalties))

  grid <- expand.grid(
    variable_penalty = variable_penalties,
    occasion_penalty = occasion_penalties,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[
    order(
      grid$variable_penalty,
      grid$occasion_penalty
    ),
    ,
    drop = FALSE
  ]
  rownames(grid) <- NULL

  projections <- vector("list", nrow(grid))
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    result <- tryCatch(
      scr_sparse_tucker3_projection(
        group_centroids = group_centroids,
        group_mass = group_mass,
        centroid_basis = centroid_basis,
        variable_basis = variable_basis,
        occasion_basis = occasion_basis,
        variable_covariance = variable_covariance,
        occasion_covariance = occasion_covariance,
        variable_penalty = grid$variable_penalty[i],
        occasion_penalty = grid$occasion_penalty[i],
        tolerance = tolerance
      ),
      error = function(error) {
        structure(
          list(message = conditionMessage(error)),
          class = "scr_sparse_tucker3_path_error"
        )
      }
    )

    projections[[i]] <- result

    if (inherits(result, "scr_sparse_tucker3_path_error")) {
      rows[[i]] <- data.frame(
        variable_penalty = grid$variable_penalty[i],
        occasion_penalty = grid$occasion_penalty[i],
        feasible = FALSE,
        n_active_variables = NA_integer_,
        n_active_occasions = NA_integer_,
        active_variables = NA_character_,
        active_occasions = NA_character_,
        weighted_residual = NA_real_,
        error = result$message,
        stringsAsFactors = FALSE
      )
      next
    }

    rows[[i]] <- data.frame(
      variable_penalty = grid$variable_penalty[i],
      occasion_penalty = grid$occasion_penalty[i],
      feasible = TRUE,
      n_active_variables = length(result$active_variables),
      n_active_occasions = length(result$active_occasions),
      active_variables = collapse_support(
        result$active_variables
      ),
      active_occasions = collapse_support(
        result$active_occasions
      ),
      weighted_residual = result$weighted_residual,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL

  structure(
    list(
      summary = summary,
      projections = projections
    ),
    class = "scr_sparse_tucker3_path"
  )
}


collapse_support <- function(index) {
  if (!length(index)) {
    return("")
  }

  paste(index, collapse = ",")
}


validate_sparse_basis_inputs <- function(
  basis,
  covariance,
  penalty,
  tolerance
) {
  if (
    !is.matrix(basis) ||
      !is.numeric(basis) ||
      nrow(basis) < 1L ||
      ncol(basis) < 1L ||
      anyNA(basis) ||
      any(!is.finite(basis)) ||
      qr(basis)$rank < ncol(basis)
  ) {
    stop("basis must be a finite full-column-rank numeric matrix.")
  }

  validate_positive_definite_matrix(
    covariance,
    "covariance"
  )

  if (
    !identical(
      dim(covariance),
      c(nrow(basis), nrow(basis))
    )
  ) {
    stop("covariance dimensions must match basis rows.")
  }

  if (
    !is.numeric(penalty) ||
      length(penalty) != 1L ||
      is.na(penalty) ||
      !is.finite(penalty) ||
      penalty < 0
  ) {
    stop("penalty must be a non-negative finite number.")
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

  invisible(TRUE)
}


validate_sparse_tucker3_projection_inputs <- function(
  group_centroids,
  group_mass,
  centroid_basis,
  variable_basis,
  occasion_basis,
  variable_covariance,
  occasion_covariance,
  variable_penalty,
  occasion_penalty,
  tolerance
) {
  validate_sparse_basis_inputs(
    basis = variable_basis,
    covariance = variable_covariance,
    penalty = variable_penalty,
    tolerance = tolerance
  )
  validate_sparse_basis_inputs(
    basis = occasion_basis,
    covariance = occasion_covariance,
    penalty = occasion_penalty,
    tolerance = tolerance
  )

  if (
    !is.matrix(group_centroids) ||
      !is.numeric(group_centroids) ||
      anyNA(group_centroids) ||
      any(!is.finite(group_centroids))
  ) {
    stop("group_centroids must be a finite numeric matrix.")
  }

  groups <- nrow(group_centroids)
  expected_columns <- nrow(variable_basis) *
    nrow(occasion_basis)

  if (ncol(group_centroids) != expected_columns) {
    stop(
      "group_centroids columns must equal variables times occasions."
    )
  }

  if (
    !is.numeric(group_mass) ||
      length(group_mass) != groups ||
      anyNA(group_mass) ||
      any(!is.finite(group_mass)) ||
      any(group_mass <= 0)
  ) {
    stop(
      "group_mass must contain one positive finite value per group."
    )
  }

  if (
    !is.matrix(centroid_basis) ||
      !is.numeric(centroid_basis) ||
      nrow(centroid_basis) != groups ||
      ncol(centroid_basis) < 1L ||
      anyNA(centroid_basis) ||
      any(!is.finite(centroid_basis)) ||
      qr(centroid_basis)$rank < ncol(centroid_basis)
  ) {
    stop(
      "centroid_basis must be a finite full-rank G-by-P matrix."
    )
  }

  probabilities <- group_mass / sum(group_mass)
  if (
    max(
      abs(
        crossprod(
          probabilities,
          centroid_basis
        )
      )
    ) > sqrt(tolerance)
  ) {
    stop(
      "centroid_basis must satisfy the probability-weighted centering constraint."
    )
  }

  invisible(TRUE)
}


validate_penalty_vector <- function(values, name) {
  if (
    !is.numeric(values) ||
      length(values) < 1L ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < 0)
  ) {
    stop(
      sprintf(
        "%s must contain non-negative finite values.",
        name
      )
    )
  }

  invisible(TRUE)
}
