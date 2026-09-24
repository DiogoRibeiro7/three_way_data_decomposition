#' Update the Tucker3 Mean Block
#'
#' Estimate the Tucker3 component-mean structure for fixed posterior
#' memberships and fixed separable covariance factors.
#'
#' The posterior group centroids are centered at their probability-weighted
#' grand mean, whitened by the two covariance factors, and weighted by the
#' square root of the posterior group masses. A rank-(P,Q,R) orthogonal Tucker
#' approximation is then fitted to this weighted tensor.
#'
#' If D is the orthonormal group-mode factor of the weighted tensor and n_g is
#' the posterior mass of group g, the returned centroid basis is
#' A_g = D_g / sqrt(n_g). Consequently,
#'
#' deqn{\pi^T A = 0}
#'
#' and
#'
#' deqn{A^T diag(n_g) A = I.}
#'
#' @param group_centroids Numeric G-by-(J*K) matrix of posterior group means.
#' @param group_mass Positive posterior group masses of length G.
#' @param variable_covariance Positive-definite J-by-J covariance matrix.
#' @param occasion_covariance Positive-definite K-by-K covariance matrix.
#' @param centroid_rank Tucker3 centroid rank P.
#' @param variable_rank Tucker3 variable rank Q.
#' @param occasion_rank Tucker3 occasion rank R.
#' @param inner_max_iter Maximum HOOI iterations for the mean block.
#' @param inner_tolerance Relative reconstruction tolerance for HOOI.
#' @return A list containing the grand mean, centroid/variable/occasion bases,
#'   Tucker core, reconstructed component means, weighted residual norm, and
#'   inner-loop convergence diagnostics.
#' @export
scr_tucker3_mean_update <- function(
  group_centroids,
  group_mass,
  variable_covariance,
  occasion_covariance,
  centroid_rank,
  variable_rank,
  occasion_rank,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8
) {
  validate_tucker3_mean_update_inputs(
    group_centroids = group_centroids,
    group_mass = group_mass,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    inner_max_iter = inner_max_iter,
    inner_tolerance = inner_tolerance
  )

  groups <- nrow(group_centroids)
  n_variables <- nrow(variable_covariance)
  n_occasions <- nrow(occasion_covariance)
  probabilities <- group_mass / sum(group_mass)

  grand_vector <- as.numeric(crossprod(probabilities, group_centroids))
  grand_mean <- matrix(
    grand_vector,
    nrow = n_variables,
    ncol = n_occasions
  )

  variable_roots <- symmetric_matrix_roots(
    variable_covariance,
    "variable covariance"
  )
  occasion_roots <- symmetric_matrix_roots(
    occasion_covariance,
    "occasion covariance"
  )

  weighted_tensor <- array(
    0,
    dim = c(groups, n_variables, n_occasions)
  )

  for (g in seq_len(groups)) {
    centroid <- matrix(
      group_centroids[g, ],
      nrow = n_variables,
      ncol = n_occasions
    )
    centered <- centroid - grand_mean
    whitened <- variable_roots$inverse_root %*%
      centered %*%
      occasion_roots$inverse_root

    weighted_tensor[g, , ] <- sqrt(group_mass[g]) * whitened
  }

  weighted_norm <- sqrt(sum(weighted_tensor^2))
  if (!is.finite(weighted_norm) ||
      weighted_norm <= sqrt(.Machine$double.eps)) {
    stop("weighted centered group centroids have zero numerical norm.")
  }

  tucker <- fit_orthogonal_tucker_array(
    tensor = weighted_tensor,
    ranks = c(
      as.integer(centroid_rank),
      as.integer(variable_rank),
      as.integer(occasion_rank)
    ),
    max_iter = as.integer(inner_max_iter),
    tolerance = inner_tolerance,
    group_contrast = sqrt(group_mass)
  )

  weighted_group_basis <- tucker$factors[[1L]]
  centroid_basis <- sweep(
    weighted_group_basis,
    MARGIN = 1L,
    STATS = sqrt(group_mass),
    FUN = "/"
  )

  variable_basis <- variable_roots$root %*%
    tucker$factors[[2L]]
  occasion_basis <- occasion_roots$root %*%
    tucker$factors[[3L]]

  means <- scr_tucker3_means(
    grand_mean = grand_mean,
    centroid_basis = centroid_basis,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    core = tucker$core,
    probabilities = probabilities,
    tolerance = max(inner_tolerance, 1e-8)
  )

  list(
    grand_mean = grand_mean,
    centroid_basis = centroid_basis,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    core = tucker$core,
    means = means,
    weighted_residual = tucker$residual,
    inner_iterations = tucker$iterations,
    inner_converged = tucker$converged
  )
}


#' Fit the Tucker3 SCR Extension
#'
#' Fit the Gaussian SCR extension with centroid-, variable-, and occasion-mode
#' reduction. The covariance model remains separable,
#' Sigma_O (x) Sigma_V, while the component means follow the Tucker3 structure
#' implemented by `scr_tucker3_means()`.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Initial posterior-membership matrix.
#' @param centroid_rank Centroid-mode rank P.
#' @param variable_rank Variable-mode rank Q.
#' @param occasion_rank Occasion-mode rank R.
#' @param variable_covariance Initial positive-definite J-by-J covariance.
#' @param occasion_covariance Initial positive-definite K-by-K covariance.
#' @param tolerance Positive outer convergence tolerance.
#' @param max_iter Positive maximum number of outer iterations.
#' @param inner_max_iter Maximum HOOI iterations in each Tucker3 mean update.
#' @param inner_tolerance Relative HOOI convergence tolerance.
#' @param display Logical; print a compact convergence summary when `TRUE`.
#' @return A fitted Tucker3 SCR list containing memberships, all three loading
#'   matrices, covariance factors, Tucker core, grand mean, component means,
#'   likelihood/BIC diagnostics, and convergence traces.
#' @export
fit_scr_s3_tucker3 <- function(
  X,
  membership,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  display = FALSE
) {
  validate_scr_s3_tucker3_inputs(
    X = X,
    membership = membership,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = tolerance,
    max_iter = max_iter,
    inner_max_iter = inner_max_iter,
    inner_tolerance = inner_tolerance,
    display = display
  )

  n <- nrow(X)
  groups <- ncol(membership)
  n_variables <- nrow(variable_covariance)
  n_occasions <- nrow(occasion_covariance)

  U <- membership
  SV <- variable_covariance
  SO <- occasion_covariance

  group_mass <- colSums(U)
  probabilities <- group_mass / n
  Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")

  # The unconstrained posterior centroids are a neutral initial mean model.
  M <- Xbar

  previous_likelihood <- -Inf
  difference <- Inf
  iteration <- 0L
  likelihood <- -Inf
  likelihood_trace <- numeric(0)
  mean_residual_trace <- numeric(0)
  converged <- FALSE

  grand_mean <- NULL
  centroid_basis <- NULL
  variable_basis <- NULL
  occasion_basis <- NULL
  core <- NULL
  inner_converged <- FALSE

  while (iteration < max_iter && difference > tolerance) {
    iteration <- iteration + 1L

    covariance_update <- update_scr_separable_covariance(
      X = X,
      membership = U,
      means = M,
      n_variables = n_variables,
      n_occasions = n_occasions
    )
    SV <- covariance_update$SV
    SO <- covariance_update$SO

    mean_update <- scr_tucker3_mean_update(
      group_centroids = Xbar,
      group_mass = group_mass,
      variable_covariance = SV,
      occasion_covariance = SO,
      centroid_rank = centroid_rank,
      variable_rank = variable_rank,
      occasion_rank = occasion_rank,
      inner_max_iter = inner_max_iter,
      inner_tolerance = inner_tolerance
    )

    grand_mean <- mean_update$grand_mean
    centroid_basis <- mean_update$centroid_basis
    variable_basis <- mean_update$variable_basis
    occasion_basis <- mean_update$occasion_basis
    core <- mean_update$core
    M <- mean_update$means
    inner_converged <- mean_update$inner_converged
    mean_residual_trace <- c(
      mean_residual_trace,
      mean_update$weighted_residual
    )

    variable_roots <- symmetric_matrix_roots(
      SV,
      "variable covariance"
    )
    occasion_roots <- symmetric_matrix_roots(
      SO,
      "occasion covariance"
    )
    whitening <- kronecker(
      occasion_roots$inverse_root,
      variable_roots$inverse_root
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
      stop("Tucker3 SCR posterior normalization failed.")
    }

    U <- weighted_kernel / normalizer
    group_mass <- colSums(U)

    if (any(group_mass <= 0)) {
      stop("a Tucker3 SCR component acquired zero membership mass.")
    }

    probabilities <- group_mass / n
    Xbar <- sweep(t(U) %*% X, 1L, group_mass, "/")

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

  n_parameters <- scr_tucker3_parameter_count(
    groups = groups,
    variables = n_variables,
    occasions = n_occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )
  bic <- 2 * likelihood - log(n) * n_parameters

  if (difference < -1e-6) {
    warning(
      "Tucker3 SCR likelihood decreased before termination.",
      call. = FALSE
    )
  }

  if (display) {
    message(
      sprintf(
        "S3-Tucker3(P=%d,Q=%d,R=%d): dif=%g, iter=%d, like=%g, BIC=%g",
        centroid_rank,
        variable_rank,
        occasion_rank,
        difference,
        iteration,
        likelihood,
        bic
      )
    )
  }

  list(
    U = U,
    A = centroid_basis,
    TB = variable_basis,
    TC = occasion_basis,
    core = core,
    grand_mean = grand_mean,
    M = M,
    SO = SO,
    SV = SV,
    dif = difference,
    like = likelihood,
    bic = bic,
    it = iteration,
    probabilities = probabilities,
    converged = converged,
    inner_converged = inner_converged,
    log_likelihood_trace = likelihood_trace,
    mean_residual_trace = mean_residual_trace,
    ranks = c(
      P = as.integer(centroid_rank),
      Q = as.integer(variable_rank),
      R = as.integer(occasion_rank)
    )
  )
}


update_scr_separable_covariance <- function(
  X,
  membership,
  means,
  n_variables,
  n_occasions
) {
  n <- nrow(X)
  groups <- ncol(membership)

  # Start from identity metrics for the first flip-flop half-step.
  SV <- diag(n_variables)
  SO_accumulator <- matrix(0, n_occasions, n_occasions)

  for (i in seq_len(n)) {
    Xi <- matrix(X[i, ], nrow = n_variables, ncol = n_occasions)

    for (g in seq_len(groups)) {
      Mg <- matrix(means[g, ], nrow = n_variables, ncol = n_occasions)
      residual <- Xi - Mg
      SO_accumulator <- SO_accumulator +
        membership[i, g] * crossprod(residual, residual)
    }
  }

  SO <- SO_accumulator / (n * n_variables)
  SO <- (SO + t(SO)) / 2
  validate_positive_definite_matrix(SO, "occasion covariance")

  inverse_SO <- solve(SO)
  SV_accumulator <- matrix(0, n_variables, n_variables)

  for (i in seq_len(n)) {
    Xi <- matrix(X[i, ], nrow = n_variables, ncol = n_occasions)

    for (g in seq_len(groups)) {
      Mg <- matrix(means[g, ], nrow = n_variables, ncol = n_occasions)
      residual <- Xi - Mg
      SV_accumulator <- SV_accumulator +
        membership[i, g] * residual %*% inverse_SO %*% t(residual)
    }
  }

  SV <- SV_accumulator / (n * n_occasions)
  SV <- (SV + t(SV)) / 2
  validate_positive_definite_matrix(SV, "variable covariance")

  # Recompute SO with the updated variable metric, matching a full flip-flop
  # covariance block rather than leaving the first identity-based half-step.
  inverse_SV <- solve(SV)
  SO_accumulator[,] <- 0

  for (i in seq_len(n)) {
    Xi <- matrix(X[i, ], nrow = n_variables, ncol = n_occasions)

    for (g in seq_len(groups)) {
      Mg <- matrix(means[g, ], nrow = n_variables, ncol = n_occasions)
      residual <- Xi - Mg
      SO_accumulator <- SO_accumulator +
        membership[i, g] * crossprod(
          residual,
          inverse_SV %*% residual
        )
    }
  }

  SO <- SO_accumulator / (n * n_variables)
  SO <- (SO + t(SO)) / 2
  validate_positive_definite_matrix(SO, "occasion covariance")

  list(SV = SV, SO = SO)
}


fit_orthogonal_tucker_array <- function(
  tensor,
  ranks,
  max_iter,
  tolerance,
  group_contrast
) {
  dimensions <- dim(tensor)

  if (length(dimensions) != 3L) {
    stop("tensor must be a three-mode array.")
  }

  factors <- vector("list", 3L)

  for (mode in 1:3) {
    unfolding <- unfold_array_mode(tensor, mode)
    decomposition <- svd(
      unfolding,
      nu = ranks[mode],
      nv = 0
    )
    factors[[mode]] <- decomposition$u[
      ,
      seq_len(ranks[mode]),
      drop = FALSE
    ]
  }

  factors[[1L]] <- enforce_group_contrast(
    factors[[1L]],
    group_contrast,
    ranks[1L],
    tolerance
  )

  previous_residual <- Inf
  converged <- FALSE
  residual <- Inf

  for (iteration in seq_len(max_iter)) {
    for (mode in 1:3) {
      projected <- tensor

      for (other_mode in setdiff(1:3, mode)) {
        projected <- array_mode_product(
          projected,
          t(factors[[other_mode]]),
          other_mode
        )
      }

      unfolding <- unfold_array_mode(projected, mode)
      decomposition <- svd(
        unfolding,
        nu = ranks[mode],
        nv = 0
      )
      candidate <- decomposition$u[
        ,
        seq_len(ranks[mode]),
        drop = FALSE
      ]

      if (mode == 1L) {
        candidate <- enforce_group_contrast(
          candidate,
          group_contrast,
          ranks[1L],
          tolerance
        )
      }

      factors[[mode]] <- candidate
    }

    core <- tensor
    for (mode in 1:3) {
      core <- array_mode_product(
        core,
        t(factors[[mode]]),
        mode
      )
    }

    reconstructed <- core
    for (mode in 1:3) {
      reconstructed <- array_mode_product(
        reconstructed,
        factors[[mode]],
        mode
      )
    }

    residual <- sqrt(sum((tensor - reconstructed)^2))
    scale <- max(sqrt(sum(tensor^2)), sqrt(.Machine$double.eps))

    if (is.finite(previous_residual) &&
        abs(previous_residual - residual) / scale <= tolerance) {
      converged <- TRUE
      break
    }

    previous_residual <- residual
  }

  list(
    core = core,
    factors = factors,
    residual = residual,
    iterations = iteration,
    converged = converged
  )
}


enforce_group_contrast <- function(
  candidate,
  group_contrast,
  rank,
  tolerance
) {
  norm_squared <- sum(group_contrast^2)
  projected <- candidate -
    group_contrast %*%
      (crossprod(group_contrast, candidate) / norm_squared)

  decomposition <- qr(projected, tol = tolerance)

  if (decomposition$rank < rank) {
    stop(
      "centroid-mode update is rank deficient in the weighted contrast space."
    )
  }

  qr.Q(decomposition, complete = FALSE)[
    ,
    seq_len(rank),
    drop = FALSE
  ]
}


unfold_array_mode <- function(tensor, mode) {
  dimensions <- dim(tensor)
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  permuted <- aperm(tensor, permutation)

  matrix(
    permuted,
    nrow = dimensions[mode],
    ncol = prod(dimensions[-mode])
  )
}


array_mode_product <- function(tensor, matrix, mode) {
  dimensions <- dim(tensor)
  permutation <- c(mode, setdiff(seq_along(dimensions), mode))
  permuted <- aperm(tensor, permutation)

  unfolded <- matrix(
    permuted,
    nrow = dimensions[mode],
    ncol = prod(dimensions[-mode])
  )

  if (ncol(matrix) != nrow(unfolded)) {
    stop("mode-product matrix has incompatible dimensions.")
  }

  product <- matrix %*% unfolded
  permuted_dimensions <- c(
    nrow(matrix),
    dimensions[permutation[-1L]]
  )
  product_array <- array(product, dim = permuted_dimensions)

  aperm(product_array, order(permutation))
}


validate_tucker3_mean_update_inputs <- function(
  group_centroids,
  group_mass,
  variable_covariance,
  occasion_covariance,
  centroid_rank,
  variable_rank,
  occasion_rank,
  inner_max_iter,
  inner_tolerance
) {
  if (
    !is.matrix(group_centroids) ||
      !is.numeric(group_centroids) ||
      anyNA(group_centroids) ||
      any(!is.finite(group_centroids))
  ) {
    stop("group_centroids must be a finite numeric matrix.")
  }

  if (
    !is.numeric(group_mass) ||
      length(group_mass) != nrow(group_centroids) ||
      anyNA(group_mass) ||
      any(!is.finite(group_mass)) ||
      any(group_mass <= 0)
  ) {
    stop("group_mass must contain one positive finite value per group.")
  }

  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

  n_variables <- nrow(variable_covariance)
  n_occasions <- nrow(occasion_covariance)

  if (ncol(group_centroids) != n_variables * n_occasions) {
    stop(
      "group_centroids columns must equal variables times occasions."
    )
  }

  validate_tucker3_dimensions(
    groups = nrow(group_centroids),
    variables = n_variables,
    occasions = n_occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  if (
    !is.numeric(inner_max_iter) ||
      length(inner_max_iter) != 1L ||
      is.na(inner_max_iter) ||
      !is.finite(inner_max_iter) ||
      inner_max_iter <= 0 ||
      inner_max_iter %% 1 != 0
  ) {
    stop("inner_max_iter must be a positive integer.")
  }

  if (
    !is.numeric(inner_tolerance) ||
      length(inner_tolerance) != 1L ||
      is.na(inner_tolerance) ||
      !is.finite(inner_tolerance) ||
      inner_tolerance <= 0
  ) {
    stop("inner_tolerance must be a positive finite number.")
  }

  invisible(TRUE)
}


validate_scr_s3_tucker3_inputs <- function(
  X,
  membership,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  tolerance,
  max_iter,
  inner_max_iter,
  inner_tolerance,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 2L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop("X must be a finite numeric matrix with at least two observations.")
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
      "membership must be a finite non-negative matrix with one row per observation."
    )
  }

  row_totals <- rowSums(membership)
  if (any(row_totals <= 0) ||
      max(abs(row_totals - 1)) > 1e-10) {
    stop("each membership row must sum to one.")
  }

  validate_positive_definite_matrix(
    variable_covariance,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    occasion_covariance,
    "occasion covariance"
  )

  n_variables <- nrow(variable_covariance)
  n_occasions <- nrow(occasion_covariance)

  if (ncol(X) != n_variables * n_occasions) {
    stop(
      "ncol(X) must equal the product of covariance-factor dimensions."
    )
  }

  validate_tucker3_dimensions(
    groups = ncol(membership),
    variables = n_variables,
    occasions = n_occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  validate_tucker3_mean_update_inputs(
    group_centroids = matrix(
      0,
      nrow = ncol(membership),
      ncol = n_variables * n_occasions
    ),
    group_mass = rep(1, ncol(membership)),
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    inner_max_iter = inner_max_iter,
    inner_tolerance = inner_tolerance
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
