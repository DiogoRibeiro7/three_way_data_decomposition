#' SCR Simulation Scenario Configuration
#'
#' Return the dimensional configuration used by the Rocci-Vichi-Ranalli
#' simulation study. Scenario I differs between the public MATLAB code and the
#' published 2025 article: the MATLAB file uses four occasions, while the paper
#' reports five. This function keeps those targets explicit rather than silently
#' choosing one.
#'
#' @param scenario Either `"I"` or `"II"`.
#' @param source Either `"matlab"` for the public reference code or `"paper"`
#'   for the dimensions reported in the 2025 article.
#' @return A list containing `J`, `Q`, `K`, and `R`.
#' @export
scr_simulation_config <- function(
  scenario = c("I", "II"),
  source = c("matlab", "paper")
) {
  scenario <- match.arg(scenario)
  source <- match.arg(source)

  if (scenario == "I") {
    occasions <- if (source == "matlab") 4L else 5L

    return(list(
      scenario = scenario,
      source = source,
      J = 5L,
      Q = 2L,
      K = occasions,
      R = 2L
    ))
  }

  list(
    scenario = scenario,
    source = source,
    J = 20L,
    Q = 5L,
    K = 5L,
    R = 2L
  )
}


#' Run the SCR Baseline Simulation
#'
#' Rebuild the simulation comparison of S3, S2, and H using the validated R
#' estimators in this package.
#'
#' The public MATLAB scripts call two helper routines that are not included in
#' the authors' public repository: `preproa(X, 1)` and a custom k-means
#' initializer. Their exact behaviour is therefore not guessed here.
#' Preprocessing and membership initialization are explicit hooks. By default,
#' data are left unchanged and every random start uses a random soft partition.
#'
#' @param n Number of observations per generated data set.
#' @param groups Number of mixture components.
#' @param n_starts Number of starting points per fitted model.
#' @param dgp Integer from 1 to 4 identifying the data-generating process.
#' @param n_simulations Number of generated data sets.
#' @param scenario Either `"I"` or `"II"`.
#' @param source Either `"matlab"` or `"paper"`; see
#'   `scr_simulation_config()`.
#' @param tolerance Convergence threshold passed to S3, S2, and H.
#' @param max_iter Maximum iterations per fitted model.
#' @param preprocess Optional function applied to each generated observation
#'   matrix. `NULL` means no preprocessing.
#' @param membership_initializer Optional function with arguments
#'   `X`, `groups`, `start`, and `seed`, returning an observation-by-group
#'   membership matrix. `NULL` uses a random soft partition.
#' @param keep_data Logical; retain generated data and true memberships.
#' @return A list containing the ARI matrix, convergence/likelihood diagnostics,
#'   simulation configuration, and optionally generated data.
#' @export
run_scr_simulation <- function(
  n,
  groups,
  n_starts,
  dgp,
  n_simulations = 100L,
  scenario = c("I", "II"),
  source = c("matlab", "paper"),
  tolerance = 1e-6,
  max_iter = 1000L,
  preprocess = NULL,
  membership_initializer = NULL,
  keep_data = TRUE
) {
  validate_scr_simulation_controls(
    n = n,
    groups = groups,
    n_starts = n_starts,
    dgp = dgp,
    n_simulations = n_simulations,
    tolerance = tolerance,
    max_iter = max_iter,
    preprocess = preprocess,
    membership_initializer = membership_initializer,
    keep_data = keep_data
  )

  scenario <- match.arg(scenario)
  source <- match.arg(source)
  config <- scr_simulation_config(scenario, source)

  ari <- matrix(
    NA_real_,
    nrow = n_simulations,
    ncol = 3L,
    dimnames = list(NULL, c("S3", "S2", "H"))
  )
  diagnostics <- vector("list", n_simulations)

  if (keep_data) {
    generated_X <- array(
      NA_real_,
      dim = c(n, config$J * config$K, n_simulations)
    )
    generated_truth <- array(
      NA_real_,
      dim = c(n, groups, n_simulations)
    )
  } else {
    generated_X <- NULL
    generated_truth <- NULL
  }

  for (simulation in seq_len(n_simulations)) {
    sample <- generate_scr_scenario_sample(
      n = n,
      groups = groups,
      dgp = dgp,
      config = config,
      seed = simulation + 13L,
      preprocess = preprocess
    )

    model_results <- fit_scr_simulation_models(
      X = sample$X,
      truth = sample$U_true,
      groups = groups,
      config = config,
      n_starts = n_starts,
      simulation = simulation,
      tolerance = tolerance,
      max_iter = max_iter,
      membership_initializer = membership_initializer
    )

    ari[simulation, ] <- vapply(
      model_results,
      function(result) result$ari,
      numeric(1)
    )

    diagnostics[[simulation]] <- lapply(
      model_results,
      function(result) {
        list(
          like = result$fit$like,
          bic = result$fit$bic,
          converged = result$fit$converged,
          iterations = result$fit$it,
          start = result$start
        )
      }
    )

    if (keep_data) {
      generated_X[, , simulation] <- sample$X
      generated_truth[, , simulation] <- sample$U_true
    }
  }

  list(
    ari = ari,
    diagnostics = diagnostics,
    config = config,
    dgp = dgp,
    n = n,
    groups = groups,
    n_starts = n_starts,
    n_simulations = n_simulations,
    preprocessing = if (is.null(preprocess)) "identity" else "user-supplied",
    initialization = if (is.null(membership_initializer)) {
      "random-soft"
    } else {
      "user-supplied"
    },
    X = generated_X,
    U_true = generated_truth
  )
}


#' Legacy Scenario-I Simulation Interface
#'
#' Compatibility wrapper for `simula1()`. It now fits S3, S2, and H rather
#' than using the random placeholders present in the historical R translation.
#'
#' @param N Number of observations.
#' @param G Number of groups.
#' @param nrep Number of random starts.
#' @param dgp Data-generating process, 1 through 4.
#' @param ns Number of simulated data sets.
#' @param source Scenario-I dimensions to reproduce: `"matlab"` or `"paper"`.
#' @param ... Additional arguments passed to `run_scr_simulation()`.
#' @return A list with historical fields `ari`, `Xt`, and `Utruet`, plus
#'   the full modern result in `details`.
#' @export
simula1 <- function(
  N,
  G,
  nrep,
  dgp,
  ns = 100L,
  source = c("matlab", "paper"),
  ...
) {
  source <- match.arg(source)

  result <- run_scr_simulation(
    n = N,
    groups = G,
    n_starts = nrep,
    dgp = dgp,
    n_simulations = ns,
    scenario = "I",
    source = source,
    keep_data = TRUE,
    ...
  )

  list(
    ari = result$ari,
    Xt = result$X,
    Utruet = result$U_true,
    details = result
  )
}


#' Legacy Scenario-II Simulation Interface
#'
#' Compatibility wrapper for `simula2()`. It now fits S3, S2, and H rather
#' than using the random placeholders present in the historical R translation.
#'
#' @inheritParams simula1
#' @return A list with historical fields `ari`, `Xt`, and `Utruet`, plus
#'   the full modern result in `details`.
#' @export
simula2 <- function(
  N,
  G,
  nrep,
  dgp,
  ns = 100L,
  ...
) {
  result <- run_scr_simulation(
    n = N,
    groups = G,
    n_starts = nrep,
    dgp = dgp,
    n_simulations = ns,
    scenario = "II",
    source = "matlab",
    keep_data = TRUE,
    ...
  )

  list(
    ari = result$ari,
    Xt = result$X,
    Utruet = result$U_true,
    details = result
  )
}


generate_scr_scenario_sample <- function(
  n,
  groups,
  dgp,
  config,
  seed,
  preprocess
) {
  set.seed(seed)

  J <- config$J
  Q <- config$Q
  K <- config$K
  R <- config$R

  probabilities <- stats::runif(groups)
  probabilities <- probabilities / sum(probabilities)

  variable_factor <- matrix(stats::rnorm(J * J), nrow = J)
  variable_covariance <- crossprod(variable_factor)
  if (Q < J) {
    variable_covariance[seq_len(Q), (Q + 1L):J] <- 0
    variable_covariance[(Q + 1L):J, seq_len(Q)] <- 0
  }

  occasion_factor <- matrix(stats::rnorm(K * K), nrow = K)
  occasion_covariance <- crossprod(occasion_factor)
  if (R < K) {
    occasion_covariance[seq_len(R), (R + 1L):K] <- 0
    occasion_covariance[(R + 1L):K, seq_len(R)] <- 0
  }

  structured_covariance <- kronecker(
    occasion_covariance,
    variable_covariance
  )

  false_factor <- 0.6 * matrix(
    stats::rnorm((J * K)^2),
    nrow = J * K
  )
  false_covariance <- crossprod(false_factor)
  false_covariance <- (structured_covariance != 0) * false_covariance
  false_covariance <- (false_covariance + t(false_covariance)) / 2

  B <- diag(J)[, seq_len(Q), drop = FALSE]
  C <- diag(K)[, seq_len(R), drop = FALSE]
  eta <- matrix(
    stats::rnorm(Q * R * groups),
    nrow = Q * R,
    ncol = groups
  )
  structured_means <- t(
    20 * kronecker(C, B) %*% eta
  )
  false_means <- 20 *
    (structured_means != 0) *
    matrix(
      stats::rnorm(groups * J * K),
      nrow = groups,
      ncol = J * K
    )

  means <- if (dgp %in% c(2L, 4L)) {
    false_means
  } else {
    structured_means
  }

  covariance <- if (dgp %in% c(3L, 4L)) {
    false_covariance
  } else {
    structured_covariance
  }

  validate_positive_definite_matrix(
    covariance,
    "simulation covariance"
  )

  generated <- generate_gaussian_mixture(
    n = n,
    probabilities = probabilities,
    means = means,
    covariances = replicate(groups, covariance, simplify = FALSE)
  )

  X <- generated$X
  if (!is.null(preprocess)) {
    X <- preprocess(X)

    if (
      !is.matrix(X) ||
        !is.numeric(X) ||
        !identical(dim(X), dim(generated$X)) ||
        anyNA(X) ||
        any(!is.finite(X))
    ) {
      stop(
        "preprocess must return a finite numeric matrix with unchanged dimensions."
      )
    }
  }

  list(
    X = X,
    U_true = generated$U,
    z = generated$z,
    probabilities = probabilities,
    means = means,
    covariance = covariance
  )
}


fit_scr_simulation_models <- function(
  X,
  truth,
  groups,
  config,
  n_starts,
  simulation,
  tolerance,
  max_iter,
  membership_initializer
) {
  models <- c("S3", "S2", "H")
  results <- vector("list", length(models))
  names(results) <- models

  for (model in models) {
    best <- NULL

    for (start in seq_len(n_starts)) {
      seed <- 10L * simulation + start
      set.seed(seed)

      membership <- initialize_scr_membership(
        X = X,
        groups = groups,
        start = start,
        seed = seed,
        initializer = membership_initializer
      )

      fit <- tryCatch(
        fit_one_scr_simulation_model(
          model = model,
          X = X,
          membership = membership,
          config = config,
          tolerance = tolerance,
          max_iter = max_iter
        ),
        error = function(error) {
          structure(
            list(message = conditionMessage(error)),
            class = "scr_simulation_fit_error"
          )
        }
      )

      if (inherits(fit, "scr_simulation_fit_error")) {
        next
      }

      if (is.null(best) || fit$like > best$fit$like) {
        best <- list(fit = fit, start = start)
      }
    }

    if (is.null(best)) {
      stop(
        sprintf(
          "all %s starts failed for simulation %d.",
          model,
          simulation
        )
      )
    }

    contingency <- crossprod(
      hard_partition(truth),
      hard_partition(best$fit$U)
    )

    best$ari <- adjusted_rand_index(contingency)
    results[[model]] <- best
  }

  results
}


fit_one_scr_simulation_model <- function(
  model,
  X,
  membership,
  config,
  tolerance,
  max_iter
) {
  J <- config$J
  Q <- config$Q
  K <- config$K
  R <- config$R

  if (model == "H") {
    return(
      fit_homoscedastic_gaussian_mixture(
        X = X,
        membership = membership,
        tolerance = tolerance,
        max_iter = max_iter
      )
    )
  }

  if (model == "S2") {
    basis <- qr.Q(
      qr(
        matrix(
          stats::runif(J * K * Q * R),
          nrow = J * K,
          ncol = Q * R
        )
      )
    )[, seq_len(Q * R), drop = FALSE]

    return(
      fit_scr_s2(
        X = X,
        membership = membership,
        basis = basis,
        tolerance = tolerance,
        max_iter = max_iter
      )
    )
  }

  variable_factor <- matrix(stats::rnorm(J * J), nrow = J)
  SV <- crossprod(variable_factor)
  occasion_factor <- matrix(stats::rnorm(K * K), nrow = K)
  SO <- crossprod(occasion_factor)

  TB <- covariance_normalize_basis(
    matrix(stats::runif(J * Q), nrow = J, ncol = Q),
    SV
  )
  TC <- covariance_normalize_basis(
    matrix(stats::runif(K * R), nrow = K, ncol = R),
    SO
  )

  fit_scr_s3(
    X = X,
    membership = membership,
    variable_basis = TB,
    occasion_basis = TC,
    variable_covariance = SV,
    occasion_covariance = SO,
    tolerance = tolerance,
    max_iter = max_iter
  )
}


covariance_normalize_basis <- function(basis, covariance) {
  gram <- crossprod(basis, solve(covariance, basis))
  eig <- eigen((gram + t(gram)) / 2, symmetric = TRUE)

  if (
    any(!is.finite(eig$values)) ||
      any(eig$values <= sqrt(.Machine$double.eps))
  ) {
    stop("initial basis is rank deficient in the covariance metric.")
  }

  basis %*%
    eig$vectors %*%
    diag(1 / sqrt(eig$values), nrow = length(eig$values))
}


initialize_scr_membership <- function(
  X,
  groups,
  start,
  seed,
  initializer
) {
  if (is.null(initializer)) {
    membership <- matrix(
      stats::runif(nrow(X) * groups),
      nrow = nrow(X),
      ncol = groups
    )
    return(membership / rowSums(membership))
  }

  membership <- initializer(
    X = X,
    groups = groups,
    start = start,
    seed = seed
  )

  if (
    !is.matrix(membership) ||
      !is.numeric(membership) ||
      !identical(dim(membership), c(nrow(X), groups)) ||
      anyNA(membership) ||
      any(!is.finite(membership)) ||
      any(membership < 0)
  ) {
    stop(
      "membership_initializer must return a finite non-negative n-by-G matrix."
    )
  }

  totals <- rowSums(membership)
  if (any(totals <= 0)) {
    stop("membership_initializer returned a row with zero mass.")
  }

  membership / totals
}


validate_scr_simulation_controls <- function(
  n,
  groups,
  n_starts,
  dgp,
  n_simulations,
  tolerance,
  max_iter,
  preprocess,
  membership_initializer,
  keep_data
) {
  positive_integers <- list(
    n = n,
    groups = groups,
    n_starts = n_starts,
    n_simulations = n_simulations,
    max_iter = max_iter
  )

  for (name in names(positive_integers)) {
    value <- positive_integers[[name]]
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        is.na(value) ||
        !is.finite(value) ||
        value <= 0 ||
        value %% 1 != 0
    ) {
      stop(paste(name, "must be a positive integer."))
    }
  }

  if (groups < 2L) {
    stop("groups must be at least two.")
  }

  if (
    !is.numeric(dgp) ||
      length(dgp) != 1L ||
      is.na(dgp) ||
      !(dgp %in% 1:4)
  ) {
    stop("dgp must be one of 1, 2, 3, or 4.")
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

  if (!is.null(preprocess) && !is.function(preprocess)) {
    stop("preprocess must be NULL or a function.")
  }

  if (
    !is.null(membership_initializer) &&
      !is.function(membership_initializer)
  ) {
    stop("membership_initializer must be NULL or a function.")
  }

  if (
    !is.logical(keep_data) ||
      length(keep_data) != 1L ||
      is.na(keep_data)
  ) {
    stop("keep_data must be TRUE or FALSE.")
  }

  invisible(TRUE)
}
