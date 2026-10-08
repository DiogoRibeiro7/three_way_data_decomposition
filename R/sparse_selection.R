#' Active-Set Parameter Count for Sparse Tucker3 SCR
#'
#' Count the effective parameters of a Tucker3 mean structure conditional on
#' selected variable and occasion supports.
#'
#' The support itself is treated as fixed. This is therefore an active-set
#' parameter count and does not include an additional combinatorial search
#' penalty.
#'
#' @param groups Number of mixture components G.
#' @param variables Total number of observed variables J.
#' @param occasions Total number of observed occasions K.
#' @param centroid_rank Centroid-mode rank P.
#' @param variable_rank Variable-mode rank Q.
#' @param occasion_rank Occasion-mode rank R.
#' @param active_variables Number of active variable rows s_V.
#' @param active_occasions Number of active occasion rows s_O.
#' @param covariance_model Common covariance structure.
#' @return Integer active-set parameter count.
#' @export
scr_sparse_tucker3_parameter_count <- function(
  groups,
  variables,
  occasions,
  centroid_rank,
  variable_rank,
  occasion_rank,
  active_variables,
  active_occasions,
  covariance_model = c("separable", "nugget", "unrestricted")
) {
  covariance_model <- match.arg(covariance_model)

  validate_tucker3_dimensions(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )
  validate_model_selection_dimension(
    active_variables,
    "active_variables"
  )
  validate_model_selection_dimension(
    active_occasions,
    "active_occasions"
  )

  if (active_variables > variables) {
    stop("active_variables cannot exceed variables.")
  }
  if (active_occasions > occasions) {
    stop("active_occasions cannot exceed occasions.")
  }
  if (active_variables < variable_rank) {
    stop("active_variables cannot be smaller than variable_rank.")
  }
  if (active_occasions < occasion_rank) {
    stop("active_occasions cannot be smaller than occasion_rank.")
  }

  covariance_parameters <- scr_covariance_parameter_count(
    variables = variables,
    occasions = occasions,
    model = covariance_model
  )

  count <- groups - 1L +
    variables * occasions +
    centroid_rank * variable_rank * occasion_rank +
    centroid_rank * (groups - 1L - centroid_rank) +
    variable_rank * (active_variables - variable_rank) +
    occasion_rank * (active_occasions - occasion_rank) +
    covariance_parameters

  as.integer(count)
}


#' Select Sparse Tucker3 Penalties by BIC or ICL
#'
#' Score post-fit sparse Tucker3 projections conditionally on the fitted
#' covariance factors and mixing probabilities. Feasible penalty pairs are
#' compared using active-set BIC or ICL.
#'
#' This function does not refit the full mixture model for each support and is
#' not a penalized maximum-likelihood estimator.
#'
#' @param X Numeric observation matrix used by the fitted model.
#' @param fit Fitted Tucker3 object containing U, A, TB, TC, SV, SO,
#'   probabilities, and named ranks.
#' @param variable_penalties Non-negative variable penalties.
#' @param occasion_penalties Non-negative occasion penalties.
#' @param criterion Selection criterion, either `"BIC"` or `"ICL"`.
#' @param criterion_tolerance Non-negative numerical tie tolerance.
#' @param tolerance Numerical tolerance used by sparse projections.
#' @return An object of class `scr_sparse_tucker3_selection` containing the
#'   scored penalty grid, selected sparse projection, and posterior memberships.
#' @export
select_scr_sparse_tucker3 <- function(
  X,
  fit,
  variable_penalties,
  occasion_penalties,
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8,
  tolerance = 1e-10
) {
  criterion <- match.arg(criterion)
  validate_bic_tolerance(criterion_tolerance)
  validate_sparse_selection_fit(X, fit)

  group_mass <- colSums(fit$U)
  group_centroids <- sweep(
    t(fit$U) %*% X,
    MARGIN = 1L,
    STATS = group_mass,
    FUN = "/"
  )

  path <- scr_sparse_tucker3_path(
    group_centroids = group_centroids,
    group_mass = group_mass,
    centroid_basis = fit$A,
    variable_basis = fit$TB,
    occasion_basis = fit$TC,
    variable_covariance = fit$SV,
    occasion_covariance = fit$SO,
    variable_penalties = variable_penalties,
    occasion_penalties = occasion_penalties,
    tolerance = tolerance
  )

  groups <- ncol(fit$U)
  variables <- nrow(fit$SV)
  occasions <- nrow(fit$SO)
  ranks <- fit$ranks
  scored <- path$summary

  scored$parameters <- NA_integer_
  scored$log_likelihood <- NA_real_
  scored$entropy <- NA_real_
  scored$bic <- NA_real_
  scored$icl <- NA_real_

  memberships <- vector("list", nrow(scored))

  for (i in seq_len(nrow(scored))) {
    if (!isTRUE(scored$feasible[i])) {
      next
    }

    projection <- path$projections[[i]]

    likelihood <- scr_tucker3_nugget_loglik(
      X = X,
      means = projection$means,
      probabilities = fit$probabilities,
      variable_covariance = fit$SV,
      occasion_covariance = fit$SO,
      nugget = 0,
      return_membership = TRUE
    )

    parameters <- scr_sparse_tucker3_parameter_count(
      groups = groups,
      variables = variables,
      occasions = occasions,
      centroid_rank = ranks[["P"]],
      variable_rank = ranks[["Q"]],
      occasion_rank = ranks[["R"]],
      active_variables = length(projection$active_variables),
      active_occasions = length(projection$active_occasions),
      covariance_model = "separable"
    )

    entropy <- scr_classification_entropy(
      likelihood$membership
    )
    bic <- 2 * likelihood$log_likelihood -
      log(nrow(X)) * parameters
    icl <- bic - 2 * entropy

    scored$parameters[i] <- parameters
    scored$log_likelihood[i] <- likelihood$log_likelihood
    scored$entropy[i] <- entropy
    scored$bic[i] <- bic
    scored$icl[i] <- icl
    memberships[[i]] <- likelihood$membership
  }

  feasible <- which(
    scored$feasible &
      is.finite(scored[[tolower(criterion)]])
  )

  if (!length(feasible)) {
    stop(
      "no feasible sparse Tucker3 penalty pair produced a finite criterion.",
      call. = FALSE
    )
  }

  selected_local <- select_best_sparse_candidate(
    comparison = scored[feasible, , drop = FALSE],
    criterion = criterion,
    criterion_tolerance = criterion_tolerance
  )
  selected_index <- feasible[selected_local]
  scored$selected <- seq_len(nrow(scored)) == selected_index

  structure(
    list(
      comparison = scored,
      path = path,
      selected_projection = path$projections[[selected_index]],
      selected_membership = memberships[[selected_index]],
      selected_penalties = c(
        variable = scored$variable_penalty[selected_index],
        occasion = scored$occasion_penalty[selected_index]
      ),
      criterion = criterion,
      criterion_value =
        scored[[tolower(criterion)]][selected_index],
      criterion_tolerance = criterion_tolerance
    ),
    class = "scr_sparse_tucker3_selection"
  )
}


select_best_sparse_candidate <- function(
  comparison,
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8
) {
  criterion <- match.arg(criterion)
  validate_bic_tolerance(criterion_tolerance)

  criterion_column <- tolower(criterion)
  required <- c(
    "parameters",
    "n_active_variables",
    "n_active_occasions",
    "variable_penalty",
    "occasion_penalty",
    criterion_column
  )

  if (
    !is.data.frame(comparison) ||
      nrow(comparison) < 1L ||
      !all(required %in% names(comparison))
  ) {
    stop(
      "comparison is missing required sparse-selection columns."
    )
  }

  values <- comparison[[criterion_column]]
  if (
    anyNA(values) ||
      any(!is.finite(values)) ||
      anyNA(comparison$parameters) ||
      any(!is.finite(comparison$parameters))
  ) {
    stop("sparse-selection criterion values and parameter counts must be finite.")
  }

  best <- max(values)
  tied <- which(best - values <= criterion_tolerance)

  support_size <- comparison$n_active_variables[tied] +
    comparison$n_active_occasions[tied]

  ordering <- order(
    comparison$parameters[tied],
    support_size,
    comparison$variable_penalty[tied],
    comparison$occasion_penalty[tied]
  )

  tied[ordering[1L]]
}


validate_sparse_selection_fit <- function(X, fit) {
  required <- c(
    "U",
    "A",
    "TB",
    "TC",
    "SV",
    "SO",
    "probabilities",
    "ranks"
  )

  if (
    !is.list(fit) ||
      !all(required %in% names(fit))
  ) {
    stop(
      "fit must contain U, A, TB, TC, SV, SO, probabilities, and ranks."
    )
  }

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
    !is.matrix(fit$U) ||
      nrow(fit$U) != nrow(X) ||
      ncol(fit$U) < 2L ||
      anyNA(fit$U) ||
      any(!is.finite(fit$U)) ||
      any(fit$U < 0)
  ) {
    stop("fit$U must be a finite non-negative membership matrix matching X.")
  }

  if (max(abs(rowSums(fit$U) - 1)) > 1e-10) {
    stop("fit$U rows must sum to one.")
  }

  if (
    ncol(X) != nrow(fit$SV) * nrow(fit$SO)
  ) {
    stop("fit covariance dimensions do not match X.")
  }

  validate_positive_definite_matrix(
    fit$SV,
    "variable covariance"
  )
  validate_positive_definite_matrix(
    fit$SO,
    "occasion covariance"
  )
  validate_mixture_probabilities(
    fit$probabilities,
    ncol(fit$U)
  )

  ranks <- fit$ranks
  if (
    is.null(ranks) ||
      !all(c("P", "Q", "R") %in% names(ranks))
  ) {
    stop("fit$ranks must contain named P, Q, and R values.")
  }

  if (
    !identical(dim(fit$A), c(ncol(fit$U), as.integer(ranks[["P"]]))) ||
      !identical(dim(fit$TB), c(nrow(fit$SV), as.integer(ranks[["Q"]]))) ||
      !identical(dim(fit$TC), c(nrow(fit$SO), as.integer(ranks[["R"]])))
  ) {
    stop("fit loading dimensions are inconsistent with fit$ranks.")
  }

  invisible(TRUE)
}
