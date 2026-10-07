#' Build a Tucker3 SCR Rank Grid
#'
#' Construct the deterministic structural candidate grid used by Tucker3 SCR
#' model selection for a fixed number of mixture components.
#'
#' The centroid-mode rank is bounded by `G - 1`, not `G`. Posterior group
#' centroids are centered at their probability-weighted grand mean, so their
#' group mode lies in a contrast space of dimension at most `G - 1`.
#'
#' @param groups Number of mixture components G.
#' @param variables Number of observed variables J.
#' @param occasions Number of occasions K.
#' @param centroid_ranks Candidate centroid-mode ranks P.
#' @param variable_ranks Candidate variable-mode ranks Q.
#' @param occasion_ranks Candidate occasion-mode ranks R.
#' @return A data frame with one row per admissible `(P, Q, R)` tuple and
#'   its Tucker3 SCR parameter count.
#' @export
scr_tucker3_rank_grid <- function(
  groups,
  variables,
  occasions,
  centroid_ranks = seq_len(groups - 1L),
  variable_ranks = seq_len(variables),
  occasion_ranks = seq_len(occasions)
) {
  validate_model_selection_dimension(groups, "groups", minimum = 2L)
  validate_model_selection_dimension(variables, "variables")
  validate_model_selection_dimension(occasions, "occasions")

  centroid_ranks <- validate_rank_candidates(
    centroid_ranks,
    upper = groups - 1L,
    name = "centroid_ranks"
  )
  variable_ranks <- validate_rank_candidates(
    variable_ranks,
    upper = variables,
    name = "variable_ranks"
  )
  occasion_ranks <- validate_rank_candidates(
    occasion_ranks,
    upper = occasions,
    name = "occasion_ranks"
  )

  grid <- expand.grid(
    P = centroid_ranks,
    Q = variable_ranks,
    R = occasion_ranks,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$P, grid$Q, grid$R), , drop = FALSE]
  rownames(grid) <- NULL

  grid$parameters <- vapply(
    seq_len(nrow(grid)),
    function(i) {
      scr_tucker3_parameter_count(
        groups = groups,
        variables = variables,
        occasions = occasions,
        centroid_rank = grid$P[i],
        variable_rank = grid$Q[i],
        occasion_rank = grid$R[i]
      )
    },
    integer(1L)
  )

  grid
}


#' Select Tucker3 SCR Structure
#'
#' Fit a fixed-G Tucker3 SCR model over a joint grid of centroid, variable,
#' and occasion ranks and select the preferred structure by BIC or ICL.
#'
#' Every candidate is fitted from the same membership matrix and covariance
#' factors. This prevents the structural comparison from being confounded by
#' different initial partitions.
#'
#' The package uses the convention
#'
#' `BIC = 2 * logLik - log(n) * k`,
#'
#' so larger values are preferred. For ICL, classification entropy is
#' `H(U) = -sum(U * log(U))` and the package uses
#' `ICL = BIC - 2 * H(U)`, again with larger values preferred. Candidates
#' whose active criterion values differ by no more than `criterion_tolerance`
#' are treated as tied; the model with fewer free parameters is selected.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Initial posterior-membership matrix shared by all
#'   candidate models.
#' @param variable_covariance Initial positive-definite J-by-J covariance.
#' @param occasion_covariance Initial positive-definite K-by-K covariance.
#' @param centroid_ranks Candidate centroid-mode ranks P.
#' @param variable_ranks Candidate variable-mode ranks Q.
#' @param occasion_ranks Candidate occasion-mode ranks R.
#' @param tolerance Positive outer convergence tolerance passed to
#'   `fit_scr_s3_tucker3()`.
#' @param max_iter Positive maximum number of outer iterations.
#' @param inner_max_iter Maximum HOOI iterations in each Tucker3 mean update.
#' @param inner_tolerance Relative HOOI convergence tolerance.
#' @param criterion Model-selection criterion, either `"BIC"` or `"ICL"`.
#' @param criterion_tolerance Non-negative absolute tolerance used to identify
#'   criterion ties before preferring the lower-dimensional model.
#' @param bic_tolerance Deprecated compatibility alias for
#'   `criterion_tolerance`; when non-NULL it overrides that value.
#' @param display Logical; print candidate progress when `TRUE`.
#' @return A list with the complete comparison table, selected fitted model,
#'   selected rank tuple, criterion name, and BIC tie tolerance.
#' @export
select_scr_tucker3_model <- function(
  X,
  membership,
  variable_covariance,
  occasion_covariance,
  centroid_ranks = seq_len(ncol(membership) - 1L),
  variable_ranks = seq_len(nrow(variable_covariance)),
  occasion_ranks = seq_len(nrow(occasion_covariance)),
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8,
  bic_tolerance = NULL,
  display = FALSE
) {
  criterion <- match.arg(criterion)
  if (!is.null(bic_tolerance)) {
    criterion_tolerance <- bic_tolerance
  }

  validate_model_selection_inputs(
    X = X,
    membership = membership,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    bic_tolerance = criterion_tolerance,
    display = display
  )

  groups <- ncol(membership)
  variables <- nrow(variable_covariance)
  occasions <- nrow(occasion_covariance)

  grid <- scr_tucker3_rank_grid(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_ranks = centroid_ranks,
    variable_ranks = variable_ranks,
    occasion_ranks = occasion_ranks
  )

  validate_scr_s3_tucker3_inputs(
    X = X,
    membership = membership,
    centroid_rank = grid$P[1L],
    variable_rank = grid$Q[1L],
    occasion_rank = grid$R[1L],
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = tolerance,
    max_iter = max_iter,
    inner_max_iter = inner_max_iter,
    inner_tolerance = inner_tolerance,
    display = display
  )

  initial_membership <- membership
  initial_variable_covariance <- variable_covariance
  initial_occasion_covariance <- occasion_covariance

  fits <- vector("list", nrow(grid))
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    P <- grid$P[i]
    Q <- grid$Q[i]
    R <- grid$R[i]

    if (display) {
      message(
        sprintf(
          "Fitting Tucker3 SCR candidate P=%d, Q=%d, R=%d (%d/%d)",
          P,
          Q,
          R,
          i,
          nrow(grid)
        )
      )
    }

    fit <- tryCatch(
      fit_scr_s3_tucker3(
        X = X,
        membership = initial_membership,
        centroid_rank = P,
        variable_rank = Q,
        occasion_rank = R,
        variable_covariance = initial_variable_covariance,
        occasion_covariance = initial_occasion_covariance,
        tolerance = tolerance,
        max_iter = max_iter,
        inner_max_iter = inner_max_iter,
        inner_tolerance = inner_tolerance,
        display = FALSE
      ),
      error = function(error) {
        stop(
          sprintf(
            "Tucker3 SCR fit failed for P=%d, Q=%d, R=%d: %s",
            P,
            Q,
            R,
            conditionMessage(error)
          ),
          call. = FALSE
        )
      }
    )

    if (
      !is.numeric(fit$like) ||
        length(fit$like) != 1L ||
        !is.finite(fit$like) ||
        !is.numeric(fit$bic) ||
        length(fit$bic) != 1L ||
        !is.finite(fit$bic)
    ) {
      stop(
        sprintf(
          "Tucker3 SCR fit returned non-finite selection statistics for P=%d, Q=%d, R=%d.",
          P,
          Q,
          R
        ),
        call. = FALSE
      )
    }

    expected_bic <- 2 * fit$like -
      log(nrow(X)) * grid$parameters[i]

    if (!isTRUE(all.equal(fit$bic, expected_bic, tolerance = 1e-8))) {
      stop(
        sprintf(
          "Tucker3 SCR fit returned an inconsistent BIC for P=%d, Q=%d, R=%d.",
          P,
          Q,
          R
        ),
        call. = FALSE
      )
    }

    entropy <- scr_classification_entropy(fit$U)
    icl <- fit$bic - 2 * entropy

    fits[[i]] <- fit
    rows[[i]] <- data.frame(
      P = P,
      Q = Q,
      R = R,
      parameters = grid$parameters[i],
      log_likelihood = fit$like,
      bic = fit$bic,
      entropy = entropy,
      icl = icl,
      converged = isTRUE(fit$converged),
      inner_converged = isTRUE(fit$inner_converged),
      iterations = as.integer(fit$it),
      difference = fit$dif,
      stringsAsFactors = FALSE
    )
  }

  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL

  selected_index <- select_best_tucker3_candidate(
    comparison,
    criterion = criterion,
    criterion_tolerance = criterion_tolerance
  )
  comparison$selected <- seq_len(nrow(comparison)) == selected_index

  selected_model <- fits[[selected_index]]

  structure(
    list(
      comparison = comparison,
      selected_model = selected_model,
      selected_ranks = c(
        P = comparison$P[selected_index],
        Q = comparison$Q[selected_index],
        R = comparison$R[selected_index]
      ),
      criterion = criterion,
      criterion_value = comparison[[tolower(criterion)]][selected_index],
      criterion_tolerance = criterion_tolerance,
      bic_tolerance = criterion_tolerance
    ),
    class = "scr_tucker3_model_selection"
  )
}


select_best_tucker3_candidate <- function(
  comparison,
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8,
  bic_tolerance = NULL
) {
  criterion <- match.arg(criterion)
  if (!is.null(bic_tolerance)) {
    criterion_tolerance <- bic_tolerance
  }

  criterion_column <- tolower(criterion)
  required_columns <- c("P", "Q", "R", "parameters", criterion_column)

  if (
    !is.data.frame(comparison) ||
      nrow(comparison) < 1L ||
      !all(required_columns %in% names(comparison))
  ) {
    stop(
      "comparison must be a non-empty data frame with P, Q, R, parameters, and bic columns."
    )
  }

  validate_bic_tolerance(criterion_tolerance)

  criterion_values <- comparison[[criterion_column]]

  if (
    anyNA(criterion_values) ||
      any(!is.finite(criterion_values)) ||
      anyNA(comparison$parameters) ||
      any(!is.finite(comparison$parameters))
  ) {
    stop("comparison must contain finite BIC values and parameter counts.")
  }

  best_value <- max(criterion_values)
  tied <- which(best_value - criterion_values <= criterion_tolerance)

  rank_sum <- comparison$P[tied] +
    comparison$Q[tied] +
    comparison$R[tied]
  core_dimension <- comparison$P[tied] *
    comparison$Q[tied] *
    comparison$R[tied]

  ordering <- order(
    comparison$parameters[tied],
    rank_sum,
    core_dimension,
    comparison$P[tied],
    comparison$Q[tied],
    comparison$R[tied]
  )

  tied[ordering[1L]]
}


validate_model_selection_inputs <- function(
  X,
  membership,
  variable_covariance,
  occasion_covariance,
  bic_tolerance,
  display
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
    !is.matrix(membership) ||
      !is.numeric(membership) ||
      nrow(membership) != nrow(X) ||
      ncol(membership) < 2L ||
      anyNA(membership) ||
      any(!is.finite(membership))
  ) {
    stop(
      "membership must be a finite numeric matrix with one row per observation and at least two columns."
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
    ncol(X) !=
      nrow(variable_covariance) * nrow(occasion_covariance)
  ) {
    stop(
      "X columns must equal variables times occasions implied by the covariance factors."
    )
  }

  validate_bic_tolerance(bic_tolerance)

  if (!is.logical(display) || length(display) != 1L || is.na(display)) {
    stop("display must be TRUE or FALSE.")
  }

  invisible(TRUE)
}


validate_model_selection_dimension <- function(
  value,
  name,
  minimum = 1L
) {
  if (
    !is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value %% 1 != 0 ||
      value < minimum
  ) {
    stop(
      sprintf("%s must be an integer greater than or equal to %d.", name, minimum)
    )
  }

  invisible(TRUE)
}


validate_rank_candidates <- function(values, upper, name) {
  if (
    !is.numeric(values) ||
      length(values) < 1L ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values %% 1 != 0) ||
      any(values < 1L)
  ) {
    stop(sprintf("%s must contain positive integers.", name))
  }

  if (any(values > upper)) {
    stop(
      sprintf(
        "%s cannot contain values greater than %d.",
        name,
        upper
      )
    )
  }

  sort(unique(as.integer(values)))
}


validate_bic_tolerance <- function(bic_tolerance) {
  if (
    !is.numeric(bic_tolerance) ||
      length(bic_tolerance) != 1L ||
      is.na(bic_tolerance) ||
      !is.finite(bic_tolerance) ||
      bic_tolerance < 0
  ) {
    stop("bic_tolerance must be a non-negative finite number.")
  }

  invisible(TRUE)
}


#' Classification Entropy for Mixture Memberships
#'
#' Compute the posterior classification entropy
#' `H(U) = -sum_{i,g} u_ig log(u_ig)`.
#'
#' Zero-probability terms are interpreted by continuity as zero.
#'
#' @param membership Numeric posterior-membership matrix.
#' @return A non-negative scalar entropy.
#' @export
scr_classification_entropy <- function(membership) {
  if (
    !is.matrix(membership) ||
      !is.numeric(membership) ||
      nrow(membership) < 1L ||
      ncol(membership) < 2L ||
      anyNA(membership) ||
      any(!is.finite(membership)) ||
      any(membership < 0)
  ) {
    stop("membership must be a finite non-negative numeric matrix.")
  }

  totals <- rowSums(membership)
  if (any(totals <= 0)) {
    stop("each membership row must have positive mass.")
  }
  if (max(abs(totals - 1)) > 1e-10) {
    stop("each membership row must sum to one.")
  }

  positive <- membership > 0
  as.numeric(-sum(membership[positive] * log(membership[positive])))
}
