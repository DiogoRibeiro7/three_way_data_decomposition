#' Select Tucker3 SCR over the Number of Groups
#'
#' Extend Tucker3 structural selection from fixed G to joint selection over
#' `(G, P, Q, R)`. For each candidate number of groups, deterministic
#' multi-start memberships are generated (or supplied by a user initializer),
#' the existing within-G Tucker3 selector is run, and the best start is retained
#' under the active criterion.
#'
#' Both BIC and ICL use the package's higher-is-better convention. Ties are
#' resolved by preferring fewer free parameters and then the smaller number of
#' groups.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param groups Candidate numbers of mixture components.
#' @param variable_covariance Initial positive-definite J-by-J covariance.
#' @param occasion_covariance Initial positive-definite K-by-K covariance.
#' @param n_starts Positive number of membership starts per candidate G.
#' @param centroid_ranks Optional common candidate centroid ranks P. When NULL,
#'   every admissible rank from 1 through G - 1 is used for each G.
#' @param variable_ranks Candidate variable-mode ranks Q.
#' @param occasion_ranks Candidate occasion-mode ranks R.
#' @param criterion Model-selection criterion, either `"BIC"` or `"ICL"`.
#' @param criterion_tolerance Non-negative absolute tie tolerance.
#' @param seed Integer seed anchor used for deterministic default starts.
#' @param membership_initializer Optional function with arguments
#'   `X`, `groups`, `start`, and `seed`, returning an n-by-G membership
#'   matrix. NULL uses the package random-soft initializer.
#' @param tolerance Outer convergence tolerance.
#' @param max_iter Maximum outer iterations.
#' @param inner_max_iter Maximum Tucker3 HOOI iterations.
#' @param inner_tolerance Tucker3 inner-loop tolerance.
#' @param display Logical; emit compact progress messages.
#' @return An object of class `scr_tucker3_group_selection` containing the
#'   group-level comparison table, all retained within-G selections, the
#'   selected fitted model, and selected `G,P,Q,R`.
#' @export
select_scr_tucker3_groups <- function(
  X,
  groups,
  variable_covariance,
  occasion_covariance,
  n_starts = 1L,
  centroid_ranks = NULL,
  variable_ranks = seq_len(nrow(variable_covariance)),
  occasion_ranks = seq_len(nrow(occasion_covariance)),
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8,
  seed = 1L,
  membership_initializer = NULL,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  display = FALSE
) {
  criterion <- match.arg(criterion)

  validate_tucker3_group_selection_inputs(
    X = X,
    groups = groups,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    n_starts = n_starts,
    centroid_ranks = centroid_ranks,
    variable_ranks = variable_ranks,
    occasion_ranks = occasion_ranks,
    criterion_tolerance = criterion_tolerance,
    seed = seed,
    membership_initializer = membership_initializer,
    display = display
  )

  groups <- sort(unique(as.integer(groups)))
  variables <- nrow(variable_covariance)
  occasions <- nrow(occasion_covariance)

  retained <- vector("list", length(groups))
  names(retained) <- as.character(groups)
  rows <- vector("list", length(groups))

  for (g_index in seq_along(groups)) {
    G <- groups[g_index]

    P_candidates <- if (is.null(centroid_ranks)) {
      seq_len(G - 1L)
    } else {
      sort(unique(as.integer(centroid_ranks[centroid_ranks <= G - 1L])))
    }

    if (length(P_candidates) == 0L) {
      stop(
        sprintf(
          "no admissible centroid_ranks remain for groups = %d.",
          G
        ),
        call. = FALSE
      )
    }

    best <- NULL

    for (start in seq_len(n_starts)) {
      start_seed <- as.integer(seed + 10000L * G + start)

      if (display) {
        message(
          sprintf(
            "G=%d start=%d/%d",
            G,
            start,
            n_starts
          )
        )
      }

      set.seed(start_seed)
      membership <- initialize_scr_membership(
        X = X,
        groups = G,
        start = start,
        seed = start_seed,
        initializer = membership_initializer
      )

      selection <- tryCatch(
        select_scr_tucker3_model(
          X = X,
          membership = membership,
          variable_covariance = variable_covariance,
          occasion_covariance = occasion_covariance,
          centroid_ranks = P_candidates,
          variable_ranks = variable_ranks,
          occasion_ranks = occasion_ranks,
          tolerance = tolerance,
          max_iter = max_iter,
          inner_max_iter = inner_max_iter,
          inner_tolerance = inner_tolerance,
          criterion = criterion,
          criterion_tolerance = criterion_tolerance,
          display = FALSE
        ),
        error = function(error) {
          structure(
            list(message = conditionMessage(error)),
            class = "scr_tucker3_group_start_error"
          )
        }
      )

      if (inherits(selection, "scr_tucker3_group_start_error")) {
        next
      }

      selected_row <- selection$comparison[
        selection$comparison$selected,
        ,
        drop = FALSE
      ]

      candidate <- list(
        selection = selection,
        start = start,
        seed = start_seed,
        row = selected_row
      )

      if (
        is.null(best) ||
          is_better_group_start(
            candidate = candidate,
            incumbent = best,
            criterion = criterion,
            tolerance = criterion_tolerance
          )
      ) {
        best <- candidate
      }
    }

    if (is.null(best)) {
      stop(
        sprintf(
          "all %d starts failed for groups = %d.",
          n_starts,
          G
        ),
        call. = FALSE
      )
    }

    retained[[g_index]] <- best$selection

    row <- best$row
    rows[[g_index]] <- data.frame(
      G = G,
      P = row$P,
      Q = row$Q,
      R = row$R,
      parameters = row$parameters,
      log_likelihood = row$log_likelihood,
      bic = row$bic,
      entropy = row$entropy,
      icl = row$icl,
      converged = row$converged,
      inner_converged = row$inner_converged,
      iterations = row$iterations,
      start = as.integer(best$start),
      start_seed = as.integer(best$seed),
      stringsAsFactors = FALSE
    )
  }

  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL

  selected_index <- select_best_group_candidate(
    comparison = comparison,
    criterion = criterion,
    criterion_tolerance = criterion_tolerance
  )
  comparison$selected <- seq_len(nrow(comparison)) == selected_index

  selected_G <- comparison$G[selected_index]
  selected_selection <- retained[[as.character(selected_G)]]

  structure(
    list(
      comparison = comparison,
      selections = retained,
      selected_selection = selected_selection,
      selected_model = selected_selection$selected_model,
      selected_structure = c(
        G = selected_G,
        P = comparison$P[selected_index],
        Q = comparison$Q[selected_index],
        R = comparison$R[selected_index]
      ),
      criterion = criterion,
      criterion_value = comparison[[tolower(criterion)]][selected_index],
      criterion_tolerance = criterion_tolerance,
      seed = as.integer(seed),
      n_starts = as.integer(n_starts)
    ),
    class = "scr_tucker3_group_selection"
  )
}


is_better_group_start <- function(
  candidate,
  incumbent,
  criterion,
  tolerance
) {
  column <- tolower(criterion)
  candidate_value <- candidate$row[[column]][1L]
  incumbent_value <- incumbent$row[[column]][1L]

  if (candidate_value > incumbent_value + tolerance) {
    return(TRUE)
  }
  if (incumbent_value > candidate_value + tolerance) {
    return(FALSE)
  }

  candidate_parameters <- candidate$row$parameters[1L]
  incumbent_parameters <- incumbent$row$parameters[1L]

  if (candidate_parameters < incumbent_parameters) {
    return(TRUE)
  }
  if (candidate_parameters > incumbent_parameters) {
    return(FALSE)
  }

  candidate$start < incumbent$start
}


select_best_group_candidate <- function(
  comparison,
  criterion = c("BIC", "ICL"),
  criterion_tolerance = 1e-8
) {
  criterion <- match.arg(criterion)
  validate_bic_tolerance(criterion_tolerance)

  required <- c(
    "G",
    "P",
    "Q",
    "R",
    "parameters",
    tolower(criterion)
  )

  if (
    !is.data.frame(comparison) ||
      nrow(comparison) < 1L ||
      !all(required %in% names(comparison))
  ) {
    stop(
      "comparison must contain G, P, Q, R, parameters, and the selected criterion."
    )
  }

  values <- comparison[[tolower(criterion)]]

  if (
    anyNA(values) ||
      any(!is.finite(values)) ||
      anyNA(comparison$parameters) ||
      any(!is.finite(comparison$parameters))
  ) {
    stop("group-selection criterion values and parameter counts must be finite.")
  }

  best_value <- max(values)
  tied <- which(best_value - values <= criterion_tolerance)

  ordering <- order(
    comparison$parameters[tied],
    comparison$G[tied],
    comparison$P[tied] + comparison$Q[tied] + comparison$R[tied],
    comparison$P[tied],
    comparison$Q[tied],
    comparison$R[tied]
  )

  tied[ordering[1L]]
}


validate_tucker3_group_selection_inputs <- function(
  X,
  groups,
  variable_covariance,
  occasion_covariance,
  n_starts,
  centroid_ranks,
  variable_ranks,
  occasion_ranks,
  criterion_tolerance,
  seed,
  membership_initializer,
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

  validate_positive_integer_vector(groups, "groups")
  if (any(groups < 2L)) {
    stop("groups must contain integers greater than or equal to two.")
  }
  if (any(groups >= nrow(X))) {
    stop("groups must be smaller than the number of observations.")
  }

  validate_model_selection_dimension(n_starts, "n_starts")

  if (!is.null(centroid_ranks)) {
    if (
      !is.numeric(centroid_ranks) ||
        length(centroid_ranks) < 1L ||
        anyNA(centroid_ranks) ||
        any(!is.finite(centroid_ranks)) ||
        any(centroid_ranks < 1L) ||
        any(centroid_ranks %% 1 != 0)
    ) {
      stop("centroid_ranks must be NULL or contain positive integers.")
    }
  }

  validate_rank_candidates(
    variable_ranks,
    upper = nrow(variable_covariance),
    name = "variable_ranks"
  )
  validate_rank_candidates(
    occasion_ranks,
    upper = nrow(occasion_covariance),
    name = "occasion_ranks"
  )

  validate_bic_tolerance(criterion_tolerance)

  if (
    !is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed %% 1 != 0 ||
      seed < 0
  ) {
    stop("seed must be a non-negative integer.")
  }

  if (
    !is.null(membership_initializer) &&
      !is.function(membership_initializer)
  ) {
    stop("membership_initializer must be NULL or a function.")
  }

  if (!is.logical(display) || length(display) != 1L || is.na(display)) {
    stop("display must be TRUE or FALSE.")
  }

  invisible(TRUE)
}
