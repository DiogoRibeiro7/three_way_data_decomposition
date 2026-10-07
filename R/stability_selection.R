#' Evaluate Tucker3 SCR Subsampling Stability
#'
#' Measure clustering stability for one fixed Tucker3 SCR structure using
#' repeated subsampling without replacement. Each subsample is fitted from
#' multiple membership starts and retains the start with the largest
#' log-likelihood. Pairwise adjusted Rand indices are then computed only on
#' observations shared by both resamples, which avoids requiring an
#' out-of-sample prediction interface.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param groups Number of mixture components G.
#' @param centroid_rank Centroid-mode rank P.
#' @param variable_rank Variable-mode rank Q.
#' @param occasion_rank Occasion-mode rank R.
#' @param variable_covariance Initial positive-definite J-by-J covariance.
#' @param occasion_covariance Initial positive-definite K-by-K covariance.
#' @param n_resamples Number of subsamples.
#' @param subsample_fraction Fraction of observations retained in each
#'   subsample. Must be in (0, 1].
#' @param n_starts Number of membership starts per subsample.
#' @param seed Non-negative integer seed anchor.
#' @param membership_initializer Optional function with arguments
#'   `X`, `groups`, `start`, and `seed`.
#' @param tolerance Outer convergence tolerance.
#' @param max_iter Maximum outer iterations.
#' @param inner_max_iter Maximum Tucker3 HOOI iterations.
#' @param inner_tolerance Tucker3 inner-loop tolerance.
#' @param display Logical; emit compact progress messages.
#' @return An object of class `scr_tucker3_stability` containing summary
#'   stability statistics, pairwise ARI values, subsample indices, and
#'   per-resample fit diagnostics.
#' @export
scr_tucker3_stability <- function(
  X,
  groups,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  n_resamples = 20L,
  subsample_fraction = 0.8,
  n_starts = 1L,
  seed = 1L,
  membership_initializer = NULL,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  display = FALSE
) {
  validate_tucker3_stability_inputs(
    X = X,
    groups = groups,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    n_resamples = n_resamples,
    subsample_fraction = subsample_fraction,
    n_starts = n_starts,
    seed = seed,
    membership_initializer = membership_initializer,
    display = display
  )

  n <- nrow(X)
  resample_size <- max(
    2L,
    as.integer(floor(n * subsample_fraction))
  )

  if (resample_size <= groups) {
    stop(
      "subsample size must exceed the number of groups.",
      call. = FALSE
    )
  }

  resamples <- vector("list", n_resamples)
  successful <- logical(n_resamples)

  for (resample in seq_len(n_resamples)) {
    resample_seed <- as.integer(seed + 100000L + resample)
    set.seed(resample_seed)

    index <- sort(sample.int(n, size = resample_size, replace = FALSE))
    X_sub <- X[index, , drop = FALSE]

    if (display) {
      message(
        sprintf(
          "stability resample %d/%d",
          resample,
          n_resamples
        )
      )
    }

    best <- NULL

    for (start in seq_len(n_starts)) {
      start_seed <- as.integer(seed + 100000L * resample + start)
      set.seed(start_seed)

      membership <- initialize_scr_membership(
        X = X_sub,
        groups = groups,
        start = start,
        seed = start_seed,
        initializer = membership_initializer
      )

      fit <- tryCatch(
        fit_scr_s3_tucker3(
          X = X_sub,
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
          display = FALSE
        ),
        error = function(error) {
          structure(
            list(message = conditionMessage(error)),
            class = "scr_tucker3_stability_fit_error"
          )
        }
      )

      if (inherits(fit, "scr_tucker3_stability_fit_error")) {
        next
      }

      if (is.null(best) || fit$like > best$fit$like) {
        best <- list(
          fit = fit,
          start = start,
          seed = start_seed
        )
      }
    }

    if (is.null(best)) {
      resamples[[resample]] <- list(
        index = index,
        success = FALSE,
        resample_seed = resample_seed
      )
      next
    }

    successful[resample] <- TRUE
    resamples[[resample]] <- list(
      index = index,
      success = TRUE,
      resample_seed = resample_seed,
      start = as.integer(best$start),
      start_seed = as.integer(best$seed),
      like = best$fit$like,
      converged = isTRUE(best$fit$converged),
      inner_converged = isTRUE(best$fit$inner_converged),
      partition = hard_partition(best$fit$U)
    )
  }

  success_ids <- which(successful)
  if (length(success_ids) < 2L) {
    stop(
      "fewer than two Tucker3 stability resamples fitted successfully.",
      call. = FALSE
    )
  }

  pairs <- utils::combn(success_ids, 2L)
  pairwise <- vector("list", ncol(pairs))

  for (i in seq_len(ncol(pairs))) {
    left_id <- pairs[1L, i]
    right_id <- pairs[2L, i]
    left <- resamples[[left_id]]
    right <- resamples[[right_id]]

    shared <- intersect(left$index, right$index)

    if (length(shared) < 2L) {
      pairwise[[i]] <- data.frame(
        resample_a = left_id,
        resample_b = right_id,
        n_shared = length(shared),
        ari = NA_real_
      )
      next
    }

    left_rows <- match(shared, left$index)
    right_rows <- match(shared, right$index)

    contingency <- crossprod(
      left$partition[left_rows, , drop = FALSE],
      right$partition[right_rows, , drop = FALSE]
    )

    ari <- tryCatch(
      adjusted_rand_index(contingency),
      error = function(error) NA_real_
    )

    pairwise[[i]] <- data.frame(
      resample_a = left_id,
      resample_b = right_id,
      n_shared = length(shared),
      ari = ari
    )
  }

  pairwise <- do.call(rbind, pairwise)
  rownames(pairwise) <- NULL
  finite_ari <- pairwise$ari[is.finite(pairwise$ari)]

  if (length(finite_ari) == 0L) {
    stop(
      "no finite pairwise ARI values were available for stability evaluation.",
      call. = FALSE
    )
  }

  structure(
    list(
      mean_ari = mean(finite_ari),
      median_ari = stats::median(finite_ari),
      sd_ari = if (length(finite_ari) > 1L) {
        stats::sd(finite_ari)
      } else {
        NA_real_
      },
      n_pairs = length(finite_ari),
      n_resamples = as.integer(n_resamples),
      n_success = length(success_ids),
      n_failed = as.integer(n_resamples - length(success_ids)),
      resample_size = resample_size,
      subsample_fraction = subsample_fraction,
      structure = c(
        G = as.integer(groups),
        P = as.integer(centroid_rank),
        Q = as.integer(variable_rank),
        R = as.integer(occasion_rank)
      ),
      pairwise = pairwise,
      resamples = resamples,
      seed = as.integer(seed)
    ),
    class = "scr_tucker3_stability"
  )
}


#' Select a Stable Tucker3 SCR Structure
#'
#' Compare candidate Tucker3 ranks for fixed G using subsampling stability.
#' Larger mean pairwise ARI is preferred. Numerical ties are resolved by fewer
#' free parameters, then lower rank complexity.
#'
#' @inheritParams scr_tucker3_stability
#' @param centroid_ranks Candidate centroid-mode ranks P.
#' @param variable_ranks Candidate variable-mode ranks Q.
#' @param occasion_ranks Candidate occasion-mode ranks R.
#' @param stability_tolerance Absolute tolerance used for stability ties.
#' @return An object of class `scr_tucker3_stability_selection` containing
#'   the candidate table, all stability objects, and the selected ranks.
#' @export
select_scr_tucker3_stable <- function(
  X,
  groups,
  variable_covariance,
  occasion_covariance,
  centroid_ranks = seq_len(groups - 1L),
  variable_ranks = seq_len(nrow(variable_covariance)),
  occasion_ranks = seq_len(nrow(occasion_covariance)),
  n_resamples = 20L,
  subsample_fraction = 0.8,
  n_starts = 1L,
  seed = 1L,
  membership_initializer = NULL,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  stability_tolerance = 1e-8,
  display = FALSE
) {
  validate_bic_tolerance(stability_tolerance)

  grid <- scr_tucker3_rank_grid(
    groups = groups,
    variables = nrow(variable_covariance),
    occasions = nrow(occasion_covariance),
    centroid_ranks = centroid_ranks,
    variable_ranks = variable_ranks,
    occasion_ranks = occasion_ranks
  )

  results <- vector("list", nrow(grid))
  rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    P <- grid$P[i]
    Q <- grid$Q[i]
    R <- grid$R[i]

    candidate_seed <- as.integer(seed + 1000000L * i)

    stability <- scr_tucker3_stability(
      X = X,
      groups = groups,
      centroid_rank = P,
      variable_rank = Q,
      occasion_rank = R,
      variable_covariance = variable_covariance,
      occasion_covariance = occasion_covariance,
      n_resamples = n_resamples,
      subsample_fraction = subsample_fraction,
      n_starts = n_starts,
      seed = candidate_seed,
      membership_initializer = membership_initializer,
      tolerance = tolerance,
      max_iter = max_iter,
      inner_max_iter = inner_max_iter,
      inner_tolerance = inner_tolerance,
      display = display
    )

    results[[i]] <- stability
    rows[[i]] <- data.frame(
      P = P,
      Q = Q,
      R = R,
      parameters = grid$parameters[i],
      mean_ari = stability$mean_ari,
      median_ari = stability$median_ari,
      sd_ari = stability$sd_ari,
      n_pairs = stability$n_pairs,
      n_success = stability$n_success,
      n_failed = stability$n_failed,
      stringsAsFactors = FALSE
    )
  }

  comparison <- do.call(rbind, rows)
  rownames(comparison) <- NULL

  selected_index <- select_best_stability_candidate(
    comparison,
    stability_tolerance = stability_tolerance
  )
  comparison$selected <- seq_len(nrow(comparison)) == selected_index

  structure(
    list(
      comparison = comparison,
      stability = results,
      selected_stability = results[[selected_index]],
      selected_ranks = c(
        P = comparison$P[selected_index],
        Q = comparison$Q[selected_index],
        R = comparison$R[selected_index]
      ),
      stability_tolerance = stability_tolerance,
      seed = as.integer(seed)
    ),
    class = "scr_tucker3_stability_selection"
  )
}


select_best_stability_candidate <- function(
  comparison,
  stability_tolerance = 1e-8
) {
  validate_bic_tolerance(stability_tolerance)

  required <- c("P", "Q", "R", "parameters", "mean_ari")
  if (
    !is.data.frame(comparison) ||
      nrow(comparison) < 1L ||
      !all(required %in% names(comparison))
  ) {
    stop(
      "comparison must contain P, Q, R, parameters, and mean_ari."
    )
  }

  if (
    anyNA(comparison$mean_ari) ||
      any(!is.finite(comparison$mean_ari)) ||
      anyNA(comparison$parameters) ||
      any(!is.finite(comparison$parameters))
  ) {
    stop("stability values and parameter counts must be finite.")
  }

  best <- max(comparison$mean_ari)
  tied <- which(best - comparison$mean_ari <= stability_tolerance)

  rank_sum <- comparison$P[tied] +
    comparison$Q[tied] +
    comparison$R[tied]

  ordering <- order(
    comparison$parameters[tied],
    rank_sum,
    comparison$P[tied],
    comparison$Q[tied],
    comparison$R[tied]
  )

  tied[ordering[1L]]
}


validate_tucker3_stability_inputs <- function(
  X,
  groups,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  n_resamples,
  subsample_fraction,
  n_starts,
  seed,
  membership_initializer,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 3L ||
      ncol(X) < 1L ||
      anyNA(X) ||
      any(!is.finite(X))
  ) {
    stop("X must be a finite numeric matrix with at least three observations.")
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

  validate_tucker3_dimensions(
    groups = groups,
    variables = nrow(variable_covariance),
    occasions = nrow(occasion_covariance),
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  if (groups >= nrow(X)) {
    stop("groups must be smaller than the number of observations.")
  }

  validate_model_selection_dimension(n_resamples, "n_resamples", minimum = 2L)
  validate_model_selection_dimension(n_starts, "n_starts")

  if (
    !is.numeric(subsample_fraction) ||
      length(subsample_fraction) != 1L ||
      is.na(subsample_fraction) ||
      !is.finite(subsample_fraction) ||
      subsample_fraction <= 0 ||
      subsample_fraction > 1
  ) {
    stop("subsample_fraction must be in (0, 1].")
  }

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
