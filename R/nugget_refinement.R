#' Alternating Nugget Refinement for Tucker3 SCR
#'
#' Refine a fixed-rank Tucker3 SCR solution by alternating between the validated
#' separable Tucker3 fitter and conditional nugget profiling. This is a
#' conservative profile-refinement algorithm, not an exact EM algorithm for the
#' nonseparable covariance model.
#'
#' Each proposed round:
#'
#' 1. refits the separable Tucker3 block from the current memberships;
#' 2. profiles the nugget parameter conditional on that fit;
#' 3. uses the profiled nugget memberships for the next round.
#'
#' The outer objective is the profiled observed-data nugget log-likelihood. A
#' proposed round that decreases this objective beyond `objective_tolerance`
#' is rejected and the last accepted state is returned.
#'
#' @param X Numeric observation matrix.
#' @param membership Initial membership matrix.
#' @param centroid_rank Tucker3 centroid rank P.
#' @param variable_rank Tucker3 variable rank Q.
#' @param occasion_rank Tucker3 occasion rank R.
#' @param variable_covariance Initial variable covariance.
#' @param occasion_covariance Initial occasion covariance.
#' @param refinement_rounds Maximum number of accepted/proposed refinement
#'   rounds after the initial Tucker3 fit and nugget profile. Zero returns the
#'   one-shot fit plus profile.
#' @param objective_tolerance Non-negative tolerance for accepting small
#'   decreases and declaring convergence.
#' @param tolerance Outer Tucker3 convergence tolerance.
#' @param max_iter Maximum Tucker3 iterations.
#' @param inner_max_iter Maximum HOOI iterations.
#' @param inner_tolerance HOOI tolerance.
#' @param nugget_upper Optional upper bound passed to nugget profiling.
#' @param nugget_tolerance Optimization tolerance for nugget profiling.
#' @param nugget_boundary_tolerance Boundary tolerance for tau = 0.
#' @param nugget_expansion_factor Automatic nugget-bound expansion factor.
#' @param nugget_max_expansions Maximum automatic nugget-bound expansions.
#' @param display Logical; emit refinement progress.
#' @return An object of class `scr_tucker3_nugget_refinement`.
#' @export
refine_scr_tucker3_nugget <- function(
  X,
  membership,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  refinement_rounds = 5L,
  objective_tolerance = 1e-6,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  nugget_upper = NULL,
  nugget_tolerance = 1e-6,
  nugget_boundary_tolerance = 1e-8,
  nugget_expansion_factor = 4,
  nugget_max_expansions = 6L,
  display = FALSE
) {
  validate_nugget_refinement_inputs(
    X = X,
    membership = membership,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    refinement_rounds = refinement_rounds,
    objective_tolerance = objective_tolerance,
    display = display
  )

  initial_fit <- fit_scr_s3_tucker3(
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
    display = FALSE
  )

  initial_profile <- profile_scr_tucker3_nugget(
    X = X,
    fit = initial_fit,
    upper = nugget_upper,
    tolerance = nugget_tolerance,
    boundary_tolerance = nugget_boundary_tolerance,
    expansion_factor = nugget_expansion_factor,
    max_expansions = nugget_max_expansions
  )

  accepted <- list(
    fit = initial_fit,
    profile = initial_profile,
    membership = initial_profile$membership
  )

  tau_trace <- initial_profile$nugget
  loglik_trace <- initial_profile$log_likelihood
  bic_trace <- initial_profile$bic
  icl_trace <- initial_profile$icl
  accepted_rounds <- 0L
  proposed_rounds <- 0L
  rollback <- NULL
  converged <- refinement_rounds == 0L

  if (refinement_rounds > 0L) {
    for (round in seq_len(refinement_rounds)) {
      proposed_rounds <- round

      if (display) {
        message(
          sprintf(
            "nugget refinement round %d/%d",
            round,
            refinement_rounds
          )
        )
      }

      proposal_fit <- fit_scr_s3_tucker3(
        X = X,
        membership = accepted$membership,
        centroid_rank = centroid_rank,
        variable_rank = variable_rank,
        occasion_rank = occasion_rank,
        variable_covariance = accepted$fit$SV,
        occasion_covariance = accepted$fit$SO,
        tolerance = tolerance,
        max_iter = max_iter,
        inner_max_iter = inner_max_iter,
        inner_tolerance = inner_tolerance,
        display = FALSE
      )

      proposal_profile <- profile_scr_tucker3_nugget(
        X = X,
        fit = proposal_fit,
        upper = nugget_upper,
        tolerance = nugget_tolerance,
        boundary_tolerance = nugget_boundary_tolerance,
        expansion_factor = nugget_expansion_factor,
        max_expansions = nugget_max_expansions
      )

      decision <- evaluate_nugget_refinement_step(
        previous_loglik = accepted$profile$log_likelihood,
        proposed_loglik = proposal_profile$log_likelihood,
        objective_tolerance = objective_tolerance
      )

      if (!decision$accept) {
        rollback <- list(
          round = as.integer(round),
          previous_log_likelihood = accepted$profile$log_likelihood,
          proposed_log_likelihood = proposal_profile$log_likelihood,
          decrease = decision$decrease
        )
        break
      }

      improvement <- proposal_profile$log_likelihood -
        accepted$profile$log_likelihood

      accepted <- list(
        fit = proposal_fit,
        profile = proposal_profile,
        membership = proposal_profile$membership
      )
      accepted_rounds <- accepted_rounds + 1L

      tau_trace <- c(tau_trace, proposal_profile$nugget)
      loglik_trace <- c(loglik_trace, proposal_profile$log_likelihood)
      bic_trace <- c(bic_trace, proposal_profile$bic)
      icl_trace <- c(icl_trace, proposal_profile$icl)

      if (improvement <= objective_tolerance) {
        converged <- TRUE
        break
      }
    }
  }

  structure(
    list(
      fit = accepted$fit,
      profile = accepted$profile,
      membership = accepted$membership,
      tau_trace = tau_trace,
      log_likelihood_trace = loglik_trace,
      bic_trace = bic_trace,
      icl_trace = icl_trace,
      accepted_rounds = as.integer(accepted_rounds),
      proposed_rounds = as.integer(proposed_rounds),
      converged = converged,
      rollback = rollback,
      objective_tolerance = objective_tolerance,
      ranks = c(
        P = as.integer(centroid_rank),
        Q = as.integer(variable_rank),
        R = as.integer(occasion_rank)
      )
    ),
    class = "scr_tucker3_nugget_refinement"
  )
}


evaluate_nugget_refinement_step <- function(
  previous_loglik,
  proposed_loglik,
  objective_tolerance
) {
  if (
    !is.numeric(previous_loglik) ||
      length(previous_loglik) != 1L ||
      !is.finite(previous_loglik) ||
      !is.numeric(proposed_loglik) ||
      length(proposed_loglik) != 1L ||
      !is.finite(proposed_loglik)
  ) {
    stop("refinement likelihood values must be finite scalars.")
  }

  if (
    !is.numeric(objective_tolerance) ||
      length(objective_tolerance) != 1L ||
      is.na(objective_tolerance) ||
      !is.finite(objective_tolerance) ||
      objective_tolerance < 0
  ) {
    stop("objective_tolerance must be a non-negative finite number.")
  }

  decrease <- previous_loglik - proposed_loglik

  list(
    accept = decrease <= objective_tolerance,
    decrease = max(decrease, 0)
  )
}


validate_nugget_refinement_inputs <- function(
  X,
  membership,
  centroid_rank,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  refinement_rounds,
  objective_tolerance,
  display
) {
  validate_scr_s3_tucker3_inputs(
    X = X,
    membership = membership,
    centroid_rank = centroid_rank,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = 1e-6,
    max_iter = 1L,
    inner_max_iter = 1L,
    inner_tolerance = 1e-8,
    display = FALSE
  )

  validate_model_selection_dimension(
    refinement_rounds,
    "refinement_rounds",
    minimum = 0L
  )

  if (
    !is.numeric(objective_tolerance) ||
      length(objective_tolerance) != 1L ||
      is.na(objective_tolerance) ||
      !is.finite(objective_tolerance) ||
      objective_tolerance < 0
  ) {
    stop("objective_tolerance must be a non-negative finite number.")
  }

  if (!is.logical(display) || length(display) != 1L || is.na(display)) {
    stop("display must be TRUE or FALSE.")
  }

  invisible(TRUE)
}
