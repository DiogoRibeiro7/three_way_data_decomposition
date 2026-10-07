#' Compare S3 with Tucker3 Centroid Reduction
#'
#' Fit the validated S3 baseline and the Tucker3 centroid-mode extension to
#' the same data from the same initial membership matrix and covariance
#' factors. Variable and occasion ranks are held fixed, so the comparison
#' isolates the additional centroid-mode reduction introduced by Tucker3.
#'
#' Tucker3 selects only the centroid rank P by BIC. The S3 and Tucker3 fits use
#' the same BIC convention, `2 * logLik - log(n) * k`, so larger values are
#' preferred.
#'
#' @param X Numeric observation matrix with observations in rows.
#' @param membership Initial posterior-membership matrix shared by both models.
#' @param variable_rank Variable-mode rank Q used by both models.
#' @param occasion_rank Occasion-mode rank R used by both models.
#' @param variable_covariance Initial positive-definite variable covariance.
#' @param occasion_covariance Initial positive-definite occasion covariance.
#' @param centroid_ranks Candidate Tucker3 centroid ranks P. Defaults to all
#'   admissible ranks from 1 through G - 1.
#' @param variable_basis Optional initial S3 variable basis. When `NULL`, a
#'   deterministic covariance-normalized coordinate basis is used.
#' @param occasion_basis Optional initial S3 occasion basis. When `NULL`, a
#'   deterministic covariance-normalized coordinate basis is used.
#' @param truth Optional true membership matrix for ARI evaluation.
#' @param tolerance Positive outer convergence tolerance.
#' @param max_iter Positive maximum number of outer iterations.
#' @param inner_max_iter Maximum HOOI iterations for Tucker3 mean updates.
#' @param inner_tolerance Tucker3 inner-loop convergence tolerance.
#' @param bic_tolerance Absolute tolerance used to identify BIC ties.
#' @param display Logical; emit compact progress messages when `TRUE`.
#' @return An object of class `scr_tucker3_s3_comparison` with a two-row model
#'   table, both fitted models, the Tucker3 rank-selection object, and the
#'   BIC-preferred model.
#' @export
compare_scr_tucker3_s3 <- function(
  X,
  membership,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  centroid_ranks = seq_len(ncol(membership) - 1L),
  variable_basis = NULL,
  occasion_basis = NULL,
  truth = NULL,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  bic_tolerance = 1e-8,
  display = FALSE
) {
  validate_tucker3_s3_comparison_inputs(
    X = X,
    membership = membership,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    truth = truth,
    bic_tolerance = bic_tolerance,
    display = display
  )

  groups <- ncol(membership)
  variables <- nrow(variable_covariance)
  occasions <- nrow(occasion_covariance)

  centroid_ranks <- validate_rank_candidates(
    centroid_ranks,
    upper = groups - 1L,
    name = "centroid_ranks"
  )

  if (is.null(variable_basis)) {
    variable_basis <- covariance_normalize_basis(
      diag(variables)[, seq_len(variable_rank), drop = FALSE],
      variable_covariance
    )
  }

  if (is.null(occasion_basis)) {
    occasion_basis <- covariance_normalize_basis(
      diag(occasions)[, seq_len(occasion_rank), drop = FALSE],
      occasion_covariance
    )
  }

  validate_scr_s3_inputs(
    X = X,
    membership = membership,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = tolerance,
    max_iter = max_iter,
    display = display
  )

  s3_fit <- fit_scr_s3(
    X = X,
    membership = membership,
    variable_basis = variable_basis,
    occasion_basis = occasion_basis,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    tolerance = tolerance,
    max_iter = max_iter,
    display = FALSE
  )

  tucker3_selection <- select_scr_tucker3_model(
    X = X,
    membership = membership,
    variable_covariance = variable_covariance,
    occasion_covariance = occasion_covariance,
    centroid_ranks = centroid_ranks,
    variable_ranks = variable_rank,
    occasion_ranks = occasion_rank,
    tolerance = tolerance,
    max_iter = max_iter,
    inner_max_iter = inner_max_iter,
    inner_tolerance = inner_tolerance,
    bic_tolerance = bic_tolerance,
    display = FALSE
  )
  tucker3_fit <- tucker3_selection$selected_model
  selected_P <- unname(tucker3_selection$selected_ranks[["P"]])

  s3_parameters <- scr_s3_parameter_count_internal(
    groups = groups,
    variables = variables,
    occasions = occasions,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  tucker_row <- tucker3_selection$comparison[
    tucker3_selection$comparison$selected,
    ,
    drop = FALSE
  ]

  s3_ari <- NA_real_
  tucker3_ari <- NA_real_
  if (!is.null(truth)) {
    truth_hard <- hard_partition(truth)
    s3_ari <- adjusted_rand_index(
      crossprod(truth_hard, hard_partition(s3_fit$U))
    )
    tucker3_ari <- adjusted_rand_index(
      crossprod(truth_hard, hard_partition(tucker3_fit$U))
    )
  }

  comparison <- data.frame(
    model = c("S3", "Tucker3"),
    P = c(groups - 1L, selected_P),
    Q = c(variable_rank, variable_rank),
    R = c(occasion_rank, occasion_rank),
    parameters = c(s3_parameters, tucker_row$parameters),
    log_likelihood = c(s3_fit$like, tucker3_fit$like),
    bic = c(s3_fit$bic, tucker3_fit$bic),
    converged = c(
      isTRUE(s3_fit$converged),
      isTRUE(tucker3_fit$converged)
    ),
    inner_converged = c(NA, isTRUE(tucker3_fit$inner_converged)),
    iterations = c(
      as.integer(s3_fit$it),
      as.integer(tucker3_fit$it)
    ),
    ari = c(s3_ari, tucker3_ari),
    stringsAsFactors = FALSE
  )

  preferred_index <- select_preferred_benchmark_model(
    comparison,
    bic_tolerance = bic_tolerance
  )
  comparison$selected <- seq_len(nrow(comparison)) == preferred_index

  if (display) {
    message(
      sprintf(
        "S3 BIC=%g; Tucker3(P=%d) BIC=%g; preferred=%s",
        s3_fit$bic,
        selected_P,
        tucker3_fit$bic,
        comparison$model[preferred_index]
      )
    )
  }

  structure(
    list(
      comparison = comparison,
      preferred_model = comparison$model[preferred_index],
      s3_fit = s3_fit,
      tucker3_fit = tucker3_fit,
      tucker3_selection = tucker3_selection,
      bic_tolerance = bic_tolerance
    ),
    class = "scr_tucker3_s3_comparison"
  )
}


#' Run a Tucker3-versus-S3 Simulation Benchmark
#'
#' Reuse the SCR simulation scenarios to compare the centroid-reduced Tucker3
#' extension directly against S3. The historical S3/S2/H reproduction remains
#' unchanged; this function is a separate research-extension benchmark.
#'
#' For every simulated data set, both models are evaluated over the same set of
#' membership starts. Each model retains the start with the largest
#' log-likelihood, matching the existing multi-start simulation convention.
#'
#' @param n Number of observations per generated data set.
#' @param groups Number of mixture components.
#' @param n_starts Number of shared membership starts.
#' @param dgp Data-generating process, 1 through 4.
#' @param n_simulations Number of generated data sets.
#' @param scenario Either `"I"` or `"II"`.
#' @param source Either `"matlab"` or `"paper"`.
#' @param variable_rank Optional Q. Defaults to the scenario configuration.
#' @param occasion_rank Optional R. Defaults to the scenario configuration.
#' @param centroid_ranks Optional Tucker3 P candidates.
#' @param tolerance Outer convergence tolerance.
#' @param max_iter Maximum outer iterations.
#' @param inner_max_iter Maximum Tucker3 HOOI iterations.
#' @param inner_tolerance Tucker3 inner-loop tolerance.
#' @param bic_tolerance Absolute BIC tie tolerance.
#' @param preprocess Optional preprocessing hook, as in `run_scr_simulation()`.
#' @param membership_initializer Optional membership initializer.
#' @param keep_data Logical; retain generated observations and truth.
#' @return An object of class `scr_tucker3_benchmark` containing per-simulation
#'   model results, ARI values, BIC values, selected P values, and summaries.
#' @export
run_scr_tucker3_benchmark <- function(
  n,
  groups,
  n_starts,
  dgp,
  n_simulations = 100L,
  scenario = c("I", "II"),
  source = c("matlab", "paper"),
  variable_rank = NULL,
  occasion_rank = NULL,
  centroid_ranks = NULL,
  tolerance = 1e-6,
  max_iter = 1000L,
  inner_max_iter = 50L,
  inner_tolerance = 1e-8,
  bic_tolerance = 1e-8,
  preprocess = NULL,
  membership_initializer = NULL,
  keep_data = FALSE
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
  validate_bic_tolerance(bic_tolerance)

  scenario <- match.arg(scenario)
  source <- match.arg(source)
  config <- scr_simulation_config(scenario, source)

  if (is.null(variable_rank)) {
    variable_rank <- config$Q
  }
  if (is.null(occasion_rank)) {
    occasion_rank <- config$R
  }
  validate_model_selection_dimension(variable_rank, "variable_rank")
  validate_model_selection_dimension(occasion_rank, "occasion_rank")

  if (variable_rank > config$J) {
    stop("variable_rank cannot exceed the scenario variable dimension.")
  }
  if (occasion_rank > config$K) {
    stop("occasion_rank cannot exceed the scenario occasion dimension.")
  }

  if (is.null(centroid_ranks)) {
    centroid_ranks <- seq_len(groups - 1L)
  }
  centroid_ranks <- validate_rank_candidates(
    centroid_ranks,
    upper = groups - 1L,
    name = "centroid_ranks"
  )

  ari <- matrix(
    NA_real_,
    nrow = n_simulations,
    ncol = 2L,
    dimnames = list(NULL, c("S3", "Tucker3"))
  )
  bic <- matrix(
    NA_real_,
    nrow = n_simulations,
    ncol = 2L,
    dimnames = list(NULL, c("S3", "Tucker3"))
  )
  likelihood <- matrix(
    NA_real_,
    nrow = n_simulations,
    ncol = 2L,
    dimnames = list(NULL, c("S3", "Tucker3"))
  )
  selected_P <- rep(NA_integer_, n_simulations)
  preferred_model <- rep(NA_character_, n_simulations)
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

    best_s3 <- NULL
    best_tucker3 <- NULL

    for (start in seq_len(n_starts)) {
      seed <- 10L * simulation + start
      set.seed(seed)

      membership <- initialize_scr_membership(
        X = sample$X,
        groups = groups,
        start = start,
        seed = seed,
        initializer = membership_initializer
      )

      paired <- tryCatch(
        compare_scr_tucker3_s3(
          X = sample$X,
          membership = membership,
          variable_rank = variable_rank,
          occasion_rank = occasion_rank,
          variable_covariance = diag(config$J),
          occasion_covariance = diag(config$K),
          centroid_ranks = centroid_ranks,
          truth = sample$U_true,
          tolerance = tolerance,
          max_iter = max_iter,
          inner_max_iter = inner_max_iter,
          inner_tolerance = inner_tolerance,
          bic_tolerance = bic_tolerance,
          display = FALSE
        ),
        error = function(error) {
          structure(
            list(message = conditionMessage(error)),
            class = "scr_tucker3_benchmark_fit_error"
          )
        }
      )

      if (inherits(paired, "scr_tucker3_benchmark_fit_error")) {
        next
      }

      if (is.null(best_s3) || paired$s3_fit$like > best_s3$fit$like) {
        best_s3 <- list(
          fit = paired$s3_fit,
          start = start,
          ari = paired$comparison$ari[paired$comparison$model == "S3"]
        )
      }

      if (
        is.null(best_tucker3) ||
          paired$tucker3_fit$like > best_tucker3$fit$like
      ) {
        best_tucker3 <- list(
          fit = paired$tucker3_fit,
          start = start,
          ari = paired$comparison$ari[
            paired$comparison$model == "Tucker3"
          ],
          selection = paired$tucker3_selection
        )
      }
    }

    if (is.null(best_s3) || is.null(best_tucker3)) {
      stop(
        sprintf(
          "all paired S3/Tucker3 starts failed for simulation %d.",
          simulation
        )
      )
    }

    ari[simulation, ] <- c(best_s3$ari, best_tucker3$ari)
    bic[simulation, ] <- c(best_s3$fit$bic, best_tucker3$fit$bic)
    likelihood[simulation, ] <- c(
      best_s3$fit$like,
      best_tucker3$fit$like
    )
    selected_P[simulation] <- unname(
      best_tucker3$selection$selected_ranks[["P"]]
    )

    s3_parameters <- scr_s3_parameter_count_internal(
      groups = groups,
      variables = config$J,
      occasions = config$K,
      variable_rank = variable_rank,
      occasion_rank = occasion_rank
    )
    tucker_row <- best_tucker3$selection$comparison[
      best_tucker3$selection$comparison$selected,
      ,
      drop = FALSE
    ]

    model_table <- data.frame(
      model = c("S3", "Tucker3"),
      P = c(groups - 1L, selected_P[simulation]),
      Q = c(variable_rank, variable_rank),
      R = c(occasion_rank, occasion_rank),
      parameters = c(s3_parameters, tucker_row$parameters),
      bic = bic[simulation, ],
      stringsAsFactors = FALSE
    )
    preferred_index <- select_preferred_benchmark_model(
      model_table,
      bic_tolerance = bic_tolerance
    )
    preferred_model[simulation] <- model_table$model[preferred_index]

    diagnostics[[simulation]] <- list(
      S3 = list(
        converged = best_s3$fit$converged,
        iterations = best_s3$fit$it,
        start = best_s3$start
      ),
      Tucker3 = list(
        converged = best_tucker3$fit$converged,
        inner_converged = best_tucker3$fit$inner_converged,
        iterations = best_tucker3$fit$it,
        start = best_tucker3$start,
        selected_P = selected_P[simulation]
      )
    )

    if (keep_data) {
      generated_X[, , simulation] <- sample$X
      generated_truth[, , simulation] <- sample$U_true
    }
  }

  summary <- data.frame(
    model = c("S3", "Tucker3"),
    mean_ari = colMeans(ari),
    sd_ari = apply(
      ari,
      2L,
      function(values) {
        if (length(values) > 1L) stats::sd(values) else NA_real_
      }
    ),
    mean_bic = colMeans(bic),
    mean_log_likelihood = colMeans(likelihood),
    preferred_count = c(
      sum(preferred_model == "S3"),
      sum(preferred_model == "Tucker3")
    ),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      ari = ari,
      bic = bic,
      log_likelihood = likelihood,
      selected_P = selected_P,
      preferred_model = preferred_model,
      diagnostics = diagnostics,
      summary = summary,
      config = config,
      variable_rank = variable_rank,
      occasion_rank = occasion_rank,
      centroid_ranks = centroid_ranks,
      dgp = dgp,
      n = n,
      groups = groups,
      n_starts = n_starts,
      n_simulations = n_simulations,
      X = generated_X,
      U_true = generated_truth
    ),
    class = "scr_tucker3_benchmark"
  )
}


scr_s3_parameter_count_internal <- function(
  groups,
  variables,
  occasions,
  variable_rank,
  occasion_rank
) {
  validate_tucker3_dimensions(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = groups - 1L,
    variable_rank = variable_rank,
    occasion_rank = occasion_rank
  )

  as.integer(
    groups - 1L +
      variables * occasions +
      (groups - 1L) * variable_rank * occasion_rank +
      (variables - variable_rank) * variable_rank +
      (occasions - occasion_rank) * occasion_rank +
      variables * (variables + 1L) / 2L +
      occasions * (occasions + 1L) / 2L -
      1L
  )
}


select_preferred_benchmark_model <- function(
  comparison,
  bic_tolerance
) {
  required <- c("model", "parameters", "bic")
  if (
    !is.data.frame(comparison) ||
      nrow(comparison) < 1L ||
      !all(required %in% names(comparison))
  ) {
    stop(
      "comparison must contain model, parameters, and bic columns."
    )
  }
  validate_bic_tolerance(bic_tolerance)

  if (
    anyNA(comparison$bic) ||
      any(!is.finite(comparison$bic)) ||
      anyNA(comparison$parameters) ||
      any(!is.finite(comparison$parameters))
  ) {
    stop("benchmark comparison statistics must be finite.")
  }

  best_bic <- max(comparison$bic)
  tied <- which(best_bic - comparison$bic <= bic_tolerance)
  tied[order(comparison$parameters[tied], comparison$model[tied])[1L]]
}


validate_tucker3_s3_comparison_inputs <- function(
  X,
  membership,
  variable_rank,
  occasion_rank,
  variable_covariance,
  occasion_covariance,
  truth,
  bic_tolerance,
  display
) {
  if (
    !is.matrix(X) ||
      !is.numeric(X) ||
      nrow(X) < 2L ||
      ncol(X) < 4L ||
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
      "membership must be a finite non-negative n-by-G numeric matrix."
    )
  }

  if (max(abs(rowSums(membership) - 1)) > 1e-10) {
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

  variables <- nrow(variable_covariance)
  occasions <- nrow(occasion_covariance)

  if (ncol(X) != variables * occasions) {
    stop(
      "X columns must equal variables times occasions implied by the covariance factors."
    )
  }

  validate_model_selection_dimension(variable_rank, "variable_rank")
  validate_model_selection_dimension(occasion_rank, "occasion_rank")

  if (variable_rank > variables) {
    stop("variable_rank cannot exceed variables.")
  }
  if (occasion_rank > occasions) {
    stop("occasion_rank cannot exceed occasions.")
  }

  if (!is.null(truth)) {
    if (
      !is.matrix(truth) ||
        !is.numeric(truth) ||
        !identical(dim(truth), dim(membership)) ||
        anyNA(truth) ||
        any(!is.finite(truth)) ||
        any(truth < 0)
    ) {
      stop(
        "truth must be a finite non-negative membership matrix matching membership."
      )
    }
    if (any(rowSums(truth) <= 0)) {
      stop("truth must assign positive mass to every observation.")
    }
  }

  validate_bic_tolerance(bic_tolerance)

  if (!is.logical(display) || length(display) != 1L || is.na(display)) {
    stop("display must be TRUE or FALSE.")
  }

  invisible(TRUE)
}
