#' Build the SCR Experiment Grid
#'
#' Construct the factorial simulation design used by the Rocci-Vichi-Ranalli
#' study. The historical MATLAB scripts split this design across six files for
#' G = 3, 5, and 7 and two scenarios. Here the factors are represented
#' explicitly in one data frame.
#'
#' @param scenario One or both of `"I"` and `"II"`.
#' @param source Either `"matlab"` or `"paper"`; only Scenario I differs
#'   between these sources.
#' @param groups Integer vector of component counts.
#' @param sample_sizes Optional integer vector. When `NULL`, Scenario I uses
#'   300 and 500, while Scenario II uses 500 and 1000.
#' @param n_starts Integer vector of random-start counts.
#' @param dgps Integer vector of data-generating processes.
#' @param n_simulations Number of simulated data sets per grid cell.
#' @return A data frame with one row per experiment cell.
#' @export
scr_experiment_grid <- function(
  scenario = c("I", "II"),
  source = c("matlab", "paper"),
  groups = c(3L, 5L, 7L),
  sample_sizes = NULL,
  n_starts = c(1L, 3L),
  dgps = 1:4,
  n_simulations = 250L
) {
  source <- match.arg(source)
  scenario <- unique(match.arg(scenario, choices = c("I", "II"), several.ok = TRUE))

  validate_positive_integer_vector(groups, "groups")
  if (any(groups < 2L)) {
    stop("groups must contain integers greater than or equal to two.")
  }

  validate_positive_integer_vector(n_starts, "n_starts")
  validate_positive_integer_vector(dgps, "dgps")
  if (any(!(dgps %in% 1:4))) {
    stop("dgps must contain only 1, 2, 3, and 4.")
  }

  validate_positive_integer_vector(n_simulations, "n_simulations")
  if (length(n_simulations) != 1L) {
    stop("n_simulations must be a single positive integer.")
  }

  rows <- list()
  row_index <- 1L

  for (scenario_value in scenario) {
    sizes <- sample_sizes
    if (is.null(sizes)) {
      sizes <- if (scenario_value == "I") c(300L, 500L) else c(500L, 1000L)
    }
    validate_positive_integer_vector(sizes, "sample_sizes")

    design <- expand.grid(
      n = as.integer(sizes),
      n_starts = as.integer(n_starts),
      dgp = as.integer(dgps),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    # Match the historical script ordering: sample size changes slowest,
    # starts next, DGP fastest.
    design <- design[
      order(
        match(design$n, sizes),
        match(design$n_starts, n_starts),
        match(design$dgp, dgps)
      ),
      ,
      drop = FALSE
    ]
    rownames(design) <- NULL

    for (group_value in groups) {
      block <- data.frame(
        scenario = scenario_value,
        source = source,
        groups = as.integer(group_value),
        n = design$n,
        n_starts = design$n_starts,
        dgp = design$dgp,
        n_simulations = as.integer(n_simulations),
        stringsAsFactors = FALSE
      )

      block$setting_id <- seq_len(nrow(block))
      block$experiment_id <- sprintf(
        "%s-G%d-N%d-S%d-D%d",
        scenario_value,
        group_value,
        block$n,
        block$n_starts,
        block$dgp
      )

      rows[[row_index]] <- block
      row_index <- row_index + 1L
    }
  }

  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid
}


#' Run an SCR Experiment Grid
#'
#' Evaluate every row of an experiment grid using `run_scr_simulation()`.
#'
#' @param grid Data frame returned by `scr_experiment_grid()`, or a compatible
#'   data frame containing the required columns.
#' @param tolerance Convergence threshold passed to each model fit.
#' @param max_iter Maximum iterations passed to each model fit.
#' @param preprocess Optional preprocessing hook passed to
#'   `run_scr_simulation()`.
#' @param membership_initializer Optional membership-initialization hook.
#' @param keep_data Logical; retain generated observations inside every result.
#' @param continue_on_error Logical; when `TRUE`, preserve failed cells as
#'   error records and continue with the remaining grid.
#' @param progress Logical; emit one message per grid cell.
#' @return An object of class `scr_experiment_results` containing the grid,
#'   cell-level runs, and one tidy ARI summary row per model and grid cell.
#' @export
run_scr_experiment_grid <- function(
  grid,
  tolerance = 1e-6,
  max_iter = 1000L,
  preprocess = NULL,
  membership_initializer = NULL,
  keep_data = FALSE,
  continue_on_error = TRUE,
  progress = TRUE
) {
  validate_scr_experiment_grid(grid)

  if (!is.logical(continue_on_error) ||
      length(continue_on_error) != 1L ||
      is.na(continue_on_error)) {
    stop("continue_on_error must be TRUE or FALSE.")
  }
  if (!is.logical(progress) || length(progress) != 1L || is.na(progress)) {
    stop("progress must be TRUE or FALSE.")
  }

  runs <- vector("list", nrow(grid))
  summary_rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    setting <- grid[i, , drop = FALSE]

    if (progress) {
      message(
        sprintf(
          "[%d/%d] %s",
          i,
          nrow(grid),
          setting$experiment_id
        )
      )
    }

    run <- tryCatch(
      run_scr_simulation(
        n = setting$n,
        groups = setting$groups,
        n_starts = setting$n_starts,
        dgp = setting$dgp,
        n_simulations = setting$n_simulations,
        scenario = setting$scenario,
        source = setting$source,
        tolerance = tolerance,
        max_iter = max_iter,
        preprocess = preprocess,
        membership_initializer = membership_initializer,
        keep_data = keep_data
      ),
      error = function(error) {
        if (!continue_on_error) {
          stop(error)
        }

        structure(
          list(message = conditionMessage(error)),
          class = "scr_experiment_error"
        )
      }
    )

    runs[[i]] <- run

    if (inherits(run, "scr_experiment_error")) {
      summary_rows[[i]] <- data.frame(
        experiment_id = setting$experiment_id,
        scenario = setting$scenario,
        source = setting$source,
        groups = setting$groups,
        n = setting$n,
        n_starts = setting$n_starts,
        dgp = setting$dgp,
        model = c("S3", "S2", "H"),
        mean_ari = NA_real_,
        sd_ari = NA_real_,
        median_ari = NA_real_,
        q25_ari = NA_real_,
        q75_ari = NA_real_,
        n_success = 0L,
        n_requested = setting$n_simulations,
        error = run$message,
        stringsAsFactors = FALSE
      )
      next
    }

    summary_rows[[i]] <- summarize_scr_ari(
      run = run,
      setting = setting
    )
  }

  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL

  structure(
    list(
      grid = grid,
      runs = runs,
      summary = summary
    ),
    class = "scr_experiment_results"
  )
}


#' Summarize an SCR Experiment Grid
#'
#' Extract the tidy model-by-setting ARI summary from an experiment result.
#'
#' @param x An object returned by `run_scr_experiment_grid()`.
#' @return A data frame with ARI summaries.
#' @export
scr_experiment_summary <- function(x) {
  if (!inherits(x, "scr_experiment_results")) {
    stop("x must be an scr_experiment_results object.")
  }

  x$summary
}


summarize_scr_ari <- function(run, setting) {
  models <- colnames(run$ari)

  do.call(
    rbind,
    lapply(
      models,
      function(model) {
        values <- run$ari[, model]
        finite <- values[is.finite(values)]

        data.frame(
          experiment_id = setting$experiment_id,
          scenario = setting$scenario,
          source = setting$source,
          groups = setting$groups,
          n = setting$n,
          n_starts = setting$n_starts,
          dgp = setting$dgp,
          model = model,
          mean_ari = if (length(finite)) mean(finite) else NA_real_,
          sd_ari = if (length(finite) > 1L) stats::sd(finite) else NA_real_,
          median_ari = if (length(finite)) stats::median(finite) else NA_real_,
          q25_ari = if (length(finite)) {
            unname(stats::quantile(finite, 0.25, names = FALSE))
          } else {
            NA_real_
          },
          q75_ari = if (length(finite)) {
            unname(stats::quantile(finite, 0.75, names = FALSE))
          } else {
            NA_real_
          },
          n_success = length(finite),
          n_requested = nrow(run$ari),
          error = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    )
  )
}


validate_scr_experiment_grid <- function(grid) {
  required <- c(
    "scenario",
    "source",
    "groups",
    "n",
    "n_starts",
    "dgp",
    "n_simulations",
    "experiment_id"
  )

  if (!is.data.frame(grid) || nrow(grid) == 0L) {
    stop("grid must be a non-empty data frame.")
  }

  missing <- setdiff(required, names(grid))
  if (length(missing)) {
    stop(
      paste(
        "grid is missing required columns:",
        paste(missing, collapse = ", ")
      )
    )
  }

  if (any(!(grid$scenario %in% c("I", "II")))) {
    stop("grid scenario values must be I or II.")
  }
  if (any(!(grid$source %in% c("matlab", "paper")))) {
    stop("grid source values must be matlab or paper.")
  }

  validate_positive_integer_vector(grid$groups, "grid$groups")
  validate_positive_integer_vector(grid$n, "grid$n")
  validate_positive_integer_vector(grid$n_starts, "grid$n_starts")
  validate_positive_integer_vector(grid$dgp, "grid$dgp")
  validate_positive_integer_vector(
    grid$n_simulations,
    "grid$n_simulations"
  )

  if (any(grid$groups < 2L)) {
    stop("grid$groups values must be at least two.")
  }
  if (any(!(grid$dgp %in% 1:4))) {
    stop("grid$dgp values must be between 1 and 4.")
  }

  if (
    !is.character(grid$experiment_id) ||
      anyNA(grid$experiment_id) ||
      any(grid$experiment_id == "")
  ) {
    stop("grid$experiment_id must contain non-empty strings.")
  }
  if (anyDuplicated(grid$experiment_id)) {
    stop("grid$experiment_id values must be unique.")
  }

  invisible(TRUE)
}


validate_positive_integer_vector <- function(x, name) {
  if (
    !is.numeric(x) ||
      length(x) == 0L ||
      anyNA(x) ||
      any(!is.finite(x)) ||
      any(x <= 0) ||
      any(x %% 1 != 0)
  ) {
    stop(paste(name, "must contain positive integers."))
  }

  invisible(TRUE)
}
