test_that("SCR experiment grid reproduces the six legacy driver designs", {
  grid_i <- scr_experiment_grid(
    scenario = "I",
    source = "matlab"
  )
  grid_ii <- scr_experiment_grid(
    scenario = "II",
    source = "matlab"
  )

  expect_equal(nrow(grid_i), 48L)
  expect_equal(nrow(grid_ii), 48L)

  expect_equal(sort(unique(grid_i$groups)), c(3L, 5L, 7L))
  expect_equal(sort(unique(grid_i$n)), c(300L, 500L))
  expect_equal(sort(unique(grid_ii$n)), c(500L, 1000L))
  expect_equal(sort(unique(grid_i$n_starts)), c(1L, 3L))
  expect_equal(sort(unique(grid_i$dgp)), 1:4)
  expect_true(all(grid_i$n_simulations == 250L))

  expect_equal(
    grid_i$setting_id[grid_i$groups == 3L],
    seq_len(16)
  )
})

test_that("experiment IDs are stable and unique", {
  grid <- scr_experiment_grid(
    scenario = c("I", "II"),
    source = "paper",
    groups = c(3L, 5L)
  )

  expect_false(anyDuplicated(grid$experiment_id))
  expect_true(
    "I-G3-N300-S1-D1" %in% grid$experiment_id
  )
  expect_true(
    "II-G5-N1000-S3-D4" %in% grid$experiment_id
  )
})

test_that("experiment grid runner returns tidy model summaries", {
  grid <- scr_experiment_grid(
    scenario = "I",
    source = "matlab",
    groups = 2L,
    sample_sizes = 40L,
    n_starts = 1L,
    dgps = 1L,
    n_simulations = 1L
  )

  result <- run_scr_experiment_grid(
    grid = grid,
    max_iter = 5L,
    keep_data = FALSE,
    progress = FALSE
  )

  expect_s3_class(result, "scr_experiment_results")
  expect_equal(length(result$runs), 1L)
  expect_equal(nrow(result$summary), 3L)
  expect_equal(result$summary$model, c("S3", "S2", "H"))
  expect_true(all(is.finite(result$summary$mean_ari)))
  expect_equal(result$summary$n_requested, rep(1L, 3L))
  expect_equal(scr_experiment_summary(result), result$summary)
})

test_that("failed experiment cells can be retained explicitly", {
  grid <- scr_experiment_grid(
    scenario = "I",
    source = "matlab",
    groups = 2L,
    sample_sizes = 2L,
    n_starts = 1L,
    dgps = 1L,
    n_simulations = 1L
  )

  result <- run_scr_experiment_grid(
    grid = grid,
    max_iter = 2L,
    keep_data = FALSE,
    continue_on_error = TRUE,
    progress = FALSE
  )

  expect_equal(nrow(result$summary), 3L)
  expect_true(all(result$summary$n_success == 0L))
  expect_true(all(!is.na(result$summary$error)))
})

test_that("experiment grid validates duplicate identifiers", {
  grid <- scr_experiment_grid(
    scenario = "I",
    groups = 3L,
    sample_sizes = 30L,
    n_starts = 1L,
    dgps = 1L,
    n_simulations = 1L
  )

  duplicate <- rbind(grid, grid)

  expect_error(
    run_scr_experiment_grid(
      duplicate,
      progress = FALSE
    ),
    "must be unique"
  )
})
