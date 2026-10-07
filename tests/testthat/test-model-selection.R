test_that("Tucker3 rank grid is deterministic and structurally valid", {
  grid <- scr_tucker3_rank_grid(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_ranks = c(2, 1, 2),
    variable_ranks = c(3, 1),
    occasion_ranks = c(2, 1)
  )

  expect_equal(nrow(grid), 8L)
  expect_equal(grid$P, c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L))
  expect_true(all(grid$P <= 2L))
  expect_true(all(grid$Q <= 3L))
  expect_true(all(grid$R <= 2L))
  expect_true(all(grid$parameters > 0L))

  expect_error(
    scr_tucker3_rank_grid(
      groups = 3,
      variables = 3,
      occasions = 2,
      centroid_ranks = 3
    ),
    "centroid_ranks cannot contain values greater than 2"
  )
  expect_error(
    scr_tucker3_rank_grid(
      groups = 3,
      variables = 3,
      occasions = 2,
      variable_ranks = 4
    ),
    "variable_ranks cannot contain values greater than 3"
  )
  expect_error(
    scr_tucker3_rank_grid(
      groups = 3,
      variables = 3,
      occasions = 2,
      occasion_ranks = 3
    ),
    "occasion_ranks cannot contain values greater than 2"
  )
})


test_that("Tucker3 BIC ties prefer the lower-dimensional model", {
  comparison <- data.frame(
    P = c(1L, 2L, 1L),
    Q = c(1L, 1L, 2L),
    R = c(1L, 1L, 1L),
    parameters = c(20L, 25L, 22L),
    bic = c(100, 100 + 5e-9, 99)
  )

  selected <- threeway:::select_best_tucker3_candidate(
    comparison,
    bic_tolerance = 1e-8
  )

  expect_equal(selected, 1L)
})


test_that("Tucker3 model selection evaluates every candidate consistently", {
  set.seed(2020)

  groups <- 3L
  variables <- 2L
  occasions <- 2L
  observations_per_group <- 24L
  labels <- rep(seq_len(groups), each = observations_per_group)
  n <- length(labels)

  means <- rbind(
    c(-2, 0, -1, 0),
    c(0, 2, 0, 1),
    c(2, 0, 1, 0)
  )
  X <- means[labels, , drop = FALSE] +
    matrix(
      rnorm(n * variables * occasions, sd = 0.45),
      nrow = n,
      ncol = variables * occasions
    )

  membership <- matrix(
    0.025,
    nrow = n,
    ncol = groups
  )
  membership[cbind(seq_len(n), labels)] <- 0.95

  selection <- select_scr_tucker3_model(
    X = X,
    membership = membership,
    variable_covariance = diag(variables),
    occasion_covariance = diag(occasions),
    centroid_ranks = c(1, 2),
    variable_ranks = 1,
    occasion_ranks = 1,
    tolerance = 1e-5,
    max_iter = 8,
    inner_max_iter = 10,
    inner_tolerance = 1e-7,
    bic_tolerance = 1e-8
  )

  expect_s3_class(selection, "scr_tucker3_model_selection")
  expect_equal(nrow(selection$comparison), 2L)
  expect_equal(selection$comparison$P, c(1L, 2L))
  expect_true(all(is.finite(selection$comparison$log_likelihood)))
  expect_true(all(is.finite(selection$comparison$bic)))
  expect_equal(sum(selection$comparison$selected), 1L)

  selected_row <- selection$comparison[
    selection$comparison$selected,
    ,
    drop = FALSE
  ]

  expect_equal(
    unname(selection$selected_ranks),
    unname(as.integer(selected_row[1L, c("P", "Q", "R")]))
  )
  expect_equal(selection$selected_model$bic, selected_row$bic)
  expect_equal(selection$selected_model$like, selected_row$log_likelihood)
  expect_equal(
    isTRUE(selection$selected_model$converged),
    selected_row$converged
  )
  expect_equal(
    isTRUE(selection$selected_model$inner_converged),
    selected_row$inner_converged
  )
  expect_equal(
    as.integer(selection$selected_model$it),
    selected_row$iterations
  )

  expected_index <- threeway:::select_best_tucker3_candidate(
    selection$comparison,
    bic_tolerance = selection$bic_tolerance
  )
  expect_true(selection$comparison$selected[expected_index])
})
