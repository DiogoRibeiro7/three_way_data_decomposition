test_that("full-support sparse parameter count matches Tucker3 count", {
  sparse <- scr_sparse_tucker3_parameter_count(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    active_variables = 3,
    active_occasions = 2
  )

  dense <- scr_tucker3_parameter_count(
    groups = 3,
    variables = 3,
    occasions = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1
  )

  expect_equal(sparse, dense)
})


test_that("smaller sparse supports reduce active-set parameter count", {
  full <- scr_sparse_tucker3_parameter_count(
    groups = 3,
    variables = 4,
    occasions = 3,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    active_variables = 4,
    active_occasions = 3
  )

  sparse <- scr_sparse_tucker3_parameter_count(
    groups = 3,
    variables = 4,
    occasions = 3,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    active_variables = 2,
    active_occasions = 1
  )

  expect_lt(sparse, full)
})


make_exact_sparse_selection_fixture <- function() {
  probabilities <- c(0.3, 0.3, 0.4)
  group_mass <- c(30, 30, 40)

  raw_A <- matrix(c(1, -1, 0), ncol = 1)
  A <- scr_centroid_basis(raw_A, probabilities)
  A <- A / sqrt(
    as.numeric(crossprod(A, group_mass * A))
  )

  B <- matrix(c(1, 0.2, 0), ncol = 1)
  B <- B / sqrt(sum(B^2))

  C <- matrix(c(1, 0.2), ncol = 1)
  C <- C / sqrt(sum(C^2))

  core <- array(3, dim = c(1, 1, 1))
  grand <- matrix(0, nrow = 3, ncol = 2)

  means <- scr_tucker3_means(
    grand_mean = grand,
    centroid_basis = A,
    variable_basis = B,
    occasion_basis = C,
    core = core,
    probabilities = probabilities
  )

  labels <- rep(1:3, times = group_mass)
  X <- means[labels, , drop = FALSE]

  U <- matrix(0, nrow = nrow(X), ncol = 3)
  U[cbind(seq_len(nrow(X)), labels)] <- 1

  fit <- list(
    U = U,
    A = A,
    TB = B,
    TC = C,
    SV = diag(3),
    SO = diag(2),
    probabilities = probabilities,
    ranks = c(P = 1L, Q = 1L, R = 1L)
  )

  list(X = X, fit = fit, means = means)
}


test_that("zero sparse penalties reproduce an exact fitted mean structure", {
  fixture <- make_exact_sparse_selection_fixture()

  result <- select_scr_sparse_tucker3(
    X = fixture$X,
    fit = fixture$fit,
    variable_penalties = 0,
    occasion_penalties = 0
  )

  expect_equal(
    result$selected_projection$means,
    fixture$means,
    tolerance = 1e-9
  )
  expect_equal(
    result$selected_penalties,
    c(variable = 0, occasion = 0)
  )
})


test_that("sparse penalty selection scores feasible cells by BIC and ICL", {
  fixture <- make_exact_sparse_selection_fixture()

  result <- select_scr_sparse_tucker3(
    X = fixture$X,
    fit = fixture$fit,
    variable_penalties = c(0, 0.25, 2),
    occasion_penalties = c(0, 0.25, 2),
    criterion = "BIC"
  )

  feasible <- result$comparison$feasible

  expect_true(any(feasible))
  expect_true(any(!feasible))
  expect_true(all(is.finite(result$comparison$bic[feasible])))
  expect_true(all(is.finite(result$comparison$icl[feasible])))
  expect_true(all(is.finite(result$comparison$parameters[feasible])))
  expect_false(any(result$comparison$selected[!feasible]))
  expect_equal(sum(result$comparison$selected), 1L)

  expect_equal(
    rowSums(result$selected_membership),
    rep(1, nrow(fixture$X)),
    tolerance = 1e-12
  )
})


test_that("sparse selector supports ICL explicitly", {
  fixture <- make_exact_sparse_selection_fixture()

  result <- select_scr_sparse_tucker3(
    X = fixture$X,
    fit = fixture$fit,
    variable_penalties = c(0, 0.25),
    occasion_penalties = c(0, 0.25),
    criterion = "ICL"
  )

  expect_equal(result$criterion, "ICL")
  expect_true(is.finite(result$criterion_value))
  expect_equal(sum(result$comparison$selected), 1L)
})


test_that("sparse selector tie-breaking prefers lower active complexity", {
  comparison <- data.frame(
    parameters = c(20L, 18L, 18L),
    n_active_variables = c(3L, 2L, 2L),
    n_active_occasions = c(2L, 2L, 1L),
    variable_penalty = c(0, 0.1, 0.2),
    occasion_penalty = c(0, 0.1, 0.2),
    bic = c(100, 100 + 5e-9, 100 + 4e-9),
    icl = c(90, 90, 90)
  )

  selected <- scr3way:::select_best_sparse_candidate(
    comparison,
    criterion = "BIC",
    criterion_tolerance = 1e-8
  )

  expect_equal(selected, 3L)
})


test_that("sparse parameter count validates impossible active sets", {
  expect_error(
    scr_sparse_tucker3_parameter_count(
      groups = 3,
      variables = 4,
      occasions = 3,
      centroid_rank = 1,
      variable_rank = 2,
      occasion_rank = 1,
      active_variables = 1,
      active_occasions = 2
    ),
    "active_variables cannot be smaller"
  )
})
