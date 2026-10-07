test_that("paired S3/Tucker3 comparison fixes Q and R", {
  set.seed(2201)

  groups <- 3L
  variables <- 2L
  occasions <- 2L
  n <- 72L
  labels <- rep(seq_len(groups), each = n / groups)

  means <- rbind(
    c(-2, 0, -1, 0),
    c(0, 2, 0, 1),
    c(2, 0, 1, 0)
  )
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(n * variables * occasions, sd = 0.4), nrow = n)

  truth <- matrix(0, nrow = n, ncol = groups)
  truth[cbind(seq_len(n), labels)] <- 1

  membership <- truth * 0.9 + 0.1 / groups
  membership <- membership / rowSums(membership)

  result <- compare_scr_tucker3_s3(
    X = X,
    membership = membership,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(variables),
    occasion_covariance = diag(occasions),
    centroid_ranks = c(1, 2),
    truth = truth,
    max_iter = 8,
    inner_max_iter = 10
  )

  expect_s3_class(result, "scr_tucker3_s3_comparison")
  expect_equal(result$comparison$model, c("S3", "Tucker3"))
  expect_true(all(result$comparison$Q == 1L))
  expect_true(all(result$comparison$R == 1L))
  expect_true(result$comparison$P[2] %in% c(1L, 2L))
  expect_true(all(is.finite(result$comparison$bic)))
  expect_true(all(is.finite(result$comparison$ari)))
  expect_equal(sum(result$comparison$selected), 1L)
  expect_true(result$preferred_model %in% c("S3", "Tucker3"))
})


test_that("S3 benchmark parameter count matches the fitted BIC", {
  set.seed(2202)

  n <- 60L
  groups <- 2L
  variables <- 2L
  occasions <- 2L
  X <- matrix(rnorm(n * variables * occasions), nrow = n)
  membership <- matrix(runif(n * groups), nrow = n)
  membership <- membership / rowSums(membership)

  result <- compare_scr_tucker3_s3(
    X = X,
    membership = membership,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(variables),
    occasion_covariance = diag(occasions),
    centroid_ranks = 1,
    max_iter = 5,
    inner_max_iter = 5
  )

  s3 <- result$comparison[result$comparison$model == "S3", ]
  expected <- 2 * result$s3_fit$like -
    log(n) * s3$parameters

  expect_equal(result$s3_fit$bic, expected)
})


test_that("Tucker3 benchmark is reproducible under fixed scenario seeds", {
  first <- run_scr_tucker3_benchmark(
    n = 40,
    groups = 2,
    n_starts = 1,
    dgp = 1,
    n_simulations = 1,
    scenario = "I",
    source = "matlab",
    variable_rank = 1,
    occasion_rank = 1,
    centroid_ranks = 1,
    max_iter = 5,
    inner_max_iter = 5,
    keep_data = FALSE
  )

  second <- run_scr_tucker3_benchmark(
    n = 40,
    groups = 2,
    n_starts = 1,
    dgp = 1,
    n_simulations = 1,
    scenario = "I",
    source = "matlab",
    variable_rank = 1,
    occasion_rank = 1,
    centroid_ranks = 1,
    max_iter = 5,
    inner_max_iter = 5,
    keep_data = FALSE
  )

  expect_s3_class(first, "scr_tucker3_benchmark")
  expect_equal(first$ari, second$ari)
  expect_equal(first$bic, second$bic)
  expect_equal(first$log_likelihood, second$log_likelihood)
  expect_equal(first$selected_P, second$selected_P)
  expect_equal(first$preferred_model, second$preferred_model)
  expect_equal(nrow(first$summary), 2L)
})


test_that("existing SCR baseline simulation remains a three-model comparison", {
  result <- run_scr_simulation(
    n = 40,
    groups = 2,
    n_starts = 1,
    dgp = 1,
    n_simulations = 1,
    scenario = "I",
    source = "matlab",
    max_iter = 5,
    keep_data = FALSE
  )

  expect_equal(colnames(result$ari), c("S3", "S2", "H"))
})
