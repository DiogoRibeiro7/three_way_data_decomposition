test_that("Tucker3 stability is deterministic under a fixed seed", {
  set.seed(2801)
  labels <- rep(1:2, each = 25)
  means <- rbind(c(-2, 0, -1, 0), c(2, 0, 1, 0))
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(50 * 4, sd = 0.35), nrow = 50)

  first <- scr_tucker3_stability(
    X = X,
    groups = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_resamples = 3,
    subsample_fraction = 0.8,
    n_starts = 1,
    seed = 42,
    max_iter = 4,
    inner_max_iter = 4
  )

  second <- scr_tucker3_stability(
    X = X,
    groups = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_resamples = 3,
    subsample_fraction = 0.8,
    n_starts = 1,
    seed = 42,
    max_iter = 4,
    inner_max_iter = 4
  )

  expect_s3_class(first, "scr_tucker3_stability")
  expect_equal(first$mean_ari, second$mean_ari)
  expect_equal(first$pairwise, second$pairwise)
  expect_equal(
    lapply(first$resamples, function(x) x$index),
    lapply(second$resamples, function(x) x$index)
  )
})


test_that("ARI stability comparison is invariant to label permutations", {
  partition_a <- rbind(
    c(1, 0),
    c(1, 0),
    c(0, 1),
    c(0, 1)
  )
  partition_b <- partition_a[, 2:1, drop = FALSE]

  contingency <- crossprod(partition_a, partition_b)

  expect_equal(adjusted_rand_index(contingency), 1)
})


test_that("fixed-rank stability records expected diagnostics", {
  set.seed(2802)
  labels <- rep(1:2, each = 24)
  means <- rbind(c(-2, 0, -1, 0), c(2, 0, 1, 0))
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(48 * 4, sd = 0.35), nrow = 48)

  result <- scr_tucker3_stability(
    X = X,
    groups = 2,
    centroid_rank = 1,
    variable_rank = 1,
    occasion_rank = 1,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    n_resamples = 3,
    subsample_fraction = 0.75,
    n_starts = 1,
    seed = 7,
    max_iter = 4,
    inner_max_iter = 4
  )

  expect_true(is.finite(result$mean_ari))
  expect_true(result$n_pairs >= 1L)
  expect_equal(result$n_success + result$n_failed, 3L)
  expect_equal(result$resample_size, floor(48 * 0.75))
  expect_equal(
    unname(result$structure),
    c(2L, 1L, 1L, 1L)
  )
})


test_that("stability selector uses lower complexity for ties", {
  comparison <- data.frame(
    P = c(2L, 1L, 1L),
    Q = c(1L, 1L, 2L),
    R = c(1L, 1L, 1L),
    parameters = c(25L, 20L, 22L),
    mean_ari = c(0.8, 0.8 + 5e-9, 0.7)
  )

  expect_equal(
    threeway:::select_best_stability_candidate(
      comparison,
      stability_tolerance = 1e-8
    ),
    2L
  )
})


test_that("stability selector evaluates candidate ranks", {
  set.seed(2803)
  labels <- rep(1:3, each = 15)
  means <- rbind(
    c(-2, 0, -1, 0),
    c(0, 2, 0, 1),
    c(2, 0, 1, 0)
  )
  X <- means[labels, , drop = FALSE] +
    matrix(rnorm(45 * 4, sd = 0.35), nrow = 45)

  result <- select_scr_tucker3_stable(
    X = X,
    groups = 3,
    variable_covariance = diag(2),
    occasion_covariance = diag(2),
    centroid_ranks = c(1, 2),
    variable_ranks = 1,
    occasion_ranks = 1,
    n_resamples = 3,
    subsample_fraction = 0.8,
    n_starts = 1,
    seed = 3,
    max_iter = 3,
    inner_max_iter = 3
  )

  expect_s3_class(result, "scr_tucker3_stability_selection")
  expect_equal(result$comparison$P, c(1L, 2L))
  expect_equal(sum(result$comparison$selected), 1L)
  expect_true(all(is.finite(result$comparison$mean_ari)))
})


test_that("stability controls reject invalid values", {
  X <- matrix(rnorm(40), nrow = 10, ncol = 4)

  expect_error(
    scr_tucker3_stability(
      X = X,
      groups = 2,
      centroid_rank = 1,
      variable_rank = 1,
      occasion_rank = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2),
      n_resamples = 1
    ),
    "n_resamples"
  )

  expect_error(
    scr_tucker3_stability(
      X = X,
      groups = 2,
      centroid_rank = 1,
      variable_rank = 1,
      occasion_rank = 1,
      variable_covariance = diag(2),
      occasion_covariance = diag(2),
      subsample_fraction = 0
    ),
    "subsample_fraction"
  )
})
