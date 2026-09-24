test_that("centroid basis is orthonormal and probability centered", {
  probabilities <- c(0.2, 0.3, 0.5)
  raw <- matrix(
    c(
      1, 0,
      0, 1,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE
  )

  A <- scr_centroid_basis(raw, probabilities)

  expect_equal(crossprod(A), diag(2), tolerance = 1e-10)
  expect_equal(
    as.numeric(crossprod(probabilities, A)),
    c(0, 0),
    tolerance = 1e-10
  )
})

test_that("Tucker3 means match the explicit index formula", {
  probabilities <- c(0.25, 0.35, 0.40)
  A <- scr_centroid_basis(
    matrix(c(1, 0, -1), ncol = 1),
    probabilities
  )
  B <- matrix(c(1, 2), ncol = 1)
  C <- matrix(c(3, 4), ncol = 1)
  core <- array(2, dim = c(1, 1, 1))
  grand <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    ncol = 2
  )

  result <- scr_tucker3_means(
    grand_mean = grand,
    centroid_basis = A,
    variable_basis = B,
    occasion_basis = C,
    core = core,
    probabilities = probabilities
  )

  expected <- matrix(0, nrow = 3, ncol = 4)
  for (g in 1:3) {
    index <- 1L
    for (k in 1:2) {
      for (j in 1:2) {
        expected[g, index] <- grand[j, k] +
          A[g, 1] * B[j, 1] * C[k, 1] * core[1, 1, 1]
        index <- index + 1L
      }
    }
  }

  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("probability-weighted Tucker3 mean equals the grand mean", {
  probabilities <- c(0.2, 0.3, 0.5)
  A <- scr_centroid_basis(
    matrix(
      c(
        1, 0,
        0, 1,
        1, 1
      ),
      nrow = 3,
      byrow = TRUE
    ),
    probabilities
  )

  B <- diag(2)
  C <- matrix(c(1, 0), ncol = 1)
  core <- array(
    c(0.5, -0.2, 0.1, 0.7),
    dim = c(2, 2, 1)
  )
  grand <- matrix(c(1, 2, 3, 4), nrow = 2)

  means <- scr_tucker3_means(
    grand,
    A,
    B,
    C,
    core,
    probabilities
  )

  expect_equal(
    as.numeric(crossprod(probabilities, means)),
    as.vector(grand),
    tolerance = 1e-10
  )
})

test_that("Tucker3 parameter count reduces to published S3 at full centroid rank", {
  groups <- 5L
  variables <- 8L
  occasions <- 4L
  q <- 3L
  r <- 2L
  p <- groups - 1L

  tucker3 <- scr_tucker3_parameter_count(
    groups = groups,
    variables = variables,
    occasions = occasions,
    centroid_rank = p,
    variable_rank = q,
    occasion_rank = r
  )

  published_s3 <- groups - 1L +
    variables * occasions +
    (groups - 1L) * q * r +
    (variables - q) * q +
    (occasions - r) * r +
    variables * (variables + 1L) / 2L +
    occasions * (occasions + 1L) / 2L -
    1L

  expect_equal(tucker3, as.integer(published_s3))
})

test_that("centroid reduction lowers the model dimension when G is large", {
  full <- scr_tucker3_parameter_count(
    groups = 10,
    variables = 8,
    occasions = 4,
    centroid_rank = 9,
    variable_rank = 3,
    occasion_rank = 2
  )
  reduced <- scr_tucker3_parameter_count(
    groups = 10,
    variables = 8,
    occasions = 4,
    centroid_rank = 2,
    variable_rank = 3,
    occasion_rank = 2
  )

  expect_lt(reduced, full)

  details <- scr_tucker3_parameter_count(
    groups = 10,
    variables = 8,
    occasions = 4,
    centroid_rank = 2,
    variable_rank = 3,
    occasion_rank = 2,
    details = TRUE
  )

  expect_equal(details[["centroid_subspace"]], 14L)
  expect_equal(details[["total"]], reduced)
})

test_that("Tucker3 mean structure rejects an uncentered centroid basis", {
  expect_error(
    scr_tucker3_means(
      grand_mean = matrix(0, 2, 2),
      centroid_basis = matrix(c(1, 0, 0), ncol = 1),
      variable_basis = matrix(c(1, 0), ncol = 1),
      occasion_basis = matrix(c(1, 0), ncol = 1),
      core = array(1, dim = c(1, 1, 1)),
      probabilities = c(0.2, 0.3, 0.5)
    ),
    "centering constraint"
  )
})
