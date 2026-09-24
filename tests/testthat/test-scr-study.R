test_that("simulation configurations keep MATLAB and paper scenario I distinct", {
  matlab <- scr_simulation_config("I", "matlab")
  paper <- scr_simulation_config("I", "paper")
  scenario_two <- scr_simulation_config("II", "paper")

  expect_equal(
    unname(unlist(matlab[c("J", "Q", "K", "R")])),
    c(5, 2, 4, 2)
  )
  expect_equal(
    unname(unlist(paper[c("J", "Q", "K", "R")])),
    c(5, 2, 5, 2)
  )
  expect_equal(
    unname(unlist(scenario_two[c("J", "Q", "K", "R")])),
    c(20, 5, 5, 2)
  )
})

test_that("simulation preprocessing is explicit and dimension preserving", {
  result <- run_scr_simulation(
    n = 40,
    groups = 2,
    n_starts = 1,
    dgp = 1,
    n_simulations = 1,
    scenario = "I",
    source = "matlab",
    max_iter = 5,
    preprocess = function(X) X,
    keep_data = FALSE
  )

  expect_equal(dim(result$ari), c(1L, 3L))
  expect_named(as.list(result$ari[1, ]), c("S3", "S2", "H"))
  expect_equal(result$preprocessing, "user-supplied")
  expect_equal(result$initialization, "random-soft")
  expect_null(result$X)
  expect_null(result$U_true)
})

test_that("simula1 no longer returns placeholder model results", {
  result <- simula1(
    N = 40,
    G = 2,
    nrep = 1,
    dgp = 1,
    ns = 1,
    source = "matlab",
    max_iter = 5
  )

  expect_equal(dim(result$ari), c(1L, 3L))
  expect_equal(dim(result$Xt), c(40L, 20L, 1L))
  expect_equal(dim(result$Utruet), c(40L, 2L, 1L))
  expect_true(all(is.finite(result$ari)))
  expect_named(
    result$details$diagnostics[[1]],
    c("S3", "S2", "H")
  )
})

test_that("custom membership initializers are normalized before fitting", {
  initializer <- function(X, groups, start, seed) {
    matrix(
      rep(c(2, 1), length.out = nrow(X) * groups),
      nrow = nrow(X),
      ncol = groups,
      byrow = TRUE
    )
  }

  result <- run_scr_simulation(
    n = 40,
    groups = 2,
    n_starts = 1,
    dgp = 1,
    n_simulations = 1,
    scenario = "I",
    source = "matlab",
    max_iter = 5,
    membership_initializer = initializer,
    keep_data = FALSE
  )

  expect_equal(result$initialization, "user-supplied")
  expect_true(all(is.finite(result$ari)))
})

test_that("simulation controls reject invalid configurations", {
  expect_error(
    run_scr_simulation(
      n = 20,
      groups = 1,
      n_starts = 1,
      dgp = 1,
      n_simulations = 1
    ),
    "groups must be at least two"
  )

  expect_error(
    run_scr_simulation(
      n = 20,
      groups = 2,
      n_starts = 1,
      dgp = 5,
      n_simulations = 1
    ),
    "dgp must be one of"
  )

  expect_error(
    run_scr_simulation(
      n = 20,
      groups = 2,
      n_starts = 1,
      dgp = 1,
      n_simulations = 1,
      preprocess = 1
    ),
    "NULL or a function"
  )
})
