test_that("deterministic MATLAB regression inputs are installed", {
  root <- scr_matlab_fixture_root()
  input_dir <- file.path(root, "inputs")

  expect_true(nzchar(root))
  expect_true(dir.exists(input_dir))

  X <- read_scr_numeric_fixture(file.path(input_dir, "X.csv"))
  U <- read_scr_numeric_fixture(file.path(input_dir, "U.csv"))
  TB2 <- read_scr_numeric_fixture(file.path(input_dir, "TB_s2.csv"))
  TB3 <- read_scr_numeric_fixture(file.path(input_dir, "TB_s3.csv"))
  TC3 <- read_scr_numeric_fixture(file.path(input_dir, "TC_s3.csv"))

  expect_equal(dim(X), c(12L, 4L))
  expect_equal(dim(U), c(12L, 2L))
  expect_equal(rowSums(U), rep(1, 12), tolerance = 1e-12)
  expect_equal(dim(TB2), c(4L, 2L))
  expect_equal(dim(TB3), c(2L, 1L))
  expect_equal(dim(TC3), c(2L, 1L))
})

test_that("H matches generated MATLAB reference outputs", {
  skip_if_not(
    scr_matlab_reference_outputs_available(),
    "MATLAB reference outputs have not been generated."
  )

  root <- scr_matlab_fixture_root()
  input_dir <- file.path(root, "inputs")
  output_dir <- file.path(root, "outputs")

  X <- read_scr_numeric_fixture(file.path(input_dir, "X.csv"))
  U0 <- read_scr_numeric_fixture(file.path(input_dir, "U.csv"))

  fit <- fit_homoscedastic_gaussian_mixture(
    X,
    U0,
    tolerance = 1e-8,
    max_iter = 10000L
  )

  expect_equal(
    fit$U,
    read_scr_numeric_fixture(file.path(output_dir, "H_U.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    fit$Mmu,
    read_scr_numeric_fixture(file.path(output_dir, "H_Mmu.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    fit$Sig,
    read_scr_numeric_fixture(file.path(output_dir, "H_Sig.csv")),
    tolerance = 1e-7
  )

  scalars <- as.numeric(
    read_scr_numeric_fixture(file.path(output_dir, "H_scalars.csv"))
  )
  expect_equal(fit$dif, scalars[1], tolerance = 1e-7)
  expect_equal(fit$like, scalars[2], tolerance = 1e-7)
  expect_equal(fit$bic, scalars[3], tolerance = 1e-7)
})

test_that("S2 matches invariant MATLAB reference outputs", {
  skip_if_not(
    scr_matlab_reference_outputs_available(),
    "MATLAB reference outputs have not been generated."
  )

  root <- scr_matlab_fixture_root()
  input_dir <- file.path(root, "inputs")
  output_dir <- file.path(root, "outputs")

  X <- read_scr_numeric_fixture(file.path(input_dir, "X.csv"))
  U0 <- read_scr_numeric_fixture(file.path(input_dir, "U.csv"))
  TB <- read_scr_numeric_fixture(file.path(input_dir, "TB_s2.csv"))

  fit <- fit_scr_s2(
    X,
    U0,
    TB,
    tolerance = 1e-8,
    max_iter = 10000L
  )

  expect_equal(
    fit$U,
    read_scr_numeric_fixture(file.path(output_dir, "S2_U.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    fit$M,
    read_scr_numeric_fixture(file.path(output_dir, "S2_M.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    fit$SV,
    read_scr_numeric_fixture(file.path(output_dir, "S2_SV.csv")),
    tolerance = 1e-7
  )

  scalars <- as.numeric(
    read_scr_numeric_fixture(file.path(output_dir, "S2_scalars.csv"))
  )
  expect_equal(fit$like, scalars[1], tolerance = 1e-7)
  expect_equal(fit$bic, scalars[2], tolerance = 1e-7)
})

test_that("S3 matches invariant MATLAB reference outputs", {
  skip_if_not(
    scr_matlab_reference_outputs_available(),
    "MATLAB reference outputs have not been generated."
  )

  root <- scr_matlab_fixture_root()
  input_dir <- file.path(root, "inputs")
  output_dir <- file.path(root, "outputs")

  X <- read_scr_numeric_fixture(file.path(input_dir, "X.csv"))
  U0 <- read_scr_numeric_fixture(file.path(input_dir, "U.csv"))
  TB <- read_scr_numeric_fixture(file.path(input_dir, "TB_s3.csv"))
  TC <- read_scr_numeric_fixture(file.path(input_dir, "TC_s3.csv"))
  SV <- read_scr_numeric_fixture(file.path(input_dir, "SV_s3.csv"))
  SO <- read_scr_numeric_fixture(file.path(input_dir, "SO_s3.csv"))

  fit <- fit_scr_s3(
    X,
    U0,
    TB,
    TC,
    SV,
    SO,
    tolerance = 1e-8,
    max_iter = 10000L
  )

  expect_equal(
    fit$U,
    read_scr_numeric_fixture(file.path(output_dir, "S3_U.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    fit$M,
    read_scr_numeric_fixture(file.path(output_dir, "S3_M.csv")),
    tolerance = 1e-7
  )
  expect_equal(
    kronecker(fit$SO, fit$SV),
    read_scr_numeric_fixture(file.path(output_dir, "S3_Sigma.csv")),
    tolerance = 1e-7
  )

  scalars <- as.numeric(
    read_scr_numeric_fixture(file.path(output_dir, "S3_scalars.csv"))
  )
  expect_equal(fit$like, scalars[1], tolerance = 1e-7)
  expect_equal(fit$bic, scalars[2], tolerance = 1e-7)
})
