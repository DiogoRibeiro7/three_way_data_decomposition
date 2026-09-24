scr_matlab_fixture_root <- function() {
  system.file(
    "extdata",
    "matlab_reference",
    package = "threeway"
  )
}


read_scr_numeric_fixture <- function(path) {
  as.matrix(
    utils::read.csv(
      path,
      header = FALSE,
      check.names = FALSE
    )
  )
}


scr_matlab_reference_outputs_available <- function() {
  root <- scr_matlab_fixture_root()
  output_dir <- file.path(root, "outputs")

  required <- c(
    "H_U.csv",
    "H_Mmu.csv",
    "H_Sig.csv",
    "H_scalars.csv",
    "S2_U.csv",
    "S2_SV.csv",
    "S2_M.csv",
    "S2_scalars.csv",
    "S3_U.csv",
    "S3_M.csv",
    "S3_Sigma.csv",
    "S3_scalars.csv",
    "REFERENCE_SOURCE.txt"
  )

  nzchar(root) &&
    dir.exists(output_dir) &&
    all(file.exists(file.path(output_dir, required)))
}
