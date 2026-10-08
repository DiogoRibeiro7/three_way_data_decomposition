scr_matlab_fixture_root <- function() {
  system.file(
    "extdata",
    "matlab_reference",
    package = "scr3way"
  )
}


scr_matlab_reference_outputs_available <- function() {
  root <- scr_matlab_fixture_root()

  if (!nzchar(root)) {
    return(FALSE)
  }

  all(
    file.exists(
      file.path(
        root,
        "outputs",
        c("mixhom.csv", "t2mixt.csv", "t3mixs.csv")
      )
    )
  )
}
