# CRAN submission: scr3way

This file tracks preparation for the first CRAN release. It is intentionally
excluded from the source package.

## Release target

- Package: scr3way
- Version: 0.1.0
- Submission type: first submission

## Required checks before submission

- [ ] `R CMD build .`
- [ ] `R CMD check --as-cran scr3way_0.1.0.tar.gz`
- [ ] R release: 0 ERRORs, 0 WARNINGs, justified NOTES only
- [ ] R-devel: 0 ERRORs, 0 WARNINGs, justified NOTES only
- [ ] Windows / win-builder checks
- [ ] Additional platform checks where available
- [ ] URLs and DOI checked
- [ ] Examples run within CRAN time limits
- [ ] Vignettes build without internet access
- [ ] Package size and installed size reviewed
- [ ] Copyright/licensing audit complete

## Current development notes

The repository contains development-only validation material under `rossi/`,
`inst/matlab/`, and `inst/extdata/matlab_reference/`. These paths are
excluded from the CRAN source package. The submitted package will distribute
the maintained R implementation and cite the original SCR methodology.

The final release PR will update the version to 0.1.0 and replace this checklist
with the actual check environments and results.
