# CRAN submission: scr3way 0.1.0

## Submission

This is the first submission of `scr3way` to CRAN.

The package implements simultaneous clustering and dimensionality reduction
for three-way data and cites the underlying SCR methodology in DESCRIPTION.

## Test environments

The release candidate has been checked on:

- Ubuntu latest, R release
- Ubuntu latest, R-devel
- Windows latest, R release
- macOS latest, R release

All four environments completed with:

- 0 ERRORs
- 0 WARNINGs
- 1 NOTE: New submission

## R CMD check

The only NOTE reported by the release-candidate matrix is the expected
first-submission note:

`New submission`

No additional check notes are currently present.

## Additional notes

- This is a new submission.
- The package contains no compiled code.
- Examples and vignettes do not require internet access.
- Generated `man/` documentation and `NAMESPACE` are committed and checked
  for roxygen consistency in CI.
- Development-only validation material under `rossi/` and `inst/matlab/`
  is excluded from the CRAN source package.
- The package distributes the maintained R implementation. External MATLAB
  material is used only as scientific/development reference material and is
  not distributed in the package tarball.
- The package cites Rocci, Vichi and Ranalli (2025),
  doi:10.1007/s00180-024-01478-1.

## Final pre-submission checks

- [x] Four-platform release-candidate CI passes
- [ ] Source tarball preflight passes
- [ ] `R CMD build .`
- [ ] `R CMD check --as-cran scr3way_0.1.0.tar.gz`
- [ ] win-builder / Windows CRAN-style check
- [ ] URLs and DOI checked
- [ ] package size reviewed
- [ ] final tarball contents reviewed
