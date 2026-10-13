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

## Source-package preflight

The CRAN source-package preflight completed successfully.

Built source package:

- file: `scr3way_0.1.0.tar.gz`
- size: 95,964 bytes
- SHA-256:
  `ac16c9af137edf8d16a802041b4d432efa3cc77cc8a6f387024dd5fbeabd0ba5`

The preflight also verified:

- `R CMD build` completes successfully;
- the built tarball passes `R CMD check --as-cran`;
- the only NOTE is `New submission`;
- package homepage and BugReports URLs are reachable;
- the cited DOI is reachable;
- development-only repository material is excluded from the source tarball;
- required package code, manuals, tests, vignette, and deterministic regression
  inputs are present.

## R CMD check

The only NOTE currently reported is the expected first-submission note:

`New submission`

No ERRORs or WARNINGs are present.

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
- [x] Source tarball preflight passes
- [x] `R CMD build .`
- [x] `R CMD check --as-cran scr3way_0.1.0.tar.gz`
- [ ] win-builder / Windows CRAN-style check
- [x] URLs and DOI checked
- [x] package size reviewed
- [x] final tarball contents reviewed
- [ ] external macOS builder check
- [ ] final CRAN submission
