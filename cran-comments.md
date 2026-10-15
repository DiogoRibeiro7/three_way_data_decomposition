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

## Final source-package fingerprint

The final source package attached to the published GitHub release
`V0.1.0` completed the release check successfully.

An earlier preflight build had a different archive SHA-256 because the source
tarball was rebuilt, but the complete contents manifest is unchanged. The
release asset below is the submission artifact.

Built source package:

- file: `scr3way_0.1.0.tar.gz`
- size: 96,300 bytes
- SHA-256:
  `c7f75ead6b6717036d485f2ed14a68bef6066c9da535196a6b71215a1563645a`

The preflight verified:

- `R CMD build` completes successfully;
- the built tarball passes `R CMD check --as-cran`;
- the only NOTE is `New submission`;
- package homepage and BugReports URLs are reachable;
- the cited DOI is reachable;
- development-only repository material is excluded from the source tarball;
- required package code, manuals, tests, vignette, deterministic regression
  inputs, and package citation metadata are present.

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
- [ ] tag `v0.1.0`
- [ ] inspect draft GitHub release
- [ ] final CRAN submission
