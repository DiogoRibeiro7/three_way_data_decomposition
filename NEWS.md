# scr3way 0.1.0

First CRAN release candidate.

## Package foundation

- Renamed the package from `threeway` to `scr3way` in preparation for CRAN.
- Consolidated the SCR baseline, Tucker3 extension, structural model selection,
  stability diagnostics, covariance extensions, sparse projections, and robust
  Student-t likelihood under one package identity.
- Reworked the README around the current scientific and package capabilities.
- Added CRAN-oriented package metadata, build exclusions, CI, and release
  preparation files.

## Research implementation

- S3, S2, and homoscedastic Gaussian-mixture baselines.
- Tucker3 centroid-mode reduction and direct S3 comparison.
- BIC, ICL, group-count, rank, and stability selection.
- Kronecker separability diagnostics and nugget covariance profiling/refinement.
- Sparse variable/occasion loading projection and active-set selection.
- Student-t likelihood and latent robustness weights.
