# scr3way

**scr3way** is an R package for simultaneous clustering and dimensionality reduction of three-way data.

The package is built around the Simultaneous Clustering and Reduction (SCR) methodology of Roberto Rocci, Maurizio Vichi, and Monia Ranalli:

> Rocci, R., Vichi, M. & Ranalli, M. (2025). *Mixture models for simultaneous classification and reduction of three-way data*. Computational Statistics, 40, 469–507. <https://doi.org/10.1007/s00180-024-01478-1>

The repository began as an R reconstruction and validation effort against the authors' public MATLAB reference implementation. It now provides a tested SCR baseline together with research extensions for Tucker3 structure, model selection, covariance departures, sparsity, stability, and robust heavy-tailed modelling.

## Models

The baseline comparison contains three models:

- **S3** — three-way SCR with a Tucker2 mean structure and separable common covariance;
- **S2** — two-way SCR applied to vectorized three-way observations;
- **H** — homoscedastic Gaussian mixture with unrestricted component means.

The package also contains a **Tucker3 SCR extension** that reduces the group/centroid mode in addition to the variable and occasion modes.

## Main capabilities

### SCR estimation and validation

- `fit_scr_s3()` — three-way S3 baseline;
- `fit_scr_s2()` — vectorized S2 baseline;
- `fit_homoscedastic_gaussian_mixture()` — H comparator;
- `run_scr_simulation()` and `run_scr_experiment_grid()` — reproducible simulation infrastructure;
- hard partitions and adjusted Rand index utilities;
- MATLAB-oriented numerical regression infrastructure retained for development validation.

### Tucker3 extension

- `fit_scr_s3_tucker3()`;
- `scr_tucker3_mean_update()`;
- explicit parameter counting and rank constraints;
- direct S3-versus-Tucker3 benchmarking.

### Structural model selection

- joint rank selection over `(P, Q, R)`;
- outer selection over the number of groups `G`;
- BIC and ICL criteria;
- deterministic tie-breaking;
- label-invariant subsampling stability selection.

### Covariance extensions

- nearest Kronecker covariance approximation and separability diagnostics;
- Kronecker-plus-nugget covariance
  `Sigma_O %x% Sigma_V + tau I`;
- conditional nugget profile likelihood;
- conservative alternating nugget refinement.

### Sparse discriminating subspaces

- row-group soft thresholding of variable and occasion loadings;
- covariance-metric re-orthonormalization;
- sparse Tucker3 mean projection;
- deterministic support paths;
- conditional active-set BIC/ICL support selection.

### Robust modelling

- multivariate Student-t Tucker3 likelihood;
- posterior memberships and Mahalanobis distances;
- latent robustness weights for downweighting remote observations;
- Student-t parameter-count support.

## Installation

The package is currently in development and is being prepared for CRAN.

```r
# install.packages("remotes")
remotes::install_github("DiogoRibeiro7/three_way_data_decomposition")
```

Then load it with:

```r
library(scr3way)
```

## Minimal example

```r
library(scr3way)

set.seed(1)

X <- matrix(rnorm(80 * 4), nrow = 80)
membership <- matrix(runif(80 * 2), nrow = 80)
membership <- membership / rowSums(membership)

fit <- fit_scr_s3_tucker3(
  X = X,
  membership = membership,
  centroid_rank = 1,
  variable_rank = 1,
  occasion_rank = 1,
  variable_covariance = diag(2),
  occasion_covariance = diag(2),
  max_iter = 20
)

fit$ranks
fit$bic
```

## Repository structure

- `R/` — package implementation;
- `tests/testthat/` — deterministic unit and numerical regression tests;
- `docs/scr-research-program.md` — detailed scientific roadmap and methodology notes;
- `rossi/` — legacy/reference R material retained in the repository but excluded from the CRAN source package;
- `inst/matlab/` and MATLAB reference fixtures — development-only validation tooling, excluded from the CRAN source package.

The public MATLAB code is used only as an external scientific reference. It is not distributed as part of the `scr3way` package tarball.

## CRAN status

The repository is being converted into a CRAN-ready package under the name `scr3way`.

Before the first `0.1.0` submission, the remaining release work includes:

- generated and audited `man/` documentation;
- examples for exported user-facing functions;
- introductory and methodology vignettes;
- clean `R CMD check --as-cran` results on current R release and R-devel;
- Windows and multi-platform checks;
- final `cran-comments.md`.

## Documentation

The detailed reconstruction history, model definitions, validation rules, and research extensions are documented in [docs/scr-research-program.md](docs/scr-research-program.md).

## License

The maintained R package code is released under the MIT License. See [LICENSE](LICENSE).
