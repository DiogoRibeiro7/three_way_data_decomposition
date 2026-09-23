# Three-Way Data Decomposition

R implementation, validation framework, and research extensions for simultaneous clustering and dimensionality reduction of three-way data.

The statistical baseline is the Simultaneous Clustering and Reduction (SCR) model developed by Roberto Rocci, Maurizio Vichi, and Monia Ranalli in *Mixture models for simultaneous classification and reduction of three-way data*, Computational Statistics 40, 469–507 (2025), DOI: 10.1007/s00180-024-01478-1.

This repository is not intended to be a generic collection of tensor algorithms. Tensor decompositions are supporting numerical machinery. The primary objective is to reproduce the SCR methodology faithfully in R, validate it against the authors' MATLAB reference implementation, and then develop statistically meaningful extensions.

## Statistical baseline

The baseline comparison contains three models:

- **S3**: three-way SCR, with Tucker2 structure on the group means and a Kronecker-structured common covariance;
- **S2**: two-way SCR applied to the vectorised three-way observations;
- **H**: ordinary homoscedastic Gaussian mixture model.

The legacy R files under `rossi/` are a port of the authors' public MATLAB code in `moniar412/SCR3waydata`. They are currently treated as reference material, not yet as validated package code.

## Repository layers

- `R/`: tested package API and shared numerical infrastructure;
- `tests/testthat/`: deterministic tests and numerical regression tests;
- `rossi/`: legacy R translation of the original SCR reference implementation;
- simulation and reproducibility code: to be rebuilt from the legacy scripts after the baseline port has been validated.

The detailed model-to-code map, known porting discrepancies, and research roadmap are documented in [docs/scr-research-program.md](docs/scr-research-program.md).

## Research programme

Development proceeds in two stages.

First, the repository will establish a faithful R reproduction of S3, S2, H, the two simulation scenarios, and the ARI experiments reported in the paper. Numerical equivalence with the MATLAB reference implementation is the acceptance criterion.

Second, the SCR formulation will be used as a baseline for extensions. Candidate directions include centroid-mode reduction through Tucker3, joint model selection over clustering and reduction dimensions, controlled departures from exact Kronecker covariance structure, and robust or sparse variants. Extensions will be implemented only after the baseline is reproducible.

## Package status

The package infrastructure is being modernised incrementally. Existing wrappers around `rTensor` provide supporting tensor decompositions, but these are not the scientific contribution of the project.

The SCR algorithms themselves are still under reconstruction and should not yet be treated as validated production implementations.

## Reference

Rocci, R., Vichi, M. & Ranalli, M. (2025). *Mixture models for simultaneous classification and reduction of three-way data*. Computational Statistics, 40, 469–507. https://doi.org/10.1007/s00180-024-01478-1

## License

See [LICENSE](LICENSE).
