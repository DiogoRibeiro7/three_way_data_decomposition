# SCR research programme

## Purpose

This repository is an R reconstruction and extension of the Simultaneous Clustering and Reduction (SCR) methodology for three-way data developed by Roberto Rocci, Maurizio Vichi, and Monia Ranalli.

The reference article is:

> Rocci, R., Vichi, M. & Ranalli, M. (2025). *Mixture models for simultaneous classification and reduction of three-way data*. Computational Statistics, 40, 469–507. DOI: 10.1007/s00180-024-01478-1.

The authors also provide a public MATLAB reference implementation at:

`moniar412/SCR3waydata`

The legacy files currently stored under `rossi/` correspond almost one-to-one with the filenames in that MATLAB repository. They are therefore scientifically important reference material. They must not be discarded or treated as unrelated exploratory scripts.

## Statistical model

Let each observation contain (J) variables measured on (K) occasions and let (G) denote the number of mixture components.

The three-way SCR model assumes a homoscedastic Gaussian mixture with a common covariance matrix constrained by

[
\Sigma = \Sigma_O \otimes \Sigma_V,
]

where (\Sigma_V) is the variable covariance and (\Sigma_O) is the occasion covariance.

The group means are constrained through a Tucker2 representation,

[
\mu_g = \mu + (C \otimes B)\eta_g,
]

where (B \in \mathbb{R}^{J \times Q}), (C \in \mathbb{R}^{K \times R}), and (\eta_g \in \mathbb{R}^{QR}). Thus the group means are represented in a (QR)-dimensional discriminating subspace rather than the original (JK)-dimensional space.

The baseline comparison used in the paper is:

| Model | Legacy function | Mean structure | Covariance structure |
| --- | --- | --- | --- |
| S3 | `t3mixs()` | Tucker2, (C \otimes B) | (\Sigma_O \otimes \Sigma_V) |
| S2 | `t2mixt()` | reduced (QR)-dimensional subspace after vectorisation | unrestricted common covariance on the vectorised observations |
| H | `mixhom()` | unrestricted component means | unrestricted common covariance |

The model is fitted by an EM-like block-coordinate algorithm. Posterior memberships, mixing proportions, covariance factors, loadings, and latent group scores are updated iteratively.

## Legacy-code map

### Core estimators

| R file | MATLAB reference | Statistical role | Target disposition |
| --- | --- | --- | --- |
| `rossi/t3mixs.R` | `t3mixs.m` | S3 three-way SCR estimator | rebuild as validated package core |
| `rossi/t2mixt.R` | `t2mixt.m` | S2 two-way SCR comparator | rebuild as validated package core |
| `rossi/mixhorm.R` | `mixhom.m` | H homoscedastic Gaussian mixture comparator | rebuild as validated package core |
| `rossi/genmixhet.R` | `genmixhet.m` | Gaussian-mixture data generator | rebuild as simulation infrastructure |
| `rossi/ftoh.R` | `ftoh.m` | convert fuzzy memberships to hard partition | validate, then promote to package utility |
| `rossi/mrand.R` | `mrand.m` | adjusted/modified Rand index used by simulations | validate against a standard ARI implementation |

### Simulation design

The original MATLAB repository separates the simulation into two dimensional regimes and then runs the full factorial experiment for (G \in \{3,5,7\}).

| File family | Role |
| --- | --- |
| `simula1.*` | scenario I: few variables / lower-dimensional setting |
| `simula2.*` | scenario II: many variables / higher-noise setting |
| `simulaG3ari*.R` | experiment grid for (G=3) |
| `simulaG5ari*.R` | experiment grid for (G=5) |
| `simulaG7ari*.R` | experiment grid for (G=7) |

The paper compares S3, S2, and H over four data-generating processes:

1. structured means and structured covariance;
2. unstructured means and structured covariance;
3. structured means and unstructured covariance;
4. unstructured means and unstructured covariance.

The experimental factors additionally include sample size, number of groups, and number of random starts. Each setting in the published study uses 250 generated samples.

## Reproduction blockers already identified

The current R files are not yet a faithful executable reproduction of the MATLAB code. The following issues must be resolved before numerical results from the R implementation can be trusted.

### 1. S3 and S2 are placeholders in the R simulation drivers

The current `simula1.R` and `simula2.R` do not call the translated S3 and S2 estimators in the model-comparison loop. Random cluster assignments and dummy likelihood values are used instead.

The MATLAB reference code calls `t3mixs()`, `t2mixt()`, and `mixhom()` directly.

This is the highest-priority reproduction defect.

### 2. Quadratic-form translation in S2 and S3 needs correction

The MATLAB code squares the transformed residuals after matrix multiplication. The current R translations of `t2mixt()` and `t3mixs()` contain expressions where operator precedence can instead square only the final matrix factor.

The R implementation must explicitly compute the transformed residual matrix first and then square it elementwise.

### 3. ARI contingency construction must be checked

The MATLAB simulation constructs the contingency matrix from transposed hard-membership matrices. The current R simulation code does not consistently preserve that orientation.

ARI tests must compare the legacy implementation against a trusted reference implementation.

### 4. Scenario-I dimensional discrepancy

The public MATLAB `simula1.m` currently sets (J=5), (Q=2), (K=4), and (R=2).

The 2025 article reports scenario I as (J=5), (Q=2), (K=5), and (R=2).

This discrepancy must be resolved explicitly. Until then, the repository must distinguish:

- reproduction of the public MATLAB code; and
- reproduction of the published 2025 simulation tables.

They are not automatically the same target.

### 5. Missing or implicit MATLAB dependencies

The simulation code uses preprocessing and initialisation routines that are not self-contained in the public reference repository. The R implementation must make every such dependency explicit and tested rather than silently substitute a different operation.

### 6. Simulation scripts execute work on source

Several legacy `simulaG*ari*.R` files call their simulation function at file load time. Reproducibility scripts must be pure definitions with execution controlled by an explicit runner.

### 7. Repeated helper names

Each experiment driver defines its own `determine_r()`. The rebuilt simulation layer should represent the experiment grid as data rather than encode it through repeated procedural helpers.

## Acceptance criteria for the baseline

The baseline is considered reproduced only when all of the following hold:

1. S3, S2, and H have deterministic unit tests on small synthetic examples.
2. Individual parameter-update steps are tested where practical.
3. R and MATLAB produce matching or numerically equivalent likelihood trajectories on fixed inputs and fixed initial values.
4. Hard partitions and ARI values agree under deterministic seeds.
5. Both simulation scenarios can be run without placeholders.
6. The (G=3,5,7) experiment grid is generated from a single configuration object.
7. Published simulation settings are represented explicitly and separately from the historical MATLAB-code settings where they differ.
8. Reproduction outputs store the seed, configuration, package version, and model convergence diagnostics.

## Package architecture

The intended package architecture is model-centred rather than file-centred.

```text
R/
  scr_fit.R
  scr_s3.R
  scr_s2.R
  mixture_h.R
  memberships.R
  simulation.R
  metrics.R
  model_selection.R

tests/testthat/
  test-s3.R
  test-s2.R
  test-h.R
  test-memberships.R
  test-simulation.R
  test-matlab-regression.R

inst/
  extdata/
    matlab_reference/
  simulations/
    scenario_i.R
    scenario_ii.R

rossi/
  # retained temporarily as immutable port/reference material
```

The legacy `rossi/` directory should remain intact until the package implementation has numerical regression coverage against it and, where possible, against the MATLAB source.

## Research extensions

The faithful baseline is a prerequisite, not the final research goal.

### Tucker3 centroid reduction

The published paper notes the possibility of replacing the Tucker2 mean structure with a Tucker3 structure that also reduces the centroid/group mode:

[
\mu_{gjk}
=
\mu_{jk}
+
\sum_{p=1}^{P}
\sum_{q=1}^{Q}
\sum_{r=1}^{R}
a_{gp} b_{jq} c_{kr}\eta_{pqr}.
]

This becomes especially relevant as (G) grows. The extension requires identifiability constraints, an estimation algorithm, parameter counting, model selection, and comparison against S3.

### Joint structural model selection

A useful general interface should allow comparison over

[
(G,Q,R)
]

and, for Tucker3 extensions,

[
(G,P,Q,R).
]

BIC is already present in the legacy code, but model selection should be made explicit and tested. ICL and stability criteria are natural additional candidates.

### Controlled departures from exact Kronecker covariance

S3 assumes exact covariance separability. A broader family can introduce a parameterised or penalised departure from

[
\Sigma_O \otimes \Sigma_V
]

while preserving the computational advantages of the three-way structure.

### Sparse discriminating subspaces

Penalties or structured sparsity on (B) and (C) can move the method from latent dimensionality reduction toward explicit selection of discriminating variables and occasions.

### Robust SCR

Heavy-tailed or contaminated component distributions can be studied after the Gaussian baseline is stable. Robustness should be integrated with the simultaneous reduction problem rather than implemented as an unrelated mixture variant.

## Development rule

No research extension should be evaluated against the current placeholder simulations.

The order is:

1. reconstruct the published baseline;
2. validate against the MATLAB implementation;
3. reproduce the reported experimental design;
4. establish a unified model API;
5. introduce extensions one at a time;
6. compare extensions against the validated S3/S2/H baseline.

This ordering keeps software modernisation, scientific reproduction, and methodological novelty logically separate.
