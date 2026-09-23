# Three-Way Data Decomposition

Research code for multilinear data analysis, tensor decomposition, and related methods for three-way and higher-order data.

The repository currently contains implementations and experiments around tensor representations, CP and Tucker decompositions, higher-order SVD, multilinear principal component analysis, tensor support vector methods, and kernel-based clustering. The long-term goal is to turn the useful numerical routines into a small, tested R package while preserving exploratory material separately.

## Scope

The current code covers several related areas:

- tensor construction and matricisation;
- Kronecker and Khatri-Rao products;
- CANDECOMP/PARAFAC decomposition;
- Tucker decomposition and HOSVD;
- multilinear principal component analysis;
- tensor reconstruction;
- tensor-based classification experiments;
- kernel clustering experiments.

The repository is being modernised incrementally. Some files are still exploratory scripts rather than stable library code.

## Current structure

Most implementation files currently live at the repository root. This reflects the historical development of the project rather than the intended final architecture.

In particular, some files mix reusable functions with executable examples or manual checks. Those concerns will be separated before the algorithms are extended.

## Modernisation plan

The cleanup will proceed in small changes:

1. define the supported statistical and numerical scope;
2. separate reusable functions from experiments and examples;
3. introduce a conventional R package layout;
4. add deterministic unit tests for tensor algebra and decompositions;
5. add automated package checks and linting;
6. review numerical stability, convergence criteria, and input validation;
7. document the mathematical assumptions and references for each method;
8. add reproducible examples and benchmarks.

Algorithmic changes will be kept separate from repository-structure changes so that numerical behaviour can be reviewed independently.

## References

The implementations draw on standard results from multilinear algebra and tensor decomposition, including work on CP/PARAFAC, Tucker decomposition, HOSVD, and multilinear principal component analysis.

More detailed references will be attached to the corresponding methods as the package structure is introduced.

## License

See [LICENSE](LICENSE).
