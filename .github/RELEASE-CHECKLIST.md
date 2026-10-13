# scr3way 0.1.0 release checklist

The package source is frozen at the following CRAN source-package fingerprint:

- file: `scr3way_0.1.0.tar.gz`
- size: **96,300 bytes**
- SHA-256:
  `c95faddc36642397f7c6c4666930bfc6b5290ab2f73532accdfc122933e41cff`

Do not change files included in the source package after this point without
regenerating the fingerprint and rerunning all release checks.

## Manual release sequence

1. Run **External-CRAN-checks** with target `win-release`.
2. Run **External-CRAN-checks** with target `mac-release`.
3. Review the external builder results. Resolve every ERROR/WARNING and any
   unexpected NOTE before proceeding.
4. Create tag `v0.1.0` from the exact reviewed `main` commit.
5. Wait for **Draft-release** to finish.
6. Inspect the draft GitHub release and attached files.
7. Verify the attached source tarball SHA-256 is exactly:
   `c95faddc36642397f7c6c4666930bfc6b5290ab2f73532accdfc122933e41cff`.
8. Keep the GitHub release as draft until external checks are accepted.
9. Submit the exact checked source tarball to CRAN manually.
10. After CRAN acceptance, publish the GitHub release and move development
    version to `0.1.0.9000` or the next planned development version.

## Completed gates

- [x] Ubuntu / R release
- [x] Ubuntu / R-devel
- [x] Windows / R release
- [x] macOS / R release
- [x] committed roxygen documentation
- [x] vignette build
- [x] runnable public API examples
- [x] source tarball audit
- [x] built-tarball `R CMD check --as-cran`
- [x] URL and DOI probes
- [x] formal `inst/CITATION`
- [ ] win-builder
- [ ] external macOS builder
- [ ] tag `v0.1.0`
- [ ] draft release reviewed
- [ ] CRAN submission
