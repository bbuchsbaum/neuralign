# Backend correctness reproductions

Run these scripts from the repository root with the optional backend packages
installed:

```sh
Rscript tools/repro/repro-kema-semantics.R
Rscript tools/repro/repro-cone-accuracy.R
```

`repro-kema-semantics.R` verifies that neuralign exposes manifoldalign's KEMA
training-score blocks directly, in `components x observations` orientation.

`repro-cone-accuracy.R` distinguishes adapter orientation from estimator
quality. The frozen fixture is an exact row permutation, so its oracle
reconstruction error is zero. CONE remains disabled until its assignment meets
the corresponding mandatory regression test.
