# Release evidence

This file records what must be true for a `neuralign` release. A configured
workflow is not a passing receipt, and a local macOS check is not evidence for
Linux or Windows.

## 0.1.0 gate ledger

| Gate | Required evidence | Current status |
|---|---|---|
| Full local tests | testthat result on release commit | 1,814 pass, 4 expected-condition warnings, 7 availability skips, 0 failures on macOS/R 4.5.1 |
| Local package check | exact E/W/N counts | `Status: OK` (0 errors, 0 warnings, 0 notes); unavailable `dkge` reported as INFO |
| Clean source install | built tarball installs and public smoke test passes | full-vignette source tarball installed into an empty task library; smoke passed |
| Minimal Suggests | optional backends absent; package check passes | workflow added; remote receipt pending |
| Cross-platform | Linux, macOS, Windows jobs pass on release commit | workflow added; remote receipt pending |
| Backend conformance | KEMA/GRASP semantics and accuracy run without opt-in flags | 70 focused backend assertions pass locally |
| Ecosystem | current `dkge` provider registers and fits | scheduled workflow added; remote receipt pending |
| Documentation | roxygen output current; examples/check pass | Rd consistency, examples, installed vignettes, and vignette rebuild pass locally |
| Tracker | completed issues closed and JSONL exported | implementation issues closed; CI/release items await remote receipts |
| Git | clean release commit, `v0.1.0` tag, local/remote SHA agreement | pending |

## Local commands

Run from a clean checkout with a valid UTF-8 locale:

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
R CMD build .
R CMD check --no-manual neuralign_0.1.0.tar.gz
R CMD INSTALL neuralign_0.1.0.tar.gz
Rscript -e 'library(neuralign); x <- list(s1 = diag(4), s2 = diag(4)); stopifnot(inherits(fit_alignment(x, reference = "s1", compute_quality = FALSE), "AlignmentResult"))'
```

The backend reproductions are:

```sh
Rscript tools/repro/repro-kema-semantics.R
Rscript tools/repro/repro-cone-accuracy.R
```

## Tag protocol

1. Update this ledger with exact local counts and remote workflow links/status.
2. Confirm `DESCRIPTION` and the `NEWS.md` heading agree.
3. Export tracker state and stage only the intended release slice.
4. Commit and push `main` without force.
5. Create annotated tag `v0.1.0` at that exact commit and push the tag.
6. Verify local `HEAD`, `origin/main`, and the peeled tag SHA agree.

If remote CI has not yet run, the tag may trigger it, but the ledger must keep
those lanes marked pending until receipts exist. Do not describe the release as
cross-platform certified before then.
