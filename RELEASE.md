# Release evidence

This file records what must be true for a `neuralign` release. A configured
workflow is not a passing receipt, and a local macOS check is not evidence for
Linux or Windows.

## 0.2.0 gate ledger

| Gate | Required evidence | Current status |
|---|---|---|
| Full local tests | testthat result on release commit | 1,940 pass, 4 expected-condition warnings, 7 availability skips, 0 failures on macOS/R 4.5.1 |
| Local package check | exact E/W/N counts | `Status: OK` on neuralign 0.2.0: 0 errors, 0 warnings, 0 notes; unavailable `dkge` reported as INFO under the workflow's optional-Suggests policy |
| Clean source install | built tarball installs and public smoke test passes | full-vignette source tarball installed into an empty task library; smoke passed |
| Minimal Suggests | optional backends absent; package check passes | externally blocked: current run 31436052357 failed before runner assignment (zero steps) |
| Cross-platform | Linux, macOS, Windows jobs pass on release commit | externally blocked: all three hosted-runner jobs in current run 31436052357 failed before runner assignment |
| Backend conformance | KEMA/GRASP semantics and accuracy run without opt-in flags | 70 focused backend assertions pass locally |
| Ecosystem | current `dkge` provider registers and fits | externally blocked: current dispatched run 31437226597 failed before runner assignment (zero steps) |
| Documentation | roxygen output current; examples/check pass | Rd consistency, examples, installed vignettes, and vignette rebuild pass locally |
| Tracker | completed issues closed and JSONL exported | implementation issues closed; CI issue remains in progress for hosted-runner receipts |
| Git | clean release commit, `v0.2.0` tag, local/remote SHA agreement | strict 0.2.0 implementation and policy are clean and pushed; tag remains withheld pending remote gates |

## Remote CI receipt

GitHub accepted and enabled both workflow files. On 2026-08-10, [R CMD check
run 31436052357](https://github.com/bbuchsbaum/neuralign/actions/runs/31436052357)
on strict-provider commit `bc34fcc` assigned zero steps to all six jobs: Linux,
macOS, Windows, clean install, minimal Suggests, and backend conformance. The
separately dispatched [ecosystem run
31437226597](https://github.com/bbuchsbaum/neuralign/actions/runs/31437226597)
on the documentation-only policy successor `7fe6ae0` likewise assigned zero
steps to its `dkge` provider-contract job. This reproduces the prior
account-level hosted-runner or billing gate; it is not package or workflow
execution evidence. Do not tag until a fresh run starts steps and passes.

## Local commands

Run from a clean checkout with a valid UTF-8 locale:

```sh
Rscript -e 'devtools::test(stop_on_failure = TRUE)'
R CMD build .
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual neuralign_0.2.0.tar.gz
R CMD INSTALL neuralign_0.2.0.tar.gz
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
5. Create annotated tag `v0.2.0` at that exact commit and push the tag.
6. Verify local `HEAD`, `origin/main`, and the peeled tag SHA agree.

If remote CI has not yet run, the tag may trigger it, but the ledger must keep
those lanes marked pending until receipts exist. Do not describe the release as
cross-platform certified before then.
