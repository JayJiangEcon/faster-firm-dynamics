# Build and run guide

Run the commands from the repository root. The package includes complete
Fortran source and the six fixed input files; no separate model-source package
is needed for these examples.

## Requirements

- Python 3.9 or later.
- GNU Fortran with OpenMP support. The recorded environment used GNU Fortran
  14.2.0 with `-O3 -funroll-loops -fopenmp -m64` on Apple M3.
- NumPy for conditional measurement-noise minimization.
- SciPy for the optional outer parameter search.

Direct evaluation and compilation orchestration use the Python standard
library. Use your existing scientific Python environment for NumPy and SciPy.
The integrated source and fixed-input identities are in
[`integrated-manifest.json`](../results/integrated-manifest.json).

## Read the complete source

- [`src/fortran/RE/diag_firms_solve.f90`](../src/fortran/RE/diag_firms_solve.f90)
- [`src/fortran/DE/diag_firms_solve.f90`](../src/fortran/DE/diag_firms_solve.f90)
- [`src/fortran/base_lib.f90`](../src/fortran/base_lib.f90)

Each model has a small `evaluate.f90` command-line driver. Both complete solvers
include the supported direct and fused moment-basis paths. The build mode
selects the runner interface; it does not change the model's economic inputs.

## Compile and evaluate a fixed parameter vector

```bash
python3 build.py RE --output /tmp/rcc-re-direct
python3 run.py --build /tmp/rcc-re-direct \
  --parameters examples/re-anchor.json --output /tmp/rcc-re-evaluation
```

For DE, use `DE` and `examples/de-anchor.json`. Add `--fc /path/to/gfortran`
to select a compiler, or set `FC`. Use a new build and output directory each
time. `--prepare-only` copies and verifies source without compiling.

The runner defaults to `data/fixed/`. An optional `--inputs-dir` may point to
another copy of the same recorded files. Hashes are checked before execution.
Set `--threads` and `--timeout` explicitly when needed; defaults are eight
threads and 300 seconds. Inherited `RCC_*` overrides are cleared.

`build.json` binds source, compiler and binary identities. The evaluation
writes `stdout.txt`, `stderr.txt`, `run.json`, and `result.json`. The latter
contains the complete parameters, moments, and weighted quadratic objective.
A failed or timed-out evaluation retains its diagnostic record. Example
parameter vectors are fixed benchmark inputs.

The parameter orders are:

```text
RE: rhoz, sigmaz, etak, phi, gamma, sigmaA, etaf, etad, sigmapi
DE: theta, rhoz, sigmaz, etak, phi, gamma, sigmaA, etaf, etad, sigmapi
```

The objective uses moments 2–17 and the corresponding submatrix of the
included weight matrix. All seventeen simulated moments are retained in the
output. `RCC_FGMM` and `RCC_RAW_Q` both report the quadratic objective.

## Conditional measurement-noise minimization

```bash
python3 build.py RE --mode sigma-profile --output /tmp/rcc-re-profile
python3 run.py --build /tmp/rcc-re-profile \
  --parameters examples/re-policy.json --sigma-bounds 0 1 \
  --output /tmp/rcc-re-profile-evaluation
```

Use DE and `examples/de-policy.json` for the DE model. The interval `[0, 1]`
above illustrates the command syntax; choose bounds appropriate to your task.
The policy-vector files contain the benchmark coordinates excluding `sigmapi`.

One policy/path computation produces `m0`, `c1` and `c2`, with
`m(s) = m0 + s*c1 + s*s*c2`. `sigma_profile.py` forms the quartic objective,
finds real stationary points on the specified interval, and evaluates those
points and both endpoints using the moment residuals. It writes the selected
`sigma_pi`, reconstructed moments, complete parameter vector, and candidate
objectives to `result.json`. This is a floating-point conditional calculation;
its historical agreement with direct evaluation is documented in
[Numerical accuracy](accuracy.md).

## Outer parameter search

`search.py` uses SciPy differential evolution with a declared seed and an
explicit evaluation budget. Each candidate invokes the conditional evaluator:
eight outer coordinates for RE and nine for DE. It uses the parameter order
above with the last coordinate omitted.

Prepare a JSON specification containing:

| Field | Required value |
| --- | --- |
| `model` | `RE` or `DE` |
| `bounds` | Eight or nine `[lower, upper]` pairs, in parameter order |
| `sigma_pi_bounds` | A finite nonnegative interval `[lower, upper]` |
| `seed` | A nonnegative integer |
| `max_evaluations` | A positive integer |

Bounds must have strictly increasing endpoints. The caller supplies the
parameter domain; this interface does not infer one from benchmark points.

```bash
python3 search.py --build /tmp/rcc-re-profile --spec /path/to/search-spec.json \
  --output /tmp/rcc-re-search
```

The search is serial across candidates and uses OpenMP inside each model
call. It records every attempted evaluation, the current best complete vector,
and the actual stopping reason in `search.json`. Reaching the evaluation
budget is reported separately from optimizer convergence. Evaluation failures
stop the search and retain their diagnostics. There is no cross-run cache.
The best evaluated candidate is not a certificate of a global optimum.

This is a new runnable search interface, not a replay of the historical
search workloads described in the timing documentation. No new estimation
result is included with the package.

## Existing measurements

```bash
python3 report_results.py
```

This command reads historical observations and does not run the model.
Keep builds and numerical runs outside a cloud-synchronized checkout.
