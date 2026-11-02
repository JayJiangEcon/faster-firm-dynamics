# Version 0.2.0: final numerical verification

The final release workflow is **Mac input construction → remote A100 solution
and simulation → Mac numerical finishing and conditional selection**.
Expanded checks executed on November 1, 2026 cover 28 distinct policy vectors
for RE and 28 for DE. This is implementation verification, not an estimation.
[Structured results and source/build identities](../results/cuda-coverage-v020.json)
bind the final implementation.

## Final results

| Comparison | RE | DE |
|---|---:|---:|
| Distinct policy vectors | 28 | 28 |
| Binary input files compared | 336 | 336 |
| Input-file mismatches | 0 | 0 |
| Direct-output entries compared | 9,223,844 | 115,063,844 |
| Direct-output bit mismatches | 0 | 0 |
| Conditional basis/coefficient entries compared | 2,380 | 2,380 |
| Conditional basis/coefficient bit mismatches | 0 | 0 |
| Direct objectives compared | 28 | 28 |
| Direct-objective bit mismatches | 0 | 0 |
| Conditional candidates, selected sigma, moments and objective | All bit-identical | All bit-identical |

Direct-output coverage comprises complete policy, default-policy, value,
prior-value, bond-price, distribution and profit arrays; macro paths; retained
micro observations; and all 17 final moments. Conditional coverage includes
`m0`, `m+`, `m-`, `c1` and `c2`, all candidate sigmas and candidate objectives,
and the selected outputs. The three basis lanes do not count as three new
policy vectors.

## Point selection and execution identity

Every point was fixed before observing its comparison outcome. Each model uses:

- One previously verified anchor.
- Nine distinct candidates from the earlier operational search test.
- Sixteen Sobol points: eight from each independent scramble, seeds 20261102
  and 20261103.
- Two boundary checks: RE persistence endpoints or DE expectation-parameter
  endpoints. The DE solver remains DE at theta exactly zero.

The engineering box scales the packaged anchor's policy coordinates by 0.75
and 1.25, caps the persistence upper endpoint at 0.98, and spans [0, 1.5] for
DE theta. Sobol direct-evaluation sigma values span [0, 1.5]; conditional
selection uses the same interval. This is a declared verification box, not a
recommended estimation domain or a global-coverage result. Full grids, fixed
draws, stopping controls, target definitions and moments 2–17 are retained.
No point was dropped or replaced after seeing its result.

CPU references came from the packaged **accelerated** Fortran implementation,
executed independently once in direct mode and once in profile mode per point:
112 processes in total. The original slow baseline was not rerun. The observer
writes arrays after computation; at the original anchors, 26 recorded files
per model exactly matched the earlier independent references.

The same immutable CPU outputs were reused for the initial and final GPU
comparisons. The final GPU implementation was executed at all 56 points in
both modes. Comparators import neither numerical implementation. Completed
CPU-to-predecessor array comparisons also support exact predecessor-to-final
byte comparisons for part of the DE record; final moments and conditional
outputs are compared to their CPU references. Thus saved outputs are reused
without treating a previous pass as a substitute for checking changed outputs.

The environment uses GNU Fortran 14.2 on the M3 and CUDA 12.4 on an
A100-SXM4-80GB server. The local finishing library and remote CUDA binary are
bound separately. Bitwise agreement applies to these finite points and
execution identities, not every possible parameter or architecture.

## What the wider checks found and fixed

The initial extended comparison found a one-ULP difference in moment 17 at
three RE points and two DE points; the maximum absolute difference was
5.421010862427522e-20. All other compared solver and simulation arrays were
exact. A probe using 252 identical macro-output inputs isolated one differing
Mac/server log result. Separately, Python/BLAS and Fortran reduction order
changed the direct objective at 18 RE and 15 DE points, by at most
1.1102230246251565e-16 and 1.3322676295501878e-15, respectively.

These differences were corrected, not accepted under a tolerance. The CUDA
worker now returns the already computed macro output path. A small local
Fortran library computes its growth variance and the direct weighted objective
with the Mac arithmetic/library path. The same macro moment enters all three
conditional lanes before polynomial selection. Policy, distribution and
micro-simulation kernels are unchanged.

Remote intermediate moments and final pipeline moments are distinct outputs.
The reported final agreement includes local numerical finishing. Replaying
saved paths verified the fix before the owner-authorized full GPU rerun.
Five scheduler/protocol tests passed, and the file comparator still detected
the known predecessor mismatch after its equal-byte fast path was added.

## Run accounting and performance scope

An SSH host-key-verification failure interrupted one final DE attempt after
eight direct evaluations. A subsequent connection verified the unchanged key;
DE resumed with separate output directories. No host-key checks were disabled.
Completed GPU batches were retained when file comparisons were consolidated
into a separate offline pass. No CPU references were rerun during recovery.
One local helper-library link attempt lacked an SDK setting and produced no
numerical output; compilation succeeded with the configured Mac SDK.

This coverage run used array dumps and is not a performance benchmark. Earlier
10-candidate integration timings—3.282 seconds for RE and 8.827 seconds for
DE—precede the final numerical-finishing layer. They remain historical
observations, not measurements of the final implementation. See
[performance scope](benchmarks.md) and the
[initial integration record](../results/cuda-validation-v020.json).

The verification instance was terminated after the final outputs and build
archive were retrieved and checked; the account query returned no instances.
Elapsed-time compute cost estimates are in the structured coverage record.
