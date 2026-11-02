# Numerical accuracy

The comparisons below are existing results for the RE and DE numerical
implementations. They use fixed inputs and the recorded CPU/compiler
environment. No new numerical evaluation was run to prepare this package.

## Recorded agreement

| Implementation comparison | Scope | Recorded outcome |
| --- | --- | --- |
| RE sparse winner and block skipping | 28 points; 476 moment cells; complete parameter bits and both objective outputs | Zero bit mismatches |
| RE complete-array anchor | Five arrays; 112,500 entries; separate processes | Zero bit mismatches; output files equal |
| DE buffer-exchange profile | Five arrays; 2,812,500 entries | Zero bit mismatches |
| DE direct memory refinement | Six arrays; 3,375,000 entries at the fixed anchor | Zero bit mismatches; output tags equal |
| DE scalar overflow fallback | Deliberately constrained scratch capacity; 606,177 overflow blocks | Complete-array and output equality retained |
| DE profile memory refinement | Six arrays; 3,375,000 entries; two separate processes | Zero bit mismatches; output tags equal |
| RE and DE fused profiling | Three moment-basis vectors and their linear/quadratic coefficients | Predecessor/successor bitwise equality |

The five common arrays are `kbpolind`, `dfpolind`, `Vtheta`, `Vthetaold`, and
`qmat`. The DE direct-memory comparisons also include `Vndtheta`. The arrays
record policy indices, default indicators, value functions, and bond prices.
The structured evidence retains individual array sizes and hashes in
[accuracy.json](../results/accuracy.json).

The forced-overflow comparison deliberately exercises a scalar fallback;
its reduced scratch capacity is a verification configuration, not the
supported production capacity. The supported DE pipeline retains full
shortlist capacity.

## Earlier computational layers

The [README](../README.md#techniques-measured-gains-and-numerical-agreement)
explains each layer together with its timing baseline. Selected measurements
and source-record identities are in
[technique-evidence.json](../results/technique-evidence.json).

| Comparison | Coverage | Recorded outcome |
| --- | --- | --- |
| Default continuation hoisted out of the endogenous-state loop | Published parameter point for each of RE and DE; historical summary | All 17 moments bitwise equal per model; does not cover the rejected reassociated caches |
| RE 625-to-25 state reduction | 28 points, 476 moment cells | 402 cells bitwise equal; maximum moment difference 5.204170427930421e-18; maximum raw-Q difference 3.469446951953614e-17; complete ranking unchanged |
| RE reachable updates and action envelope | 28 points per comparison | 476 moments and 28 raw-Q values bitwise equal |
| RE stepwise layers downstream of collapse | Four adjacent comparisons at one anchor; eight more at two predefined points | 450,000 plus 900,000 complete-array entries; zero mismatches |
| DE reachability, including prior envelope layer | Eight predefined Sobol points, separate processes | 22,500,000 complete-array entries, 136 moments, parameters and objectives bitwise equal |
| DE sparse winner | Eight predefined points | 22,500,000 complete-array entries and 136 moments; zero mismatches |
| DE certified skipping | Two points | 5,625,000 array entries and complete outputs equal |
| DE grouped Howard updates | Five predefined points | 14,062,500 array entries and complete outputs equal |
| DE common-default calculation | One anchor | 2,812,500 array entries and complete outputs equal |

The RE stepwise comparisons start with the collapsed implementation. Their
bitwise agreement does not imply bitwise equality between the earlier
625-state representation and the collapsed representation. Complete-array
counts sum entries over points or adjacent comparisons.

Existing DE comparisons also include temporary source-mutation checks for
reachability, sparse selection, and grouped Howard updates. Their comparators
read independently produced outputs, detected the source changes, and found
agreement again after fresh execution of the unchanged source.

## Foundational and optional techniques

Sparse successor enumeration and cached screening with ordered rescoring are
components of the enclosing refactors, not new isolated numerical comparisons.
DE retains its 625-state representation while visiting 25 nonzero successors.
The current action envelope uses its bin-based near-tie guard; rescoring a
retained action in original order does not by itself prove that no winner
was excluded. The existing finite envelope and sparse-winner comparisons
supply the current implementation's evidence.

The historical RE interval guard retained all 476 moment bits over 28 stored
rows. Its recomputed raw Q differed from the stored-moment recomputation by
at most 3.47e-17, with the tested rank ordering unchanged. The historical DE
interval candidate retained all 85 moment cells and five raw-Q values against
its stated comparators at five fixed points. These are distinct predecessor
implementations, not a proof for the current guard or another architecture.

Optional moment-only simulation pruning retained 476 moment cells and all
28 raw-Q values for RE, and 85 moments and five raw-Q values for DE. The
comparison concerns moment outputs; omitted diagnostics are not reproduced.
The current package does not expose this historical optional mode. Its
performance limitations are recorded in [Performance](benchmarks.md#foundational-and-optional-techniques).

The first technique's original 17-moment comparisons concern the initial
hoist. Additional payoff reuse and removal of the repeated terminal distribution
call belong to the enclosing refactor evidence; they are not independently
validated by that hoist-only comparison.

## Which dimensions are removed

RE state collapse removes the lagged aggregate- and firm-productivity state
coordinates: four exogenous coordinates become two. The current productivity
and endogenous grids are unchanged. The state-count and rounding comparisons
above apply to this representation change.

Conditional `sigma_pi` profiling removes one coordinate from the outer
parameter search: RE 9 → 8 and DE 10 → 9. It selects `sigma_pi` conditionally
and reconstructs the complete vector. The comparisons below verify moment
reconstruction at the same complete vector; they do not measure the total
runtime of an outer search. State-space and parameter-search reductions refer
to different computational objects.

## Conditional profile and direct evaluation

The conditional-profile implementation constructs the moment basis at
`sigma_pi = 0, +1, -1`. Its coefficient mapping is documented in
[Algorithms](algorithms.md). Predecessor/successor agreement for the fused
traversal is bitwise in the recorded comparisons.

Reconstructing a moment vector from the coefficients and evaluating the same
complete parameter vector directly can produce small rounding differences:

| Model | Maximum absolute moment difference | Absolute raw-objective difference |
| --- | ---: | ---: |
| RE | 6.245004513516506e-17 | 0 |
| DE | 1.734723475976807e-18 | 2.7755575615628914e-17 |

These are different comparisons: bitwise equality between two implementations
of the same profile calculation does not imply bitwise equality between
polynomial reconstruction and a separately executed direct calculation.

## Interpretation and execution identity

Complete-array equality, finite multi-point agreement, and timing comparisons
answer different numerical questions. They are kept separate here. The
reported comparisons establish agreement for their recorded configurations;
they are not guarantees across every parameter vector, compiler, or
architecture.

[measurement-sources.json](../results/measurement-sources.json) retains the
historical source identities underlying the measurements above.

The complete packaged Fortran sources and fixed inputs have their own
[integrated manifest](../results/integrated-manifest.json). The integration
checks in [integrated-validation.json](../results/integrated-validation.json)
record successful RE/DE compilation, equality of numerical subroutine source
bodies with the recorded fused implementations, Python entry-point checks,
and analytic examples for the conditional polynomial minimizer.
No model evaluation or parameter estimation was executed for this integration.
The new Python search interface has not undergone an end-to-end model search.
Historical timing and agreement observations are not fresh measurements of
these new entry points. Source comparison and compilation alone do not
establish numerical agreement across new parameter vectors or environments.
