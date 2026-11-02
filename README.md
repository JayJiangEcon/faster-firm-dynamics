# Faster Firm Dynamics

I'm pleased to share a set of techniques developed by GPT, with my assistance,
that substantially accelerate the solution of an important firm-dynamics model
and the search for its best-fitting parameters. The application comprises the
Rational Expectation (RE) and Diagnostic Expectation (DE) versions of the model
in Bordalo et al.'s *Real Credit Cycles*, published in the *American Economic
Review* in 2026.

For a sense of the computational scale, the authors' accompanying documentation
reports approximately two weeks to generate the complete model output from
scratch on a 2023 Mac Studio with an Apple M2 Ultra chip and 128 GB of memory.
Our CPU development and measurements use an eight-core MacBook Air M3; the
CUDA benchmarks use an NVIDIA A100-SXM4-80GB.

The baseline for the speed comparison is the original, unmodified solver
running on my MacBook Air M3. At the same fixed complete parameter vectors,
recorded original-code single-evaluation times ranged from 873.49 to 1,258.27
seconds for RE and from 792.13 to 794.42 seconds for DE. The current accelerated
direct implementations recorded medians of 0.6423 seconds and 6.9618 seconds
over three evaluations, respectively. Relative to those medians, the original
observations imply approximately **1,360–1,960× for RE and 114× for DE**, before
adding any benefit from reducing the parameter-search dimension.

With a single A100 GPU, we achieve **equivalent evaluation-throughput gains
of 24,700–35,600× for RE and approximately 3,300× for DE** relative to the
original solver on my MacBook Air M3. Pointwise numerical agreement remains
supported by the recorded comparisons, including the final pipeline checks
at 28 distinct policy vectors per model.

Beyond these gains, **profiling the measurement-noise parameter reduces the outer
search from nine to eight coordinates for RE and from ten to nine for DE.**

Preserving numerical results has been central to this work. The existing
checks include bit-for-bit agreement across millions of policy, value-function,
and price-array entries. State reduction and conditional moment reconstruction
also have separately documented comparisons with rounding differences of about
10⁻¹⁸ to 10⁻¹⁷. The documentation below gives the precise agreement achieved
by each technique.

The techniques developed in this toolkit span several layers, from Fortran
code efficiency and dynamic optimization solvers to the parameter search
space. Before generative AI, a task of this scope might have been difficult
to accomplish even through collaboration among experts from multiple fields.
Given current academic conventions, I have not listed GPT as a coauthor, but
its contribution has been immense.

I have decided to make these algorithms and their implementation publicly
available. As time and opportunities permit, I hope to continue this work by:

1. Explaining the core techniques in greater detail.
2. Conducting more comprehensive validation.

I nevertheless reserve the right to publish this work as an academic paper
in the future.

**Q: How can I use this repository to accelerate my own structural model?**

A: Simply ask your AI agent:

> Given [this repository](https://github.com/JayJiangEcon/faster-firm-dynamics),
> please analyze which ideas could be applied to the structural model we are
> currently working with.

Through this project, I hope to demonstrate the substantial potential of
generative AI to accelerate algorithms and encourage others to explore it.
My hope is that the economics community can use generative AI with confidence
to conduct better and faster economic research.

The rest of this README was written by AI, and I have not fully checked every
detail. I recommend using AI to help read and use this project. The AI-written
documentation follows.

---

**Algorithms and Implementation Techniques for Accelerating Heterogeneous-Firm Models**

Version 0.2.0

This toolkit brings together algorithms and implementation techniques for
solving heterogeneous-firm dynamic models faster while preserving numerical
results. It reduces redundant computation and memory traffic in state
representation, action selection, policy iteration, and simulation.

The current implementation and measurements use the rational-expectations
(RE) and diagnostic-expectations (DE) models in:

> Bordalo, Pedro, Nicola Gennaioli, Andrei Shleifer, and Stephen J. Terry.
> 2026. “Real Credit Cycles.” *American Economic Review* 116(4): 1274–1308.
> [https://doi.org/10.1257/aer.20211820](https://doi.org/10.1257/aer.20211820)

Toolkit author: Jiannan (Jay) Jiang.

GPT5.6 Sol and GPT-6 Astra developed the vast majority of the ideas, techniques, and implementations in this toolkit.

The package contains accelerated numerical kernels, source transformations,
recorded performance measurements, and numerical agreement comparisons.
It covers the RE and DE evaluators and conditional moment profiling.
Complete RE and DE Fortran solvers, a shared numerical library, and the fixed
inputs are included. The measured comparisons and their execution environments
are documented alongside each method.

## Contents

| Material | Entry point |
| --- | --- |
| Algorithm descriptions | [Core techniques below](#techniques-measured-gains-and-numerical-agreement), [Algorithms](docs/algorithms.md) |
| Existing timing measurements | [Performance](docs/benchmarks.md) |
| Numerical agreement and rounding differences | [Numerical accuracy](docs/accuracy.md) |
| Build and execution instructions | [Usage](docs/usage.md) |
| Complete integrated solver sources | [Fortran source directory](src/fortran/) |
| Conditional minimization and outer search | [Profile selection](sigma_profile.py), [CPU search](search.py), [CUDA search](gpu.py) |
| Source and fixed-input identities | [Integrated manifest](results/integrated-manifest.json) |
| Structured measurements | [CPU timings](results/benchmarks.json), [CPU accuracy](results/accuracy.json), [CUDA 0.2.0](results/cuda-coverage-v020.json) |
| Version 0.2.0: CUDA search pipeline | [Techniques and coverage below](#version-020-cuda-search-pipeline), [Pipeline guide](docs/cuda-pipeline.md), [CUDA source](src/cuda/) |

The implementation combines state reduction, action envelopes, sparse winner
selection, certified block skipping, reachable-state fixed-policy iteration,
buffer exchange, and fused conditional moment calculations. The algorithms
document their floating-point ordering, tie handling, and fallback paths.

## Core implementation

The complete integrated Fortran implementations are directly available:

- [RE solver](src/fortran/RE/diag_firms_solve.f90)
- [DE solver](src/fortran/DE/diag_firms_solve.f90)
- [Shared numerical library](src/fortran/base_lib.f90)
- [Fixed simulation inputs and targets](data/fixed/)

[build.py](build.py) compiles these files directly. Both solvers include their
adopted computational layers and the fused conditional moment calculation.
[run.py](run.py) supports direct evaluations and conditional `sigma_pi`
minimization. [search.py](search.py) connects the conditional evaluator to an
outer search over user-specified bounds and a declared evaluation budget.
The objective is the weighted quadratic moment distance, using moments 2–17.
[gpu.py](gpu.py) supplies the batch-oriented CUDA pipeline. Its deferred
optimizer updates differ from the CPU driver's immediate updates; they need
not follow the same search trajectory.

| Algorithm | Implementation location |
| --- | --- |
| Eliminate redundant Fortran computation | [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90): `default_value_reference_order`, `politeration` |
| Sparse transition traversal and cached screening with ordered rescoring | [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90): `expected_value_table`, `politeration`, `simulate` |
| RE state reduction | [RE solver](src/fortran/RE/diag_firms_solve.f90): `params`, `discretize` |
| Action envelopes, reachable-state updates, sparse winners and block skipping | [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90): `politeration` |
| Grouped Howard, common-default calculations, buffer exchange and direct action-data handling | [DE solver](src/fortran/DE/diag_firms_solve.f90): `politeration` |
| Conditional moment basis and fused traversals | [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90): `simulate`, `profile_sigma_micro_moments` |
| Conditional selection and CPU outer search | [sigma_profile.py](sigma_profile.py), [search.py](search.py) |
| CUDA batch search and worker orchestration | [gpu.py](gpu.py), [cuda_pipeline.py](cuda_pipeline.py) |
| Input-only construction | [prepare_cuda_inputs.F90](src/fortran/prepare_cuda_inputs.F90) |
| Mac numerical finishing | [finalize_numeric.f90](src/fortran/finalize_numeric.f90) |
| CUDA model solution and parallel simulation | [RE backend](src/cuda/re_native_parallel.cu), [DE backend](src/cuda/de_native_simulation_parallel.cu) |

The [integrated manifest](results/integrated-manifest.json) binds the packaged
Fortran sources and fixed inputs. The [measurement source record](results/measurement-sources.json)
identifies the implementations underlying the existing measurements. See
[Usage](docs/usage.md) for direct builds, conditional evaluation, and search.

## Quick start

To display the included timing measurements, run from this directory:

```bash
python3 report_results.py
```

This reads existing records and does not execute a model. Python utilities use
only the standard library for direct builds and evaluations. Conditional
minimization uses NumPy; the optional outer search also uses SciPy. Building
an evaluator requires GNU Fortran with OpenMP. See the [usage guide](docs/usage.md)
for complete commands. No additional model-source download is required.

## Two kinds of dimension reduction

Two techniques reduce the dimension of a computational problem, in addition
to making its individual operations faster:

| Reduction | Coordinates before → after | What is saved | Existing evidence |
| --- | --- | --- | --- |
| RE state reduction | `(A[t-1], A[t], z[t-1], z[t])` → `(A[t], z[t])`: **4 → 2 exogenous state coordinates** | Removes two lag coordinates from policy and transition computations; **625 → 25 states** with five nodes per coordinate | 28-point workload: **29.36×** wall ratio; maximum moment difference **5.20e-18** |
| Conditional profiling of `sigma_pi` | RE: **9 → 8 outer search coordinates**; DE: **10 → 9 outer search coordinates** | Solves the policy once for the remaining coordinates, then handles `sigma_pi` through a cheap conditional polynomial objective | Existing profile/direct comparisons: maximum moment differences **6.25e-17 (RE)** and **1.73e-18 (DE)** |

Here `A` is aggregate productivity and `z` is firm productivity; `t` denotes
time. The state reduction applies to RE with `theta = 0`. Both current
productivity grids and the capital/debt grids retain their full resolution.

The parameter reduction is *profiling*, not fixing `sigma_pi` to a constant:
for each choice of the remaining parameters, its conditional value is selected
within the specified interval. The full parameter vector is then reconstructed.
This package supplies the moment-basis evaluator, conditional minimizer, and
an outer parameter-search driver. The detailed mechanism is in technique 13 below. Its end-to-end
search-time benefit depends on the evaluations avoided; the measured fusion
speedups in technique 14 concern the evaluator itself.

## Techniques, measured gains, and numerical agreement

All numbers below come from existing measurements recorded on September 1–5,
2026. These CPU measurements were reused without rerunning them. The recorded
CPU environment used Apple M3 and GNU Fortran 14.2.0. Later optimized
comparisons used eight OpenMP threads; the first technique has the historical
thread-setting limitation described below. For the later comparisons,
the build configuration and fixed-input identities are in the
[measurement source record](results/measurement-sources.json).

A speedup is predecessor time divided by successor time: values above 1 mean
faster execution. **Wall time** measures elapsed runtime; **stage time** measures
one computational component. Each row has its own predecessor and workload.
Historical batch ratios and interleaved timing pairs have different exposure
to machine load. These factors must not be multiplied into an overall speedup.
An *anchor* is a fixed benchmark parameter vector. A moment is a simulated
summary statistic; Q is the numerical objective computed from those moments.

| Technique | Model | Existing timing evidence | Agreement evidence |
| --- | --- | --- | --- |
| Eliminate redundant Fortran computation | RE / DE | Hoist-only historical DE ≈792 → 305 s; RE ≈1,258 → 611 s; thread effects not isolated | Initial hoist: all 17 moments bitwise equal at each model's published point |
| Traverse only nonzero transition successors | RE / DE | 625 → 25 successor visits before RE collapse; no isolated wall factor | Included in enclosing finite refactor comparisons; DE retains 625 states |
| Cache screening scores and rescore retained actions | RE / DE | Shared component of envelope/sparse-winner gains; no additional multiplier | Original-order rescoring with the current near-tie guard; finite implementation evidence |
| Collapse redundant states | RE | 1,902.493 → 64.801 s; **29.36×** over 28 points | Maximum moment difference 5.20e-18 |
| Update only reachable states | RE / DE | RE **1.1768×**; DE anchor **2.18–2.27×** | Bitwise moments; separate complete-array comparisons |
| Construct an action envelope | RE / DE | RE **1.7748×** over 28 points | Bitwise moments and Q; tied winners retained |
| Select winners from sparse lists | RE / DE | RE **1.1730×**; DE eight-point workload **1.3656×** | Complete arrays and output comparisons |
| Skip certified debt blocks | RE / DE | RE **1.0598×**; DE two-point workload **1.0774×** | Complete arrays and output comparisons |
| Group Howard updates in four accumulators | DE | **1.1911×** median in a congested five-pair run; noisy | 14,062,500 array entries, zero mismatches |
| Reuse common-default calculations | DE | About **1.03× action-stage** gain; isolated wall gain unestablished | 2,812,500 array entries, zero mismatches |
| Exchange value buffers | DE | **1.1542×** over eight points | Complete arrays and iteration histories |
| Write and consume action data directly | DE | **1.0275×** wall; **1.0785×** action stage | 3,375,000 entries plus overflow fallback |
| Profile out one search coordinate | RE / DE | **9 → 8 / 10 → 9** outer coordinates; one policy solve for all three basis vectors | Profile reconstruction checked against direct evaluation |
| Fuse the three moment traversals | RE / DE | RE **1.7692×**; DE **1.1578×** | Basis vectors and coefficients bitwise equal |

The following sections explain what changes, why it saves work, and what was
actually compared. [Technique evidence](results/technique-evidence.json)
contains the additional earlier-stage extracts and source-record identities;
[timing records](results/benchmarks.json) and
[accuracy records](results/accuracy.json) contain the later measurements.

### 1. Eliminate redundant Fortran computation

**Principle.** Our first optimization removes repeated calculations from the
original Fortran implementation. In the policy solver, the continuation value
following default depends only on the exogenous state and the current input
value array, yet the original code recomputes it for every one of the 900
endogenous grid points. We move this calculation outside the endogenous-state
loop and reuse its result within each value update. It is recomputed when the
input value array changes. The initial transformation retains the original
625-successor summation order and also removes concurrent OpenMP writes to
the same continuation-value entry.

The broader implementation also precomputes payoff components at their actual
dependency level and skips the terminal `distribution()` call when policy
iteration already produced that distribution. The `doVFI=0` path retains the
call. Payoff reuse feeds the screening/rescoring mechanism in technique 3;
it does not replace final action arithmetic with a reassociated formula.

**Recorded effect.** The early comparison records approximately **792 → 305 s
for DE** and **1,258 → 611 s for RE** at the published parameter vectors on
the same machine. These rounded historical observations are reused, not new
measurements or interleaved estimates of an isolated speedup. The original
solver used a hard-coded 36-thread setting; the retained summary does not
specify the hoist-only successor's thread setting, so these times do not
separate the effects of hoisting and thread allocation. These observations
apply to the initial hoist, not to the entire expanded redundancy-removal group.

**How agreement was checked.** The historical record reports all 17 simulated
moments bitwise equal to the corresponding baseline for each model. This is
fixed-point evidence, not a guarantee over the parameter domain. More
aggressive early caches changed floating-point summation order and discrete
policy choices; those failed variants are excluded from this technique's
agreement and timing claims.

Implementation: `default_value_reference_order` and its calls in `politeration`
in the [RE](src/fortran/RE/diag_firms_solve.f90) and
[DE](src/fortran/DE/diag_firms_solve.f90) solvers. These integrated solvers also
contain later optimizations; the early timings do not measure the complete
current implementation. The historical evidence extract and its limits are
retained in [technique-evidence.json](results/technique-evidence.json).

### 2. Traverse only nonzero transition successors

**Principle.** In the original four-coordinate exogenous representation, the
next lagged coordinates are fixed by the current state. Only the two next
current productivity coordinates vary, giving 5 × 5 = 25 possible successors
among 625 labels. The optimized expectation and simulation loops enumerate
those 25 successors directly, retaining their original order. DE still stores
625 exogenous states; this does not perform the RE state collapse described
in technique 4. It also differs from restricting Howard updates to endogenous
actions reachable under the current policy.

**Work saved and evidence.** The affected successor loops visit 25 rather than
625 entries. This is a loop-work reduction, not a measured 25-fold wall-speed
increase. Existing comparisons of the enclosing refactors support finite
numerical agreement; the retained records do not isolate this component's
wall-time contribution. After RE state collapse, the same two-coordinate
successor enumeration covers its 25-state representation directly.

Implementation: `expected_value_table`, `default_value_reference_order`, and
successor loops in `politeration`, `distribution`, and `simulate` in the
[RE](src/fortran/RE/diag_firms_solve.f90) and
[DE](src/fortran/DE/diag_firms_solve.f90) solvers.

### 3. Cache screening scores and recompute retained actions in original order

**Principle.** Compute continuation expectations once per exogenous state and
candidate action at each policy-improvement step. Reuse factorized payoff components to
form inexpensive screening scores. The action envelope then retains a small
candidate set, including near ties under its implemented guard. Retained
candidates use the original dividend expression and ordered successor
recurrence; represented ties select the lowest global action index. Later
common-default acceleration has its own verified inverse-recurrence path.

This separates cheap candidate selection from expensive final action scoring.
A cached sum accumulated from zero cannot simply replace the recurrence that
starts from the dividend: the two can round differently and select different
discrete actions. Earlier attempts that made this substitution changed results
and are excluded from the retained evidence.

**Evidence and scope.** This shared mechanism is part of the current envelope
and sparse-winner implementations. Their finite comparisons are reported in
techniques 6 and 7; no separate wall-time multiplier is assigned to this
component. The historical interval-guard predecessor used explicit score
intervals, whereas the current envelope uses bins and a `1e-12` near-tie guard.
They are different selection constructions. Ordered rescoring preserves a
retained action's evaluation; it alone does not prove that screening retained
every possible winner over all parameter values. Historical interval timings
and their limits appear in the supporting-techniques table below.

Implementation: `expected_value_table`, payoff-component arrays, and candidate
selection/rescoring inside `politeration` in the
[RE](src/fortran/RE/diag_firms_solve.f90) and
[DE](src/fortran/DE/diag_firms_solve.f90) solvers.

### 4. Remove two RE state coordinates

**Principle.** When the diagnostic-expectations parameter `theta` is zero,
policies and transitions depend on current aggregate and firm productivity,
not their lag labels. The transform replaces 625 exogenous labels with 25
economically distinct states, a 25-fold reduction in state count. It updates
transition and simulation lookups consistently and marginalizes the initial
distribution. The 30 × 30 endogenous grid, fixed shocks, simulation size,
and convergence controls stay at their full settings. This reduction uses
the RE structure; the DE implementation retains its larger state space.

**Measured effect.** The same 28-point workload took 1,902.493 s with our earlier accelerated
625-state interval-guard implementation and 64.801 s after collapse:
**29.36×** by the ratio of historical batch totals. This is an observed batch
comparison, not an interleaved estimate of a portable speedup.

**How agreement was checked.** Across 476 moment cells, 402 were bitwise equal.
The maximum absolute difference was **5.2042e-18** in moments and **3.4694e-17**
in raw Q; the complete ranking of the 28 points was unchanged. Thus this
comparison supports numerical agreement with small observed rounding
differences, rather than bitwise identity to the 625-state computation.
Only the full-settings batch supplies these figures.

Implementation: [RE solver](src/fortran/RE/diag_firms_solve.f90), `discretize`.

### 5. Update only states reachable under the current policy

**Principle.** Howard iteration repeatedly evaluates a fixed policy. Most
actions are never selected by that policy, so their intermediate values need
not be updated on every push. The RE implementation updates the reset action
and the actions selected by continuing firms during the first 99 pushes,
then reconstructs all 900 action values on the 100th push before policy
improvement. DE refines this further: it tracks reachable pairs of a successor
exogenous state and an action, rather than a global union of actions. Both
retain the successor summation order.

**Measured effect.** RE used an average of 233.02 live actions out of 900.
Its 28-point wall total fell from **64.8007 to 55.0649 s (1.1768×)**; policy
time fell from **52.9307 to 43.3545 s (1.2209×)**. For DE, the anchor wall
comparisons were **22.160 → 10.158 s** and, with execution order reversed,
**21.924 → 9.645 s**, giving **2.18× and 2.27×**. A separate four-point
sequential workload took **110.820 → 40.042 s (2.768×)**. The different designs
are retained as separate observations.

**How agreement was checked.** RE retained all 476 moment cells, all 28 Q
values, and their ranking bit for bit. DE retained all 527 moment cells across 31 fixed points. A further eight-point
DE comparison, with points fixed in advance using a scrambled Sobol design,
matched **22,500,000 complete-array entries**, 136 moments, parameter bits,
and both objective outputs. That wider comparison spans the envelope and
reachability layers together. Later RE stepwise comparisons also checked
complete arrays downstream of state collapse, as described below.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 6. Build an action envelope before evaluating candidates

**Principle.** Actions share a common dependence on inherited debt. The
envelope uses this structure to identify the debt thresholds at which actions
can become competitive. Binary search locates thresholds; compact bins
collect potentially winning actions. The implementation preserves all
potential tied winners and recomputes retained actions using the original
scalar arithmetic order. This avoids repeatedly scanning the entire action
grid for every debt state while preserving discrete tie selection.

**Measured effect.** In the RE 28-point workload, wall time fell from
**55.0649 to 31.0254 s (1.7748×)** relative to the active-update predecessor.
Policy time fell from **43.3545 to 19.0993 s (2.2699×)**. The packaged DE
implementation also uses an envelope; this README does not assign it the RE
factor or infer an isolated DE factor from a combined comparison.

**How agreement was checked.** All 476 RE moment cells, 28 Q values, and the
complete ranking were bitwise equal. Later complete-array comparisons checked
the RE envelope against its active-update predecessor. The DE eight-point
comparison in the preceding section covers the envelope together with
successor-local updates.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 7. Select winners directly from sparse lists

**Principle.** Once the envelope produces a short candidate list, filling a
dense array and scanning it with `MAXLOC` reintroduces unnecessary work.
Sparse selection visits the candidate list directly, preserves the lowest
index among tied winners, and retains a dense fallback for nonfinite cases.
Compact bins reduce the work needed to assemble those lists.

**Measured effect.** Three interleaved RE anchor pairs reduced median wall
time from **0.7862 to 0.6703 s (1.1730×)** and policy time by **1.2894×**.
For DE, an eight-point workload fell from **151.4816 to 110.9266 s (1.3656×)**.
A separate DE profile measured **1.3606× wall** and **1.9511× action-stage**
gains; its three samples used the anchor and adjacent floating-point values
of the measurement-noise parameter, not three identical repetitions.

**How agreement was checked.** RE matched all 476 moments over 28 points and
all 112,500 entries across five arrays in a separate-process anchor comparison.
DE matched 527 moments over 31 points, and its eight-point comparison retained
**22,500,000 array entries** and 136 moments without mismatches. A temporary
DE source mutation changed moments, Q, and four arrays; a fresh execution of
the unchanged source restored equality, demonstrating comparator sensitivity.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 8. Skip debt blocks only when the incumbent is certified

**Principle.** A previous winning action can remain optimal after continuation
values change. A conservative bound compares its margin with the possible
change in competing values. A whole 30-debt-state block is skipped only when
every state passes that check. The chosen value and default decision are
still updated; uncertified blocks use the ordinary envelope calculation.

**Measured effect.** Three RE anchor pairs gave **0.6807 → 0.6423 s (1.0598×)**
wall time and **1.0431×** policy time, skipping 24.46% of blocks. DE skipped
34.46% at the anchor and 34.18% at an adjacent-float point. Its order-balanced
anchor result was **1.2057×**, its two-point workload result **1.0774×**, and
its instrumented action-stage result **1.2208×**. These are separate workloads,
not interchangeable estimates of one speedup.

**How agreement was checked.** RE retained the 28-point moment/objective
outputs and the separate anchor arrays. DE retained **5,625,000 array
entries** across the two points, with byte-equal output files. The RE
stepwise array comparisons also include the transition from sparse selection
to certified skipping.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 9. Group DE Howard updates without reassociating sums

**Principle.** Reachable pairs can be grouped by current productivity and
selected action, then processed in contiguous runs of lag states. Four
independent accumulators keep intermediate values in registers and reuse
transition data. Each accumulator still receives its 25 successor terms in
the same order; this does not replace the sum with a different reduction tree.

**Measured effect.** The initial three interleaved pairs gave **1.2158×** by
the ratio of median wall times, **1.0945×** by the median paired ratio, and
**1.3828×** in value-push time. A separate three-point sequential batch was
slower: **23.4065 → 24.4902 s (0.9558×)**. A later five-pair run under heavy
host contention had median wall times **9.66 → 8.11 s (1.1911×)**, but paired
ratios ranged from **0.7002× to 1.2604×**. These measurements show a useful
kernel-level opportunity with noisy elapsed-time gains; they do not establish
a stable wall-time factor.

**How agreement was checked.** Five predefined points, including an anchor,
an adjacent-float point, and three fixed Sobol points, retained
**14,062,500 array entries** and all output files. A temporary change to the
first accumulator was detected in values, policies, and moments; restoring
the tracked source restored byte equality.

Implementation: [DE solver](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 10. Reuse common-default calculations with independent lanes

**Principle.** Some contiguous debt runs share a capital choice and default
in every successor state. They therefore share a monotone continuation-value
recurrence. The implementation evaluates the upper candidate and works
backward through the recurrence to identify the earliest tied winner, using
exact boundary checks. A structure-of-arrays layout processes 16 independent
lanes while preserving the arithmetic within each lane. Full shortlist
capacity and a scalar fallback handle cases outside the compact fast path.

**Measured effect.** Recorded action-stage times were **1.7498 and 1.7425 s**,
compared with a fresh profiled predecessor median of **1.7989 s**: approximately
**1.03×** at that stage. These observations do not establish a reliable
isolated wall-time gain.

**How agreement was checked.** The anchor comparison retained all
**2,812,500 entries** across five arrays and every compared output file.
The maintained implementation uses full-capacity scratch storage. Separate
profile-versus-direct checks show small rounding differences, which are
reported under conditional profiling below.

Implementation: [DE solver](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 11. Exchange value buffers instead of copying them

**Principle.** Consecutive fixed-policy pushes alternate between old and new
value arrays. Exchanging the two allocatable buffers avoids copying the full
array after each sparse push. The final dense copy remains in place, with
the same reachability logic, arithmetic, and synchronization.

**Measured effect.** Eight fixed DE points took **63.8672 → 55.3335 s (1.1542×)**
in total; policy time improved by **1.1802×**. A separate conditional-profile
comparison gave **1.1320×** over three interleaved pairs. Its order-specific
ratios differed substantially, so the profile factor is less stable than the
eight-point direct-workload result.

**How agreement was checked.** The existing direct and profile comparisons
retained numerical outputs and iteration histories. The profile comparison
retained all **2,812,500 entries** across five arrays, including policy indices,
default decisions, value functions, and bond prices, as well as the moment
basis and selected profile outputs.

Implementation: [DE solver](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 12. Write and consume action data directly

**Principle.** Shortlists and dividends are written directly to their final
buffers. Winner selection consumes sparse residual data without copying it
back into temporary dense storage, and only counters that were touched are
cleared. Bounded scratch storage has an explicit scalar overflow path.
These changes reduce memory traffic without changing candidate construction
or arithmetic order.

**Measured effect.** Three timing brackets, each with predecessor runs around
the candidate, gave **1.0275× wall** and **1.0785× action-stage** gains. The
wall factor is the ratio of summed bracket-reference times to summed candidate
times. It measures this action-data refinement relative to the buffer-exchange
predecessor. The later profile agreement check supplies no additional timing
claim.

**How agreement was checked.** The direct comparison retained
**3,375,000 entries** across six complete arrays and all compared outputs.
A deliberately reduced scratch capacity forced **606,177 overflow blocks**
and still retained equality. Two separate profile processes also matched
3,375,000 entries and all compared outputs.

Implementation: [DE solver](src/fortran/DE/diag_firms_solve.f90), `politeration`.

### 13. Profile out one search coordinate with a shared policy solution

**Principle.** The measurement-noise scale `sigma_pi` enters measured profit
after policies have been solved. With other parameters and shocks fixed,
the relevant moment vector is quadratic in this scale. Let `s = sigma_pi`
and let `m0`, `m+`, and `m-` denote moment vectors evaluated at `s = 0, +1, -1`.
Then, component by component,

```text
c1 = (m+ - m-) / 2
c2 = (m+ + m-) / 2 - m0
m(s) = m0 + c1*s + c2*s*s
```

With the other parameters held fixed, a weighted squared-moment objective
is therefore a polynomial of degree at most four in `s`. Conditional selection
checks stationary points and interval endpoints within the specified bounds.
The outer search varies the other **eight RE** or **nine DE** parameters;
conditional selection supplies `sigma_pi` to reconstruct all **nine RE** or
**ten DE** parameters. Only one parameter-search coordinate is removed by
this step.

The evaluator solves the policy once and emits the three basis vectors and
coefficients. This permits repeated conditional moment calculations without
re-solving the economic model. The packaged `run.py` selects the conditional
value from this basis, and `search.py` supplies the outer parameter search.

**Measured effect.** The structural saving is one shared policy solution for
three basis evaluations. No isolated wall-time factor for shared-path reuse
is quoted here. The next section reports the separately measured improvement
from fusing traversals in an already shared-path implementation.

**How agreement was checked.** Reconstructed moments were compared with direct
evaluation of the same complete parameter vector. The maximum absolute moment
differences were **6.2450e-17 (RE)** and **1.7347e-18 (DE)**; the corresponding
raw-Q differences were **0** and **2.7756e-17**. Polynomial reconstruction
and direct evaluation can round differently even when the underlying
conditional formula is the same.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `profile_sigma_micro_moments`.

### 14. Fuse the three conditional moment traversals

**Principle.** The three noise-scale lanes traverse the same simulated path.
A fused traversal loads shared quantities once and updates three independent
sets of accumulators. Each lane preserves observation order and its original
arithmetic, so reuse does not require a different summation order.

**Measured effect.** Three interleaved pairs against the earlier profile
implementation reduced RE median wall time from **1.15 to 0.65 s (1.7692×)**
and simulation time from **0.7180 to 0.2419 s (2.9686×)**. DE median wall time
fell from **8.44 to 7.29 s (1.1578×)** and simulation time from **1.8772 to
0.8985 s (2.0892×)**. One DE wall-time pair was slower, at **0.9364×**. These
measurements precede later DE memory refinements and describe traversal fusion
relative to its own profile predecessor.

**How agreement was checked.** The before/after comparisons retained the three
basis vectors, their coefficients, and recorded conditional candidate and
selection outputs bit for bit. This is distinct from the small rounding
differences between polynomial reconstruction and direct evaluation above.

Implementation: [RE](src/fortran/RE/diag_firms_solve.f90) and [DE](src/fortran/DE/diag_firms_solve.f90), `profile_sigma_micro_moments`.

## Supporting implementation techniques and historical variants

The main table describes mechanisms in the packaged CPU solvers. The following
items clarify their implementation and earlier development. Historical variants
are not additional modes supplied by this package, and their timing ratios
must not be multiplied into the current pipeline's performance.

| Technique | Role and availability | Existing evidence and limits |
| --- | --- | --- |
| Controlled threading and contiguous memory layout | Packaged: runtime OpenMP thread count, action-contiguous expectation/price arrays, factorized payoff arrays, and loop unrolling | Replaces the original hard-coded 36 threads; recorded optimized comparisons normally use eight. These changes accompany other refactors; no isolated speed factor is assigned. |
| Interval-guard action screening | Historical predecessor; not a separate packaged mode. Retains actions whose score upper bound reaches the largest lower bound, then rescores them in original order | At the wrapper point, RE 158.049 → 36.696 s (about 4.31×); DE 178.827 → 77.886 s (about 2.30×). Historical, non-interleaved comparisons against conservative refactors; finite agreement, not a portable interval proof or a current-envelope timing. |
| Moment-only simulation pruning | Historical optional CPU variant, default off and not included as a packaged mode. Firm paths stop after period 58; the aggregate distribution still spans 302 periods | RE: all 476 moments and 28 raw-Q values matched over 28 points. A stabilized hot-host comparison reduced simulation medians 0.4345 → 0.3936 s (9.4%); the approximately 2.4% total saving is a projection. DE retained 85 moments at five points without an isolated speed gain. |

Thread count is set through `run.py --threads`; `build.py` records the compiler
and flags. Contiguous arrays reduce indirect access and organize repeated
loads without licensing reassociated floating-point reductions. The failed
aggressive compiler/cache variants are not part of these retained techniques.

The interval predecessor has separate finite evidence: the RE historical
28-row comparison retained all 476 moment bits; the DE five-point comparison
retained 85 moment cells and five raw-Q values against its stated comparators.
Those records do not establish a domain-wide guarantee for the interval bound,
and they must not be attributed to the current envelope's different guard.

Moment-only pruning removes unused diagnostics as well as unused firm periods.
The firm moments use periods 51–56 and leads through period 58, while aggregate
moments retain the original macro path and GDP construction. That variant
cannot supply the omitted diagnostic outputs. Its finite agreement is retained
as historical evidence, without changing the packaged simulation mode.

See [performance scope](docs/benchmarks.md#foundational-and-optional-techniques),
[accuracy scope](docs/accuracy.md#foundational-and-optional-techniques), and
[structured evidence](results/technique-evidence.json).

## What supports the result-preservation claim

The evidence combines fixed-input moment comparisons, complete-array
comparisons, and separately executed predecessor/successor programs. The
arrays cover policies (`kbpolind`), default indicators (`dfpolind`), value
functions (`Vtheta`, `Vthetaold`, and where applicable `Vndtheta`), and bond
prices (`qmat`). Temporary source mutations in existing comparisons checked
that the comparator detected altered calculations; fresh executions of the
unchanged source restored agreement.

For RE, a later stepwise check ran five stages in distinct processes: state
collapse, reachable updates, action envelope, sparse winner, and certified
skipping. Four adjacent comparisons at the anchor checked **450,000 array
entries** in total with zero mismatches. Two additional points fixed in advance
from independent Sobol designs added eight adjacent comparisons and
**900,000 entries**, again with zero mismatches. These comparisons begin at
the collapsed implementation; they do not convert the earlier 625-to-25-state
rounding comparison into a bitwise result.

The counts describe entries compared across points or adjacent stages, not
independent random observations. Numerical agreement is established for the
recorded inputs and execution identities. The evidence distinguishes bitwise
equality from small measured rounding differences and does not claim a
universal guarantee across parameter values, compilers, or architectures.
See [numerical accuracy](docs/accuracy.md) for the evidence map and
[performance](docs/benchmarks.md) for timing designs and other measured variants.

## Version 0.2.0: CUDA search pipeline

**A complete CUDA-oriented pipeline, from parameter search to core model solution.**

A deferred population search proposes a batch of candidates. A persistent input
worker constructs their grids, transitions and fixed-draw paths; the CUDA worker
solves their policies and distributions and runs the required simulations.
Local numerical finishing and conditional measurement-noise minimization then
return an ordered objective vector to the optimizer. The workers remain available for the next generation.
Each candidate starts an independent cold model solve.

The pipeline combines the existing GPU techniques with:

- Input construction that requires no preliminary CPU policy solve.
- Production evaluation without CPU reference solution arrays.
- One policy solve for the three conditional moment-basis evaluations.
- Mac numerical finishing that preserves the macro statistic and direct-objective arithmetic.
- Ordered batch results, explicit evaluation budgets and deferred search updates.
- SSH execution with input construction, conditional selection and search on the Mac.

The measured Mac-to-A100 route preserves the tested Mac input and conditional
selection paths. Across 28 distinct policy vectors per model, the twelve input
files per point,
all direct-output entries checked, and the conditional basis/selection outputs
were bit-identical to their CPU references. See [numerical coverage and execution
identities](docs/cuda-validation-v020.md) for the exact scope.

| Final version 0.2.0 numerical comparison | RE | DE |
|---|---:|---:|
| Different policy vectors checked | 28 | 28 |
| Direct-output entries checked | 9,223,844 | 115,063,844 |
| Direct-output bit mismatches | 0 | 0 |
| Direct objectives and conditional selections | All bit-identical | All bit-identical |

The point set includes the original anchor, nine earlier operational-test
candidates, sixteen points from two independent Sobol scrambles and two boundary
checks. The same independently executed CPU references were reused throughout.
The expanded checks found and resolved last-bit differences in numerical
finishing. [The verification report](docs/cuda-validation-v020.md) states the
point-selection rules, arithmetic changes and complete comparison scope.

Earlier small-search timings are retained as [historical integration
measurements](docs/benchmarks.md#initial-pipeline-timing-before-final-numerical-finishing).
The final coverage run is not a performance benchmark or a full-search speedup
measurement.

Start with the [execution-path guide](docs/usage.md#choose-an-execution-path),
then follow the [CUDA pipeline instructions](docs/cuda-pipeline.md). The entry
point is [`gpu.py`](gpu.py); the persistent workers and transport are managed
by [`cuda_pipeline.py`](cuda_pipeline.py).

### 1. Batch independent candidates across the GPU

Add a candidate dimension to the CUDA work grid so that multiple independent
model solves share kernel launches. Each candidate retains its own buffers,
stopping flags and cold initialization. The recorded RE 12-job comparison
fell from 3.0694 seconds with serial GPU execution to 1.2718 seconds with the
native batch and four finishing workers (**2.41×**). This is a combined backend
comparison, not an isolated launch-overhead measurement.

### 2. Retain CUDA workers and device allocations

Keep the CUDA context, workers and allocation pool alive between search
batches. This avoids repeating startup and allocation work for each generation;
it does not reuse a different candidate's solved policy. There is no isolated
speed multiplier for persistence in the final pipeline. Ordered responses and
separate candidate storage preserve evaluation identity.

### 3. Construct inputs without solving the model on the CPU

A persistent Mac worker builds parameter-dependent grids, transitions,
exogenous paths and action caches. Fixed draws and their normal transforms are
loaded once; parameter-dependent inputs are rebuilt for every candidate.
This removes the preliminary CPU policy solve from the production path.
All **336 input files per model** matched the independent CPU references in
the 28-point coverage. No separate incremental timing is claimed.

### 4. Compact transitions and cooperate on action selection

Store only the 25 possible transition successors and let GPU threads share
independent action-screening work. Compact shortlists, sparse winners and
active support reduce unnecessary action and state work. Final selection
preserves the established tie rules and ordered FP64 accumulation. In an early
RE three-job comparison, cooperative shortlisting reduced policy time from
0.384 to 0.257 seconds (**1.49×**, one observation per variant). Compact
transitions and sparse support have no separate additive speed claim.

### 5. Group repeated Howard updates

Represent repeated selected-action successor groups compactly for 99 Howard
updates, then restore the full update. Each state's 25-term accumulation order
is retained. For the recorded DE eight-job workload, Howard time fell from
709.20 to 468.43 milliseconds (**1.51×**); complete-process time fell from
3.222 to 3.087 seconds (**1.04×**). The stage gain is not the whole-model gain.

### 6. Preserve distribution order while reusing adjacency

Build stable destination adjacency and reuse allocations for distribution
updates. This reduces setup and memory-management work while preserving the
source order of each accumulation. Its numerical behavior is covered by the
complete-array comparisons; no isolated speed multiplier was measured.

### 7. Prune unused simulation work and finish candidates concurrently

Compute the observation range required by the objective and skip exactly
zero population mass. Retain every consumed observation in its original
arithmetic order. Independent candidates then finish through separate host
workers and per-thread CUDA streams. Recorded finishing times fell from
0.591 to 0.302 seconds for DE (**1.96×**) and from 0.466 to 0.244 seconds for
RE (**1.91×**). Observation pruning was already present in the comparator;
these gains measure parallel finishing, not pruning.

### 8. Share one policy solve across conditional measurement-noise evaluations

Use the same policy and macro path for the zero, plus-one and minus-one
measurement-noise lanes. Their moments define the polynomial minimized on
the Mac over the declared interval, leaving eight outer RE coordinates and
nine outer DE coordinates. The final checks cover all basis coefficients,
candidate objectives and selected outputs. No additional GPU-profile speed
multiplier has been measured.

### 9. Preserve the Mac arithmetic in a small numerical finishing step

Return the already computed 302-entry macro output path with each CUDA
response. A small resident Mac Fortran library calculates growth variance and
the direct weighted objective using the reference log implementation and
reduction order. This resolved the last-bit differences found by the wider
checks without moving policy solution or simulation back to the CPU. The
final 28-point comparisons per model include this step; its incremental
runtime has not been measured separately.

### 10. Connect batches to deferred population search

Freeze the population while proposing and evaluating one generation, then
apply updates after the ordered batch results return. Explicit budgets and
candidate ordering connect the optimizer to persistent workers. This exposes
parallel work to the GPU; it changes scheduling relative to immediate-update
search and does not promise the same optimization trajectory. No final
end-to-end search speedup is inferred from the kernel throughput figures.

The [GPU technique record](docs/gpu-v010.md#how-the-techniques-work) supplies
the historical incremental measurements. These gains use different workloads
and predecessors and must not be multiplied together. The
[final verification report](docs/cuda-validation-v020.md) covers the integrated
Mac → A100 → Mac pipeline, including complete solver and simulation outputs.

## Substantial contribution from GPT

**GPT, used through Codex, made substantial contributions to algorithm
development, code implementation, performance analysis, numerical verification
work, and technical documentation. Its contribution went well beyond language
editing.**

Jiannan (Jay) Jiang directed the work, made the scientific decisions, and
takes responsibility for this toolkit and its documentation. GPT is
acknowledged for its substantive contribution but, given current academic
conventions, is not listed as an author in the citation metadata. OpenAI is the
provider of the AI system.

## License and citation

Original contributions are provided under the [Academic Citing License](LICENSE).
The existing template requires scientific publications reporting results
obtained with the software to cite the specified reference. This is a
nonstandard license; OSI approval is not claimed. [NOTICE](NOTICE) explains
the template's origin and attribution of the underlying model software.

The reference for this toolkit is:

> Jiang, Jiannan (Jay). 2026. *Faster Firm Dynamics: Algorithms and Implementation
> Techniques for Accelerating Heterogeneous-Firm Models.* Version 0.2.0.
> Software and accompanying methods, performance, and numerical accuracy
> documentation.

[CITATION.cff](CITATION.cff) contains the corresponding machine-readable
metadata. Please also cite the original paper when using its model.
