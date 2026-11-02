# Performance measurements

The main comparison set was recorded on September 4–5, 2026; foundational
and optional historical observations below also cover September 1–2.
They are reused for this toolkit. No new numerical evaluation was run
during its preparation.

The [README technique guide](../README.md#techniques-measured-gains-and-numerical-agreement)
also covers the initial elimination of redundant Fortran computation,
state reduction, reachable updates, action envelope,
DE sparse selection, certified skipping, grouped Howard, and common-default
measurements from September 1–4. Their selected observations are stored in
[technique-evidence.json](../results/technique-evidence.json). The table below
and `report_results.py` retain the later comparison set; they are not a
complete list of all README techniques.

## Environment and measurement design

- Host: Apple M3, macOS 15.7.7.
- Compiler: GNU Fortran 14.2.0.
- Flags: `-O3 -funroll-loops -fopenmp -m64`.
- Parallel execution: eight OpenMP threads, one evaluator at a time.
- Inputs: the fixed files identified in [measurement-sources.json](../results/measurement-sources.json).

Each speedup is predecessor wall time divided by successor wall time. Each
row has its own predecessor; the ratios must not be multiplied into an
overall speedup. Source identities for the current pipelines are separate
from the historical implementation stages measured in these comparisons.

The numerical controls were 50 default iterations, 20 policy iterations,
100 Howard updates, and at most 1,500 distribution iterations, with
distribution tolerance `1e-6` and policy tolerance `1e-4`. These settings do
not by themselves identify the stopping condition of an individual run.

## Initial elimination of redundant computation

The first README technique reuses the historical hoist-only summary: DE
approximately 792 → 305 seconds and RE approximately 1,258 → 611 seconds,
with all 17 moments bitwise equal at each model's published point. These are
rounded, same-machine observations, not a fresh or interleaved benchmark.
The original solver hard-coded 36 threads; the retained hoist-only summary
does not identify the successor's thread setting. The timing comparison
therefore does not isolate hoisting from thread allocation, and the eight-thread
configuration above must not be imputed to this early successor.

The extract is stored with its source-record identity in
[technique-evidence.json](../results/technique-evidence.json). Later integrated
solver timings include additional transformations and are separate evidence.

## Foundational and optional techniques

Sparse transition traversal reduces affected loops from 625 labels to 25
nonzero successors before RE state collapse. Cached screening and ordered
rescoring also contribute to the enclosing solver comparisons. Neither has
an isolated wall-time factor in the retained evidence, so neither adds a
multiplier to the main table. Likewise, the initial hoist timings above do
not measure every payoff precomputation or removal of a repeated distribution
call now described under redundancy elimination.

The historical interval-guard wrapper comparisons were RE 158.049 → 36.696 s
(about 4.31×) and DE 178.826807 → 77.885571 s (about 2.30×), relative to
conservative refactors at the corresponding complete vector. They were not
interleaved timing experiments. These implementations preceded the packaged
action envelope and do not supply its isolated performance or proof.

Optional moment-only simulation pruning had a stabilized hot-host RE simulation
median of 0.4345 → 0.3936 s (about 9.4% less stage time). Applying that reduction
to the earlier simulation share, 7.858 / 31.025 seconds, projects approximately
2.4% less total evaluator time. That projection is not an observed isolated
end-to-end gain. Other sequential comparisons were contaminated by host drift;
the DE port established finite agreement without an isolated speed gain.
The optional variant remains default off historically and is not a mode in
this package.

Runtime-controlled threading, contiguous action arrays, factorized payoff
storage, and `-funroll-loops` are supporting implementation choices in the
packaged path. The original hard-coded 36-thread setting was removed. No
separate speedup is assigned to each choice; numerical agreement is bound to
the recorded compiler and execution environment.

Selected evidence and historical-source extracts are in
[technique-evidence.json](../results/technique-evidence.json).

## Main comparisons

| Computational change | Before | After | Wall speedup | Measurement |
| --- | ---: | ---: | ---: | --- |
| RE sparse winner | 0.786222 s | 0.670290 s | 1.172958× | Three interleaved pairs, anchor medians |
| RE certified block skipping | 0.680668 s | 0.642251 s | 1.059816× | Three interleaved pairs, anchor medians |
| DE buffer exchange, direct | 63.867202 s | 55.333482 s | 1.154223× | Totals across eight fixed points |
| DE buffer exchange, profile | 41.597714 s | 36.745880 s | 1.132038× | Totals across three pairs at one profile point |
| DE direct memory refinement | Three bracket references | Three candidate runs | 1.027493× | Ratio of summed bracket-reference time to candidate time |
| RE fused conditional profiling | 1.15 s | 0.65 s | 1.769231× | Three interleaved pairs, anchor medians |
| DE fused conditional profiling | 8.44 s | 7.29 s | 1.157750× | Three interleaved pairs, anchor medians |

Complete vectors, timing brackets, and available observations are included in
[benchmarks.json](../results/benchmarks.json). That file is read directly by
`report_results.py`.

### State and policy kernels

For RE sparse winner and block skipping, the policy-stage ratios were
1.289354× and 1.043058×. Their timing baselines are distinct. The block-skipping
comparison reused complete certified blocks rather than changing the action
selection rule within a partially certified block.

The DE buffer-exchange direct workload contains eight predetermined points.
Each point was evaluated once by each implementation, with alternating order.
The profile workload repeats one fixed profile point in three interleaved
pairs. The wider record, including the profile-array comparison, contains
24 executions, ten distinct complete parameter vectors, and 218.901657
seconds of summed evaluator wall time.

The profile workload's order strata had ratios 1.403619× and 0.715911×.
Timing noise and execution order therefore matter even though its aggregate
ratio was positive. The eight-point direct workload is finite workload
evidence, not a survey of every parameter configuration.

### Direct memory refinement

The selected DE memory refinement had wall ratio 1.027493× and action-stage
ratio 1.078491×. It was measured in three interleaved rounds with neighboring
predecessor runs providing each timing reference. All three brackets are
retained in the structured measurements.

The development comparison also measured four other memory variants. Across
the full comparison there were 39 evaluator calls, one complete parameter
vector, and 257.656748 seconds of summed process wall time. The small wall gain
is reported as measured; it is not rounded into a larger performance claim.

### Fused conditional profiling

| Model | Predecessor observations | Successor observations | Simulation-stage speedup |
| --- | --- | --- | ---: |
| RE | 1.28, 1.12, 1.15 s | 0.77, 0.65, 0.65 s | 2.968590× |
| DE | 7.22, 9.66, 8.44 s | 7.71, 7.29, 6.71 s | 2.089243× |

RE simulation medians were 0.718007 and 0.241868 seconds. DE simulation
medians were 1.877229 and 0.898521 seconds. One DE paired wall ratio was
0.936446×, so the successor was slower in that pair despite its positive
median result. The DE fusion measurements preceded the later buffer and
direct-memory changes and must not be relabeled as timings of the final
combined implementation.

## Other measured variants

| Variant | Recorded performance | Use in the current pipeline |
| --- | --- | --- |
| RE grouped Howard | Median policy 0.398623 → 0.431809 s; wall 0.652182 → 0.672778 s | Comparison only |
| DE standalone SoA direct production | Wall ratio 0.920455×; stage ratio 0.958111× | Comparison only |
| DE old static Howard control | Wall ratio 0.955616×; stage ratio 0.940966× | Control only |
| DE cleaned Howard | Wall ratio 0.996811×; stage ratio 0.967949×; unstable brackets | Comparison only |
| DE guided Howard scheduling | Wall ratio 1.003889×; stage ratio 0.990439× | Comparison only |
| DE path-only reuse | 7.73 → 7.68 s | Superseded by fused traversal |
| DE lane-last fusion | 7.73 → 16.82 s | Comparison only |

The numerical comparisons for the SoA and Howard memory variants retained
their recorded finite agreement, but performance favored the selected
pipeline. These observations concern the implementation variants, host, and
workload described here.

The development records also contain four setup failures: an RE source-pattern
adaptation before compilation, a DE compiler-environment failure before object
creation, a DE argument-parsing failure before evaluator launch, and a DE
input-read failure after 0.15 seconds in the compile directory. They produced
no numerical result. Runtime was not recorded for the first three failures.

## Reading the evidence

`source_id` links each structured measurement to a dated record hash in
[provenance.json](../results/provenance.json). Numerical values are retained
from those records. The included data provide the relevant vectors and timing
observations; hashes identify their source records and do not establish
numerical validity by themselves. See [Numerical accuracy](accuracy.md) for
the separate agreement evidence.

## Context for the preface

The records below support the first-person preface and provide additional
computational workload context. All figures come from existing measurements.
[Preface evidence](../results/preface-evidence.json) retains the numerical
extracts and the source-record identities. No new evaluator or search was run
for this documentation update.

The accompanying March 2026 computational documentation for *Real Credit
Cycles* identifies a 2023 Mac Studio with an Apple M2 Ultra, 128 GB RAM, and
macOS Sequoia 15.6.1. Its approximately two-week estimate concerns complete
model-output generation without the precomputed intermediate simulations.
### Original-code baseline on the same laptop

The headline ratios use actual original-source executions on the eight-core
MacBook Air M3 with 16 GB RAM, compared with current direct evaluations at the
same complete anchor vectors. The original solver hashes match the upstream
snapshots in the source manifest. Complete parameter bits were checked against
the original wrapper and the accelerated timing records. Fixed shocks, targets,
weights, full grids, simulation settings, and convergence controls are retained.

| Model | Original-source wall time | Current direct wall time | Original/current ratio |
| --- | ---: | ---: | ---: |
| RE | 1,258.2728 s | 0.642250791 s | 1,959.16× |
| DE | 792.1285 s | 6.961757583 s | 113.78× |

Each original time is one recorded direct evaluation. Each accelerated time is
the median of three recorded direct evaluations at the matching vector. The
RE successor is the state-collapsed action-envelope/sparse-winner/certified-
skipping implementation; the DE successor includes the promoted action-direct
memory refinement. Its three wall observations are 5.293285875, 7.492780666,
and 6.961757583 seconds. These timings use the direct evaluator, with no
conditional `sigma_pi` search.

The original source retains its configured 36 threads; the accelerated code
uses eight. GNU Fortran 14.2.0 was used on this laptop. The ratios therefore
measure the complete operational implementation improvement, including thread
configuration. The original and current timings were recorded on different
dates and are not interleaved original/current pairs. Fanless-laptop thermal
state and machine load can change the observed ratios. The same source report
also records another original RE observation of approximately 873.49 s and
another original DE observation of approximately 794.42 s. Relative to the
accelerated medians above, the original observations imply about 1,360–1,960×
for RE and 113.8–114.1× for DE. The preface includes that observed baseline
variation. This range is not a confidence interval and does not measure all
variation in the accelerated implementation.

The 1,902.492909-second RE batch belongs to our earlier accelerated 625-state
interval-guard implementation. Likewise, the 2,341.804856-second DE batch
belongs to our conservative full-action refactor. Their 61.3× and 15.0× ratios
remain development-stage comparisons; they are not the original-code baseline
for the preface. The detailed method sections retain these separate comparisons.

No matched total-search timing isolates the additional benefit of profiling
out a parameter. A reduction in coordinate count alone does not determine
that speedup for an adaptive search. The finite bitwise comparisons and the
small state-collapse/profile reconstruction differences are documented in
[Numerical accuracy](accuracy.md).
