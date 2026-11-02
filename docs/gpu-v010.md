# Historical version 0.1.0 GPU measurements

This page records the fixture-driven validation workflow measured before the
version 0.2.0 search pipeline. Statements about missing integration below
describe that earlier version. Use [the current guide](gpu.md) for execution.

# CUDA acceleration for RE and DE

The experimental FP64 backends implement the policy solve, stationary
population distribution, and required simulation outputs for both models.
RE uses its 25-state collapsed representation; DE uses 625 exogenous states.
Each has 900 endogenous states. The kernels retain the model's grids,
convergence rules, fixed draws, and observation ordering.

The interface currently consumes a CPU-exported fixture for each complete
parameter vector. The fixture supplies parameter-dependent transitions,
grids, fixed simulation inputs, and independently computed comparison arrays.
The GPU starts its own policy solve from the standard cold initialization;
CPU policy and value outputs are not used as the solver's starting solution.
This is a validation and kernel-development interface. The supplied exporter
runs a complete CPU evaluation for each new parameter vector to produce the
independent comparison arrays. The GPU runner currently requires those arrays
and compares against them after its own solve. Simply wrapping this exporter
and GPU runner inside a search would therefore solve each candidate twice.
The reported GPU timings exclude that complete CPU export, not merely file
serialization or lightweight preprocessing.

A search-ready backend needs parameter-dependent input construction separated
from CPU solving, a production evaluation mode that does not require reference
solutions, and integration with the optimizer. Input construction may stay on
the CPU; it does not have to move to the GPU to remove the duplicate solve.
The CPU direct and conditional-profile runners remain available independently.
See [Search integration](#search-integration) for the current implementation gap.

## Search integration

The GPU kernels can evaluate different parameter vectors concurrently; the
previous mixed batches establish finite evidence for that capability. A search
can exploit this when a group of candidates is available before any of their
objective values are returned. Grids, a predeclared sampling design, or a
population generated from the preceding generation have that structure.
Independent local searches can also contribute one ready candidate each.
A single adaptive local trajectory generally must wait for one evaluation
before deciding its next candidate.

The packaged `search.py` currently uses SciPy differential evolution with
`workers=1`, `updating="immediate"`, and population multiplier five. After
profiling, this gives 40 candidates for the eight-coordinate RE search and
45 for the nine-coordinate DE search. It calls the CPU profile evaluator once
per candidate. It does not call either CUDA binary.

Generation-wide batching requires a deferred population update. SciPy's
[official documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html)
specifies that only deferred updating supports parallel or vectorized
population evaluation. Changing immediate updating to deferred updating changes
which information is available when subsequent candidates are constructed;
identical seeds do not establish an identical search trajectory or endpoint.
A serial deferred implementation is the appropriate comparator for a new
batched deferred implementation. Keeping the existing immediate-update
trajectory instead requires respecting its sequential dependencies, or
parallelizing across independent search runs.

Three integration steps remain:

1. Construct parameter-dependent transitions, grids and paths without a CPU
   policy solve, and separate production evaluation from reference comparisons.
2. Connect the conditional measurement-noise moment basis and minimization to
   the GPU evaluator. Current CUDA jobs evaluate a supplied complete vector;
   they do not yet provide the profiled objective consumed by `search.py`.
3. Supply a batch objective and explicit optimizer update semantics, preserving
   candidate order, random draws, evaluation budgets and per-candidate failure
   handling. Measure performance on heterogeneous candidate batches and then
   on the declared search workload.

Per-point numerical agreement and batch throughput do not establish search
trajectory equivalence, equal evaluation counts, or end-to-end search speedup.
The cumulative throughput figures in this document remain measurements of the
GPU validation workload with pre-exported inputs. No end-to-end GPU search has
been executed or timed.

## Build and run

Requirements: Linux, an NVIDIA GPU with sufficient memory, CUDA Toolkit,
Python 3, and GNU Fortran for fixture export. The recorded A100 build used
CUDA 12.4.131 and driver 595.91.07. Source targets default to `sm_80`; set
`CUDA_ARCH` explicitly for another architecture and repeat numerical checks.

From the repository root:

```bash
CUDA_ARCH=80 bash src/cuda/build.sh RE
CUDA_ARCH=80 bash src/cuda/build.sh DE

python3 tools/export_gpu_fixture.py \
  --parameters examples/re-anchor.json \
  --build build/re-fixture --output runs/re-fixture

mkdir -p runs/re-gpu
RCC_BATCH_DUMP="$PWD/runs/re-gpu" build/cuda/re \
  --batch runs/re-fixture/job.txt --simulation-workers 4

python3 tools/compare_gpu_outputs.py \
  --cpu runs/re-fixture --gpu runs/re-gpu/job0 \
  --output runs/re-gpu/comparison.json
```

Use `examples/de-anchor.json`, separate output/build directories, and
`build/cuda/de` for DE. The exporter records source and binary identities,
complete parameter bits, and array hashes. It runs one CPU evaluation; it does
not estimate parameters. `--prepare-only` writes instrumented source without
compiling or evaluating. Fixture exports and run outputs are ignored by Git.

For batches, concatenate the `job.txt` files for the desired fixtures into a
single text file, one job per line. Each job must supply its own fixture and
all six explicit kernel parameter options. Use the generated line so that
parameters remain consistent with that fixture. Re-export after changing any
model parameter or fixed input. Batches support 1–64 jobs. Four simulation
workers are the default; eight helped the measured 48-job RE workload.
The worker count controls candidate-level finishing, not model dimensions.

Fixtures are binary little-endian streams using 32-bit integers and FP64
values. If exporting on a different host, preserve the generated files exactly
and update only the fixture directory in each job line. No numerical output
should be used to validate a different model or parameter vector.

## How the techniques work

| Technique | Mechanism | Numerical-order protection | Measured scope |
|---|---|---|---|
| Compact transitions | Store only the 25 possible successors of each exogenous state | Explicit successor indexing and ordered FP64 operations | Included in both complete backends; no isolated new timing |
| Cooperative action shortlist | A block shares independent action work and compacts candidates before final selection | Preserve tie rules and final comparisons | RE three-job policy time: 0.384 s with serial shortlist versus 0.257 s with cooperative shortlist, one observation each |
| Sparse winner and active support | Avoid work outside the exact policy-dependent dependency set | Retain required successors and the final full-state update | Included in complete-array comparisons |
| Compact Howard groups | Compute repeated selected-action successor groups in a compact layout for 99 updates, then restore the full update | Preserve each state's 25-term accumulation order | DE eight-job Howard stage: 709.20 to 468.43 ms, 1.51×; complete process: 3.222 to 3.087 s, 1.04× |
| Candidate-axis batches | Give each candidate its own state in the GPU batch dimension, sharing launches | Independent candidate stopping and buffers | RE 12 jobs: 3.069 to 1.272 s versus serial GPU, 2.41× |
| Allocation and adjacency reuse | Reuse device allocations and form stable destination adjacency | Stable source order in distribution accumulation | Included in full-backend comparisons; not an additional multiplicative claim |
| Parallel simulation finishing | Run independent candidates through separate host workers and per-thread CUDA streams | Within-candidate macro and moment loops remain ordered | DE finishing 0.591 to 0.302 s, 1.96×; RE finishing 0.466 to 0.244 s, 1.91× |
| Required-observation pruning | Compute the micro time range and macro statistics consumed by the objective; skip exactly zero population mass | Retain every required observation and its arithmetic order | Already present before the parallel-finishing change; no new incremental gain claimed |

DE-specific lag-state indexing is replaced by the actual RE transition and
macro state layout. The port does not run DE with a zero parameter in place
of the RE solver. Previously slower alternatives are not compulsory layers.
Conditional measurement-noise profiling remains a CPU toolkit feature; no
GPU-profile speedup has been measured here.

## Performance evidence

All rows below were executed on November 1, 2026 on an A100-SXM4-80GB with an
AMD EPYC 7742 host. The container's CPU quota was 13.6 CPU equivalents.
Full process time includes device-context startup, file reads, validation,
transfers and CPU finishing. It excludes CPU fixture export. The comparisons
are against preceding accelerated GPU implementations, not the original
Fortran solver. They therefore do not replace the README's M3 CPU baseline.

| Model and workload | Comparator | Selected backend | Full-process gain |
|---|---:|---:|---:|
| DE, 8 jobs, four paired observations | 3.1598 s | 2.8194 s | 1.12× |
| RE, 12 jobs, two paired observations | 3.0694 s, serial GPU | 1.2718 s, native batch / 4 workers | 2.41× |
| RE, same native 12-job binary | 1.6042 s, 1 worker | 1.2718 s, 4 workers | 1.26× |

Larger-batch observations, each measured once:

| Model | Jobs | Simulation workers | Full seconds | Evaluations/second |
|---|---:|---:|---:|---:|
| DE | 16 | 4 | 4.5143 | 3.54 |
| DE | 32 | 4 | 7.6434 | 4.19 |
| RE | 24 | 4 | 1.6295 | 14.73 |
| RE | 48 | 4 | 2.0231 | 23.73 |
| RE | 48 | 8 | 1.7064 | 28.13 |

DE batches repeat two unique vectors; RE batches repeat three. These measure
execution throughput and candidate isolation, not coverage of that many
independent parameter points. Startup, CPU scheduling and short-run variability
limit the precision of small end-to-end differences. Stage gains cannot be
multiplied without a matched combined measurement. Additional device-input
residency did not improve DE full-process timing and is not part of this build.

## Overall speed relative to the original M3 baseline

A separate short timing run on November 1, 2026 used exactly the same complete
benchmark vectors as the original-code CPU records. Four GPU processes were
executed: one single-job and one repeated-anchor batch for each model. No new
original-code CPU run was needed. All GPU jobs started their own policy solve
from the standard initialization; repeated jobs did not reuse a solved policy.

| Model | Original M3 seconds/evaluation | Accelerated M3 median seconds/evaluation | Cold GPU single job | GPU batch | Amortized GPU seconds/evaluation | Equivalent original/GPU throughput gain |
|---|---:|---:|---:|---:|---:|---:|
| RE | 873.49–1,258.27 | 0.6423 | 1.1542 s | 48 jobs / 1.6975 s | 0.03536 | 24,700–35,600× |
| DE | 792.13–794.42 | 6.9618 | 1.1592 s | 32 jobs / 7.6927 s | 0.24040 | 3,295–3,305× |

The batch rates are 28.28 RE evaluations/s and 4.16 DE evaluations/s.
The ratio is computed directly as original single-evaluation seconds divided
by observed GPU batch seconds per job; incremental speedup ratios are not
multiplied together. It compares the original serial evaluation rate with
GPU batch throughput, not a rerun of an original 48- or 32-job batch.

Keep three quantities distinct:

1. The unchanged M3-to-M3 CPU gains are approximately 1,360–1,960× for RE and
   114× for DE, using the existing three-observation accelerated medians.
2. GPU single-job cold latency was 1.1542 s for RE and 1.1592 s for DE, including
   context initialization, validation and complete array dump writes. The
   corresponding original/GPU cold-process ratios are about 757–1,090× and
   683–685×. In particular, this cold RE GPU process is slower than the recorded
   0.6423 s accelerated M3 CPU median; batching is the GPU's measured advantage.
3. Batched GPU throughput amortizes startup across independent solves. The
   24,700–35,600× and approximately 3,300× figures combine hardware, numerical
   implementation and batch parallelism. They exclude CPU fixture export and
   transfer, and have no added multiplier for reducing the search dimension.

Each of the four workloads was measured once. Original CPU observations come
from earlier dates. The batches repeat one fixed vector per model, use eight
simulation workers for RE and four for DE, and do not write the full arrays.
Single-job timings do write arrays, so their differences from batch time per
job should not be attributed solely to parallelism. These measurements are
not diverse-parameter search performance or complete parameter-to-moment
pipeline timings.

All 82 GPU evaluations passed their internal comparisons. Separately written
single-job outputs also matched 329,423 RE and 4,109,423 DE compared entries
with zero bit differences. There were no failed, stopped or discarded numerical
runs in this short timing set. The original CPU baseline was reused, not rerun.
The [structured cumulative timing record](../results/gpu-total-speed.json)
contains the observations, exact ratios, binary identities and timing boundaries.

## Numerical agreement

Independent CPU and GPU processes write separate arrays. The comparator reads
only those files, without importing either solver. The final DE mixed batch
compared **32,875,384 entries** across eight jobs with **zero bit differences**.
The final RE mixed batch compared **3,953,076 entries** across twelve jobs with
**zero bit differences**. RE's first three-point check also compared 988,269
entries with zero differences. Repetition is not additional unique coverage.

The compared objects comprise complete policy/default/value/old-value/price
arrays, stationary distribution, profit table, macro output/default/spread
arrays, all 17 moments, and the required micro investment/profit/forecast
observations. Unused micro-array slots are outside the comparison claim.
Existing comparator mutation-sensitivity evidence is reused. Exact agreement
is finite evidence for the recorded points and build; new parameters, inputs,
compilers or GPU architectures require their own checks.

The packaged computational kernels are unchanged from the measured files.
CLI parameters are now mandatory and component headers describe their role
in the complete backend. [Source identities](../results/gpu-source-identities.json)
and [aggregate measurements](../results/gpu-benchmarks.json) distinguish those
packaging edits from the measured source identities. The fixture exporter was
prepared against the packaged CPU source; its source transformation and compile
checks do not count as an additional numerical validation run.
