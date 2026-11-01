# Algorithms for faster heterogeneous-firm model solutions

This document describes computational acceleration techniques and their
implementation in the RE and DE pipelines of *Real Credit Cycles* (2026).
Complete sources are in [`src/fortran/`](../src/fortran/). Their identities are
in [`integrated-manifest.json`](../results/integrated-manifest.json).

### State reduction: two fewer RE state coordinates

The exogenous tuple `(A[t-1], A[t], z[t-1], z[t])` becomes `(A[t], z[t])`.
Here `A` denotes aggregate productivity, `z` firm productivity, and `t` time.
Under RE (`theta = 0`), the two lag coordinates can be removed from the
policy/transition representation while retaining both current grids.

The RE transformation removes redundant lagged labels from a 625-state
exogenous representation, leaving 25 current exogenous states. The DE build
retains its 625 states. Each build retains 900 endogenous grid points.
This is a model-dependent reduction; the DE pipeline does not invoke it.

The action envelope organizes candidate actions by their dependence on
inherited debt. Sparse winner selection consumes the retained candidates
without constructing a dense winner array for every state. Compact bins
reduce indexing and clearing work. The implementation preserves the lowest
global action index when represented candidate values tie, matching the
existing dense winner-selection convention. Nonfinite candidate values retain
a dense fallback.

Implementation: `params`, `discretize`, and `politeration` in the [RE](../src/fortran/RE/diag_firms_solve.f90) and [DE](../src/fortran/DE/diag_firms_solve.f90) solvers.

### Reachable-state fixed-policy evaluation

For a fixed policy, only selected endogenous actions and the default reset
state can affect subsequent policy evaluation. The RE path restricts the
first 99 Howard updates to this reachable set. The DE path retains
successor-specific reachable pairs because its transition support depends
on the full exogenous state. The hundredth update reconstructs the complete
value array before policy improvement.

The DE implementation groups tasks by current exogenous pair and selected
action, then packs contiguous lag-state runs. Four independent accumulators
can advance together while each retains its own successor-addition order.
This exposes parallel arithmetic without combining independent sums into a
different floating-point reduction.

Implementation: `politeration` in the [RE](../src/fortran/RE/diag_firms_solve.f90) and [DE](../src/fortran/DE/diag_firms_solve.f90) solvers.

### Certified block skipping

After a full policy-improvement step, the implementation stores conservative
action-margin information for the incumbent actions. At the next step it
accounts for changes in the expected-value table. A complete 30-debt-state
block is skipped only when every incumbent in that block remains certified.
Otherwise the sparse action envelope executes normally. Reused actions still
receive value recomputation and the default comparison.

This certificate concerns action stability within the numerical policy loop.
It is not an optimizer stopping rule or a claim about global parameter-search
coverage. The RE anchor skipped 24.46% of complete debt blocks; the recorded
state-level opportunity was 24.51%, so this implementation uses complete blocks.

Implementation: `politeration` in the [RE](../src/fortran/RE/diag_firms_solve.f90) and [DE](../src/fortran/DE/diag_firms_solve.f90) solvers.

### Common-default inverse evaluation

Some DE debt shortlists share the same default continuation in all reachable
successor states. On a certified same-capital debt run, the candidate values
then differ through their initial dividends while following the same sequence
of monotone fused multiply-add operations. The implementation evaluates the
largest initial dividend and works backward through the numerical threshold
needed to recover the lowest-index winning action.

A structure-of-arrays layout processes 16 lanes at a time. Explicit boundary
verification handles floating-point rounding and ties. Scalar correction and
the ordinary candidate recurrence remain available when certification or
boundary recovery fails. The supported pipeline retains full `endognum`
shortlist capacity; it does not assume that every shortlist fits in 30 entries.

Implementation: `politeration` in the [DE solver](../src/fortran/DE/diag_firms_solve.f90).

### Buffer exchange and direct memory production

The DE Howard transformation exchanges two value buffers between sparse
updates instead of copying their complete contents at every update. The final
dense update retains the necessary copy. It preserves the reachable set,
arithmetic order, and synchronization points.

The later direct-memory refinement produces shortlist indices and dividends
at their final addresses, consumes compressed-row winners directly, and clears
only touched counters. Full-capacity storage and scalar overflow handling
remain part of the supported implementation.

Implementation: `politeration` in the [DE solver](../src/fortran/DE/diag_firms_solve.f90).

### Conditional profiling: one fewer outer search coordinate

For the supported conditional profile, let `s` denote the measurement-noise
parameter `sigma_pi`. Hold all other model parameters and fixed draws
constant. The policy solution and underlying simulation path are shared
across three evaluations of the moment vector:

```text
m0 = m(0)
mp = m(+1)
mm = m(-1)
c1 = (mp - mm) / 2
c2 = (mp + mm) / 2 - m0
m(s) = m0 + c1*s + c2*s*s
```

Here `m(s)` is the vector of simulated moments at the specified `s`; `m0`,
`mp`, and `mm` are the three evaluated vectors; `c1` and `c2` are their
componentwise linear and quadratic coefficients. The representation applies
to the implemented moment construction. It is not a general polynomial
approximation for arbitrary new statistics.

Let `p` collect the parameters other than `s`. The RE vector has nine
parameters, so `p` has eight; DE has ten parameters, so `p` has nine. Let
`J` be the indices of the moments used in the objective, `m_J(p,s)` their
simulated values, `d_J`
the corresponding fixed target vector, and `W_J` its fixed weight matrix.
Define the residual vector `r` and the conditional objective `Q` by

```text
r(p,s) = d_J - m_J(p,s)
Q(p,s) = r(p,s)' W_J r(p,s)
Q_profile(p) = min over s in [lower, upper] of Q(p,s)
```

The prime denotes transpose. `lower` and `upper` are the specified bounds on
`s`. `Q` has degree at most four in `s`. Conditional selection considers its
real stationary points and both interval endpoints. The conditional `s` is
appended to `p` to recover the complete parameter vector. Thus the expensive
outer search has **8 rather than 9 RE coordinates**, or **9 rather than 10 DE
coordinates**. `sigma_pi` remains a varying model parameter.

The conditional formulation retains the specified moments, targets, weights,
and parameter bounds. The recorded reconstruction/direct comparisons quantify
its floating-point differences; the remaining parameters still require a
policy solution. A conditional minimizer does not establish a global minimum
over those remaining parameters.

The fused traversal uses separate lane-first accumulators, preserving each
lane's observation and arithmetic order. The complete Fortran solvers emit the
three basis vectors and coefficients. `sigma_profile.py` performs conditional
selection, `run.py` provides the evaluation interface, and `search.py` supplies
the optional outer parameter search. Floating-point reconstruction can differ
slightly from a fresh direct evaluation, as reported in the accuracy guide.

Implementation: `profile_sigma_micro_moments` in the [RE](../src/fortran/RE/diag_firms_solve.f90) and [DE](../src/fortran/DE/diag_firms_solve.f90) solvers; [conditional selection](../sigma_profile.py).