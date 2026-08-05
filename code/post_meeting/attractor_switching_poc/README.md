# Attractor-Switching Proof of Concept

This analysis tests a different stochastic question from the earlier
two-cell establishment pilot:

> Can a finite stochastic perturbation move the five-state system from one
> deterministic basin of attraction to the other?

The experiment is anchored at `b6 = 1.5e-6`, inside the bistable region found
by formal continuation. At this parameter value, the code refines and verifies
three positive equilibria before simulating anything:

- a stable low-tumor attractor;
- a saddle equilibrium on the basin boundary;
- a stable high-tumor attractor.

Each trajectory starts exactly at either stable attractor. The validated common
multiplicative Ito noise is applied for 60 days. Noise is then set to zero, and
the terminal stochastic state is relaxed under the deterministic ODE for 1,500
days. A trajectory counts as a switch only if this noise-free relaxation
converges to the opposite attractor in full five-state log-distance. A transient
crossing of a tumor threshold is therefore not enough.

## Run

From the repository root:

```bash
julia --project=. test/attractor_switching_poc_tests.jl
julia --project=. code/post_meeting/attractor_switching_poc/run_attractor_switching_poc.jl
```

The default proof of concept uses 100 trajectories per starting-attractor and
noise-level condition. Environment variables can change the ensemble size for
development runs:

```bash
POC_NSIMS=50 POC_RESOLUTION_REPS=30 \
  julia --project=. code/post_meeting/attractor_switching_poc/run_attractor_switching_poc.jl
```

## Validation boundary

The run is blocked unless the archived six-gate SDE validation still matches
the current model and dependency fingerprint. The new workflow additionally
checks equilibrium stability, zero-noise basin preservation, exact seed
reproducibility, nonnegative full-horizon paths, deterministic endpoint
classification, and switching-probability sensitivity to a finer time-step
configuration.

This is a numerical proof of concept. The noise amplitude has not been
calibrated to longitudinal biological data, so it must not be presented as an
estimate of an in-vivo switching frequency.
