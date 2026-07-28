# SDE validation and near-fold pilot

This directory replaces the project's preliminary stochastic workflow with an
auditable pipeline. The deterministic continuation results remain the current
main evidence; stochastic probabilities should be interpreted only after every
validation gate passes.

## Primary source and endpoint

Kreger, Roussos Torres, and MacLean define a stochastic metastasis as
successful when a simulation starts from two tumor cells and the tumor
population never drops below one cell during a 365-day observation window.
The paper explicitly calls this a liberal definition.

- Paper: https://doi.org/10.1158/2326-6066.CIR-22-0617
- Original code: https://github.com/maclean-lab/ModelingMDSCs

The pipeline records two endpoints:

1. `paper_established`: the original one-year criterion above.
2. `strict_established`: a sensitivity endpoint requiring the paper criterion
   plus tumor burden at or above `1e3` on at least 90% of saved days during the
   last 30 days.

The strict endpoint is configurable and is not claimed to come from the
original paper.

## Why the old stochastic outputs are preliminary

The legacy beta-6 solver stopped when any population became negative, but the
summary then counted a stopped path as successful whenever its saved tumor
values remained above one. It did not check the solver return code or whether
the path reached day 365. The committed legacy CSV therefore contains many
success probabilities near one while the mean "successful" tumor remains near
its initial value of two, which is consistent with premature termination rather
than a biological result.

This pipeline instead accepts only:

- a successful full-horizon solve; or
- an explicit termination caused by the paper's tumor-removal event.

Every other solve is retained as `valid=false` and excluded from biological
probability denominators.

## Stochastic model choices

- The primary noise mode is `:paper_common`: one shared Wiener process with
  multiplicative diffusion `g_i(X)=epsilon*X_i`, matching the original
  notebook's explicit scalar Wiener process.
- `epsilon=1` reproduces the paper's numerical noise scale but is not a
  biologically calibrated value.
- `:independent` noise is available for sensitivity analysis but represents a
  different stochastic model.
- The paper regression and near-fold pilot both use the unguarded Gompertz
  drift from the original CIR stochastic ODE and the formal continuation
  calculation, so the reported `b6*` and stochastic drift refer to the same
  deterministic vector field.
- The later guard from this repository's preliminary SDE scripts remains
  available as `:legacy_guarded`, but it is not used by the validated path.
- No post-step clamping is used. Steps leaving the nonnegative domain are
  rejected, and accepted-state negativity is an explicit QC failure.

## Run order

From the repository root:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Then run:

```bash
julia --project=. test/runtests.jl
SDE_PAPER_REPS=500 SDE_RESOLUTION_REPS=100 \
  julia --project=. code/post_meeting/sde_validation/validate_sde_pipeline.jl
SDE_NSIMS=100 \
  julia --project=. code/post_meeting/sde_validation/run_near_fold_pilot.jl
```

Increase the pilot to 200 replicates per point only after the first 100-path
run has zero invalid trajectories and stable numerical-resolution checks.

## Validation gates

The validation script checks:

1. `noise_scale=0` agrees with the ODE and is independent of seed;
2. a 50-path ensemble at `noise_scale=0.025` remains close to the ODE away
   from the fold;
3. identical stochastic seeds reproduce identical saved paths;
4. accepted states remain finite and nonnegative;
5. the original four-state effects (`a10=a11=b6=0`) recover the published
   one-year establishment probability `0.321711` within a 99% Wilson interval;
6. the establishment probability changes by no more than 0.10 across fixed
   `dt=0.01`, fixed `dt=0.005`, and the primary adaptive solver.

The pilot reads `validation_gates.csv`, requires all six gates and zero invalid
trajectories, and verifies a SHA-256 fingerprint over the model, validation,
pilot, `Project.toml`, and `Manifest.toml`. Any failed, missing, or stale gate
blocks the pilot.

## Pilot design

At fixed `b5=1e-4`, the first pilot uses:

| `b6 / b6*` | `b6` |
|---:|---:|
| 0.8 | 5.3410035e-7 |
| 0.9 | 6.0086289e-7 |
| 1.0 | 6.6762544e-7 |
| 1.1 | 7.3438798e-7 |
| 1.2 | 8.0115052e-7 |

Raw replicate rows include seed, solver status, termination reason, thresholds,
final states, saved extrema, nonnegativity diagnostics, and solver rejections.
The grouped summary reports valid/invalid counts and Wilson intervals.

The validated entry point requires at least 100 replicates per point and freezes
`b5=1e-4`, `noise_scale=1`, and `dtmax=0.05`. Alternative stochastic models or
numerical settings must write to a separate output set and receive their own
validation. Each `b6` condition has a fixed seed block, so increasing from 100
to 200 replicates preserves the original first 100 trajectories.

## Interpretation boundary

Starting from two tumor cells tests stochastic establishment or early fate
selection. It does not by itself demonstrate switching between deterministic
attractors. A switching claim requires a separate experiment initialized on a
continuation-derived stable branch and a full-state basin classifier.

Likewise, the deterministic fold is not a stochastic bifurcation. Finite noise
can smooth a sharp deterministic threshold into a horizon-dependent
probability curve, and the exact-fold point is only a numerical estimate.
