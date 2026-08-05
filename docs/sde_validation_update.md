# July 2026 SDE Validation Update

This note documents the corrected stochastic validation update presented to the MacLean Lab group on July 28, 2026. The source, tests, raw trajectories, summary outputs, and presentation are included for computational reproduction.

## Why the earlier result changed

The earlier near-100% stochastic "success" rate was a numerical artifact. Some trajectories terminated before day 365 after a state became negative, and the incomplete solve was classified using its last stored tumor value. That allowed failed solves to enter the success count.

The rebuilt workflow instead:

- requires a valid trajectory through day 365;
- logs the termination reason and excludes invalid paths;
- uses the original CIR paper establishment criterion as the primary endpoint;
- reports a stricter endpoint as a separate sensitivity analysis;
- checks the zero-noise limit, small-noise behavior, seed repeatability, nonnegativity, and time resolution;
- reproduces the published four-state stochastic benchmark before testing the five-state B-cell extension.

## Reported validation and pilot results

The source deck reports:

- 26 of 26 validation tests passed;
- zero invalid trajectories in the final pilot;
- a four-state benchmark estimate of 0.330 versus the published value of 0.3217;
- a small-noise mean error of 0.0395, below the stated 0.10 acceptance threshold;
- 1,000 valid near-fold trajectories across `b6 / b6* = 0.8, 0.9, 1.0, 1.1, 1.2`;
- establishment probabilities of approximately 0.200-0.255 under the paper criterion;
- establishment probabilities of approximately 0.110-0.160 under the stricter criterion.

All 95% confidence intervals overlap, and the estimated response is non-monotonic. The pilot therefore does not resolve a stochastic threshold near the deterministic fold.

## Interpretation

The formal saddle-node remains the primary result. This pilot measures tumor establishment from a two-cell initial condition; it does not test noise-induced switching between the high- and low-tumor attractors. The next mechanistic SDE experiment should initialize near the deterministic attractors and track crossing or switching events explicitly.

## Included files

- `code/post_meeting/sde_validation/EffectorBSDE.jl`: model, solver, event handling, trajectory QC, and endpoint classifiers.
- `code/post_meeting/sde_validation/validate_sde_pipeline.jl`: six-gate validation workflow.
- `code/post_meeting/sde_validation/run_near_fold_pilot.jl`: gated five-point near-fold pilot.
- `test/runtests.jl`: 26 regression and failure-classification tests.
- `results/post_meeting/sde_validation/`: validation gates, 500 paper-regression trajectories, 1,000 near-fold trajectories, grouped estimates, and the final figure.
- `presentations/post_meeting/sde_validation_update/Effector_B_cell_SDE_Meeting_Update.pptx`: original six-slide update deck shared with the lab.
- `presentations/post_meeting/sde_validation_update/Effector_B_cell_SDE_Meeting_Update.pdf`: PDF rendering included for GitHub preview.
- `results/post_meeting/sde_validation/near_fold_establishment_probability.png`: the validated result panel used in slide 5.

## Reproducibility status

The original implementation was recovered from the shared conversation's Git bundle. Its five-commit history is preserved in this repository, ending at `bcb62b2` (`Add validated near-fold SDE results`). The output files identify code state `af4c97e3ad4c024157c582b46821af344f217e44`, source/dependency fingerprint `855828c86e165bfef1eae6129ef879b8db0862fc9373955000223d41e12ca734`, and Julia 1.12.6.

The current `code/core/b6_sde_threshold_refined.jl` is a compatibility entry point that directs validated work to the new pipeline. The pre-correction implementation is retained separately as `code/core/legacy_b6_sde_threshold_refined.jl` and must not be cited as the generator of the July result.

To reproduce the checks and pilot from the repository root:

```bash
julia --project=. test/runtests.jl
SDE_PAPER_REPS=500 SDE_RESOLUTION_REPS=100 \
  julia --project=. code/post_meeting/sde_validation/validate_sde_pipeline.jl
SDE_NSIMS=200 \
  julia --project=. code/post_meeting/sde_validation/run_near_fold_pilot.jl
```

The pilot refuses to run if validation is missing, failed, or stale relative to the code/dependency fingerprint.
