# July 2026 SDE Validation Update

This note documents the corrected stochastic validation update presented to the MacLean Lab group on July 28, 2026. It separates the reported scientific result from files that are still needed for full computational reproduction.

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

- `presentations/post_meeting/sde_validation_update/Effector_B_cell_SDE_Meeting_Update.pptx`: original six-slide update deck shared with the lab.
- `presentations/post_meeting/sde_validation_update/Effector_B_cell_SDE_Meeting_Update.pdf`: PDF rendering included for GitHub preview.
- `results/post_meeting/sde_validation/near_fold_establishment_pilot.png`: the result panel extracted from slide 5 of the source deck.

## Reproducibility status

The corrected SDE implementation, validation test suite, and raw pilot CSV were not present in the available MacLean Lab workspace, the archived GitHub package, the two supplied shared conversations, or the Slack attachments reviewed for this update. They have not been reconstructed or replaced with the older stochastic script.

In particular, `code/core/b6_sde_threshold_refined.jl` predates the corrected validation pipeline and should not be cited as the generator of the July result.

To complete the reproducible package, locate and add:

1. the corrected SDE simulation and event-classification script;
2. the 26-test validation suite;
3. the raw 1,000-trajectory pilot output and summary CSV;
4. the exact noise configuration, solver settings, and random seeds.
