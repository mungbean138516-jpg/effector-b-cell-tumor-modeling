# Effector B-cell Tumor Modeling Reviewer Guide

## Project in one paragraph

This project extends a published tumor-immune dynamical system by adding an
effector B-cell compartment. It asks how B-cell-mediated tumor control and
MDSC-mediated B-cell suppression alter long-term tumor behavior. The analysis
progresses from deterministic parameter sweeps to equilibrium continuation,
local stability analysis, basin mapping, a two-parameter saddle-node boundary,
validated stochastic establishment experiments, and an attractor-switching
proof of concept. The principal result is a model-defined tipping structure,
not a clinical prediction or treatment recommendation.

## Inputs model and outputs

| Component | What is used in this repository |
| --- | --- |
| Biological system | Five modeled populations: tumor cells (`T`), MDSCs (`M`), NK cells (`N`), CTLs (`C`), and effector B cells (`B`) |
| Main interaction parameters | `b6` represents B-cell-mediated tumor control; `b5` represents MDSC suppression of B cells |
| Data input | Literature-based equations and parameter assumptions; no patient-level clinical dataset is used |
| Generated data | Deterministic and stochastic simulation trajectories, equilibrium branches, parameter grids, basin classifications, validation records, and replicate-level outcomes |
| Deterministic model | Coupled nonlinear ordinary differential equations implemented in Julia and SciML |
| Stochastic model | Multiplicative Ito SDE extension with explicit solver, seed, termination, and nonnegativity checks |
| Main outputs | CSV tables, trajectory plots, heatmaps, continuation curves, stability diagnostics, basin maps, switching probabilities, confidence intervals, and validation gates |
| Intended use | Mechanistic and computational investigation of tumor-control regimes and loss of control |
| Not supported | Patient-specific prognosis, an in-vivo event rate, a clinical indication, a treatment schedule, or a health-policy recommendation |

The clinical research experience described in the author's application
materials is separate from this repository. No protected Keck Medicine patient
information is included or used in the model.

## Research sequence

1. **Model extension.** Add effector B cells and biologically motivated
   interactions to the Tumor-MDSC-NK-CTL system.
2. **Screening.** Use long-time parameter sweeps to identify a narrow change in
   tumor behavior as `b6` increases.
3. **Dynamical explanation.** Refine equilibria, compute Jacobian spectra, and
   use pseudo-arclength continuation to distinguish a true saddle-node from a
   finite-time plotting threshold.
4. **Dependence on context.** Map basins of attraction and continue the fold in
   the `b5`-`b6` plane to test how initial state and MDSC pressure alter the
   boundary.
5. **Stochastic validation.** Rebuild the SDE outcome pipeline after finding
   that incomplete trajectories could be misclassified as successful.
6. **Switching proof of concept.** Initialize at verified attractors, apply a
   finite noise pulse, turn noise off, and classify the final basin after
   deterministic relaxation.

## Verified findings

- At fixed `b5 = 1e-4`, a simple real eigenvalue crosses zero near
  `b6* = 6.676e-7`, supporting a generic saddle-node interpretation in the
  current model and parameterization.
- Long-time behavior is not determined by `b6` alone. Initial tumor, NK-cell,
  and B-cell levels can place the system in different basins at the same
  parameter value.
- Increasing MDSC suppression of B cells shifts the fold so that stronger
  B-cell tumor-control potency is required to reach the low-tumor regime.
- The corrected two-cell SDE pilot produced one-year establishment estimates of
  approximately `0.20` to `0.255` across the tested fold neighborhood. The
  confidence intervals overlap, so this experiment does not establish a
  monotone stochastic threshold.
- In the separate attractor-switching proof of concept at `b6 = 1.5e-6`, the
  observed low-to-high switching proportion increased from `0.03` at
  `epsilon = 0.10` to `0.64` at `epsilon = 0.30`. No high-to-low event was
  observed in 100 runs per tested noise level; this does not prove that reverse
  switching is impossible.

## Why the validation correction matters

The preliminary stochastic workflow stopped when a population became negative
but could still classify the last saved tumor value as a successful outcome.
The replacement workflow requires a valid full-horizon solve or a documented
tumor-removal event, records solver and termination metadata, excludes invalid
paths from probability denominators, reproduces the published four-state
benchmark, and blocks analysis unless validation gates and source fingerprints
pass. This correction changed the interpretation of the stochastic results and
is a central example of why outcome definitions and numerical quality control
must be designed with the model rather than added after a result appears.

## Interpretation and real-world relevance

The current work is deliberately methodological. It shows how a mechanistic
model can generate apparently sharp conclusions that change after longer
simulation, equilibrium analysis, initial-condition testing, or solver-quality
checks. Its immediate contribution is therefore a reproducible framework for
reasoning about tipping behavior, uncertainty, and model validity.

Any clinical or policy use would require a separate next stage: identify a
well-defined population and decision, obtain appropriate longitudinal data,
estimate or constrain parameters, evaluate identifiability, validate on
independent data, quantify uncertainty and subgroup behavior, and involve
clinical experts before interpreting the outputs as actionable evidence.

## Recommended reading order

1. [`README.md`](../README.md) for the repository map and headline results.
2. [`presentations/post_meeting/june_update/June_update.pdf`](../presentations/post_meeting/june_update/June_update.pdf) for the deterministic bifurcation and basin story.
3. [`docs/june_update_figure_provenance.md`](june_update_figure_provenance.md) for the slide-to-code and slide-to-output audit.
4. [`docs/sde_validation_update.md`](sde_validation_update.md) for the corrected stochastic pipeline and its limitations.
5. [`docs/attractor_switching_poc.md`](attractor_switching_poc.md) for the latest proof of concept, bilingual meeting script, and likely questions.
6. [`code/post_meeting/`](../code/post_meeting/) and [`results/post_meeting/`](../results/post_meeting/) for reproducible sources and outputs.

## Reproducibility status

The repository includes a locked Julia environment, committed source and raw
replicate outputs, fixed random seeds, solver metadata, validation-gate tables,
tests, source fingerprints, and figure provenance. Earlier stochastic figures
under `results/beta5/` and `results/beta6/` are retained as preliminary history
and should not be cited as validated evidence.
