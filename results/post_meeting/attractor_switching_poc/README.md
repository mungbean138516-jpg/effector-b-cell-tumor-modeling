# Attractor-Switching POC Results

This directory contains the audited outputs from the 100-replicate-per-condition
attractor-switching proof of concept at `b6 = 1.5e-6`.

- `equilibria.csv`: refined low, saddle, and high equilibria with stability.
- `experiment_configuration.csv`: fixed numerical and experimental settings.
- `switching_replicates.csv`: one row per stochastic trajectory, including
  explicit seed, solver status, endpoint states, basin classification, Git SHA,
  source fingerprint, and runtime provenance.
- `switching_summary.csv`: directional switching probabilities with 95% Wilson
  confidence intervals.
- `validation_gates.csv`: all seven pass/fail gates.
- `time_resolution_sensitivity.csv`: paired primary and finer-step outcomes.
- `representative_trajectories.csv`: seeds used in the trajectory figure.
- `attractor_switching_probability.png`: main probability result.
- `representative_switching_trajectories.png`: one switch and two non-switching
  examples, including deterministic relaxation after the noise pulse.
- `poc_summary.txt`: plain-text run summary.

See [`docs/attractor_switching_poc.md`](../../../docs/attractor_switching_poc.md)
for the bilingual interpretation, limitations, and meeting script.
