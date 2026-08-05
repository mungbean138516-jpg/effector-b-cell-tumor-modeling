# Validated SDE pilot results

Generated on 2026-07-28 from code commit
`af4c97e3ad4c024157c582b46821af344f217e44` with Julia 1.12.6.
The model/solver/dependency source fingerprint is
`855828c86e165bfef1eae6129ef879b8db0862fc9373955000223d41e12ca734`.

## Validation outcome

All required gates passed with zero invalid trajectories:

| Gate | Result | Acceptance rule |
|---|---:|---|
| Zero noise versus ODE | normalized error `5.12e-7` | `< 1e-3` |
| Small-noise ensemble mean | normalized error `0.0395` | `< 0.10`, all 50 paths reach day 30 |
| Seed reproducibility | exact repeat | required |
| Nonnegative domain | 0 failures in 20 paths | 0 required |
| Original-paper regression | `0.330`, 99% CI `[0.278, 0.386]` | published `0.321711` inside interval |
| Time resolution | probability range `0.07` | `<= 0.10` across fixed `dt=0.01`, fixed `dt=0.005`, and adaptive SOSRI |

## Near-fold pilot

The pilot used 200 trajectories per point, `b5=1e-4`, the paper-scale common
multiplicative noise (`epsilon=1`), a 365-day horizon, and the unguarded
Gompertz drift shared by the original CIR stochastic ODE and this repository's
formal continuation calculation.

| `b6 / b6*` | Paper persistence, 95% CI | Strict sensitivity, 95% CI | Invalid |
|---:|---:|---:|---:|
| 0.8 | `0.255 [0.200, 0.320]` | `0.160 [0.116, 0.217]` | 0 |
| 0.9 | `0.245 [0.191, 0.309]` | `0.155 [0.111, 0.212]` | 0 |
| 1.0 | `0.210 [0.159, 0.272]` | `0.115 [0.078, 0.167]` | 0 |
| 1.1 | `0.250 [0.195, 0.314]` | `0.145 [0.103, 0.200]` | 0 |
| 1.2 | `0.200 [0.150, 0.261]` | `0.110 [0.074, 0.161]` | 0 |

The intervals overlap substantially. This pilot therefore does not resolve a
monotone change in one-year establishment probability across the deterministic
fold. It also does not demonstrate stochastic switching: every trajectory
starts from a two-cell inoculation state, so the experiment measures
establishment or early fate selection.

The unit noise scale reproduces the paper's numerical convention but is not
biologically calibrated. These probabilities are conditional simulation
outcomes, not clinical predictions.

## Files

- `validation_gates.csv`: validation decisions and provenance.
- `time_resolution_estimates.csv`: all three numerical-resolution estimates.
- `paper_regression_replicates.csv`: 500 raw paper-regression trajectories.
- `near_fold_pilot_replicates.csv`: 1,000 raw five-state pilot trajectories.
- `near_fold_pilot_summary.csv`: grouped probabilities and Wilson intervals.
- `near_fold_establishment_probability.png`: PI-ready summary figure.

The earlier stochastic outputs under `results/beta5/` and `results/beta6/`
remain preliminary because incomplete numerical trajectories could be counted
as successful tumors.
