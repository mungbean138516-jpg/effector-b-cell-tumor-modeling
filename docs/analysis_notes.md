# Step 0 and Step 1 Feedback Deck

## Scope
- This deck only covers the current post-meeting work.
- Step 0 = freeze the model structure.
- Step 1 = long-time dense 1D ODE sweep in `b6`.

## Step 0
- Direct and coupled CTL formulations give nearly identical tumor-level outputs near the `b6` threshold.
- Because of that, the simpler direct CTL formulation is fixed for Step 1.
- Step 1 then keeps:
  - direct CTL formulation
  - `b5 = 1e-4`
  - baseline initial condition
  - only `b6` varied

## Step 1 results to emphasize
- The trajectories are the main evidence, because they show how the tumor path changes across `b6`.
- At `b6 ≈ 6.05e-7`, the tumor looks suppressed at day 365 but rebounds by day 10,000.
- Around `6.95e-7` to `9e-7`, the system stays in a low intermediate regime over long time.
- Clear controlled behavior begins later, around `1.15e-6`.

## Key numerical anchors
- `b6 = 1e-7`: `T365 ≈ 9.75e6`, `T10000 ≈ 9.75e6`
- `b6 ≈ 6.05e-7`: `T365 ≈ 5.76e3`, `T10000 ≈ 9.74e6`
- `b6 ≈ 6.95e-7`: `T365 ≈ 2.92e3`, `T10000 ≈ 3.02e3`
- first `T10000 < 1e5`: `b6 ≈ 6.69e-7`
- first `T10000 < 1e3`: `b6 ≈ 1.15e-6`

## Main message
- Day 365 is useful but not sufficient near the threshold.
- The 10,000-day trajectories are what distinguish delayed escape from durable low-tumor behavior.
- The switch is narrow, which is why the threshold-region trajectories are worth showing prominently.
