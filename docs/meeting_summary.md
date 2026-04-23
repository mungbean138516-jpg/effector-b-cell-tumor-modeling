# Long-Term ODE Dense b6 Sweep

## Scope
- Focused script: `4:15 after meeting modifications/step1_d2/b6_dense_step0_step1_pretty.jl`
- Model kept fixed after the meeting
- Only `b6` was varied
- Baseline used `b5 = 1e-4`, direct CTL formulation, and the baseline initial condition
- Each simulation was extended to `10,000` days

## Main findings
- Increasing `b6` steadily lowers the early tumor burden.
- The strongest qualitative change is not gradual; it happens in a narrow transition window near `b6 ~ 10^-6`.
- Below that window, higher `b6` can delay growth substantially without preventing eventual long-time escape.
- Around the threshold, a distinct intermediate regime appears with long-time tumor levels around `10^3`.
- Above the threshold, the system reaches a low-tumor controlled state.

## Useful numerical anchors
- At `b6 = 1e-7`, `T100 = 1.28e4`, `T365 ≈ 9.75e6`, `T10000 ≈ 9.75e6`
- Near `b6 = 6.05e-7`, `T365 = 5.76e3`, but `T10000 = 9.74e6`
- Near `b6 = 6.95e-7`, `T365 = 2.92e3`, `T10000 = 3.02e3`
- First case with `T10000 < 1e5`: `b6 ≈ 6.69e-7`
- First case with `T10000 < 1e3`: `b6 ≈ 1.15e-6`

## Interpretation
- A 365-day readout is not enough near the threshold because some trajectories still escape later.
- The long-time ODE view separates three qualitatively different regimes:
  - Tumor dominant
  - Intermediate or dormancy-like
  - Tumor controlled
- Long-time control is associated with large increases in B and NK levels.
- CTL is highest near the transition region rather than in the deepest control state.

## Questions to ask
- What should count as "control" in this project: a short endpoint, a long endpoint, or a threshold-crossing rule?
- Should the intermediate regime be treated as biologically meaningful dormancy?
- Is the implied threshold for `b6` biologically plausible?
- Which next scan is most valuable: `b5 × b6`, initial condition sensitivity, or stochastic robustness?

## Suggested next steps
- Run a long-time `b5 × b6` sweep
- Add SDE around the threshold window only
- Test sensitivity to initial tumor seed
- Use bifurcation or continuation analysis if the threshold remains sharp
