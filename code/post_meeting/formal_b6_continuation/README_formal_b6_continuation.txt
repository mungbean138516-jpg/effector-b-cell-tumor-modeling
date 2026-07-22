Formal b6 continuation
======================

Purpose
-------
Use BifurcationKit.jl pseudo-arclength continuation to trace equilibrium
branches across b6, including turning points that ordinary parameter sweeps
cannot follow.

This is the next formal validation after:
1. dense b6 threshold sweep,
2. initial-condition robustness,
3. b5 × b6 phase diagram,
4. up/down hysteresis screen,
5. basin-of-attraction maps.

Setup
-----
From the repository root, run once:

    include("setup_bifurcation_environment.jl")

Then run:

    include("b6_formal_palc_continuation.jl")

Required existing input
-----------------------
The script looks for:

    b6_equilibrium_branches.csv

in the current folder or common repository result folders. This file supplies
high- and low-tumor equilibrium seeds that already passed Newton and residual
checks.

Main output
-----------
    formal_b6_continuation_outputs/equilibrium_branches_across_b6.png

Interpretation
--------------
- solid/filled stable branch points: max real eigenvalue < 0
- open unstable branch points: max real eigenvalue > 0
- red stars: detected fold points
- a stable–unstable–stable branch structure supports saddle-node bistability
- no detected unstable branch does not prove there is none; step sizes,
  parameter range, and continuation settings may need refinement

Package basis
-------------
The script follows the current BifurcationKit workflow:
ODEBifProblem + PALC + ContinuationPar, with fold and bifurcation detection.
