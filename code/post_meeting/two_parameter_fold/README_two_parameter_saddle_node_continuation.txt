Two-parameter saddle-node continuation in b5-b6
================================================

Scientific purpose
------------------
Continue the previously classified saddle-node bifurcation while varying:

    b5 = MDSC inhibition of effector B cells
    b6 = effector-B anti-tumor strength

The result is a formal bifurcation boundary:

    b6_fold(b5)

This upgrades the earlier simulation-defined b5 × b6 phase boundary into a
steady-state bifurcation result.

Required prior output
---------------------
The script needs:

    formal_continuation_low_seed.csv

It also uses, when available:

    zero_eigenvalue_classification.csv
    step3_control_boundary.csv

Run
---
From the repository root:

    include("setup_two_parameter_continuation.jl")

Then:

    include("continue_saddle_node_in_b5_b6.jl")

Main output
-----------
    two_parameter_fold_continuation_outputs/
        saddle_node_threshold_in_b5_b6_plane.png
        saddle_node_curve_b5_b6.csv
        two_parameter_special_points.csv

Interpretation
--------------
If the formal fold curve rises with b5:

    stronger MDSC suppression requires greater B-cell anti-tumor strength
    before the low-tumor equilibrium becomes dynamically available.

If the quadratic coefficient b20 crosses zero:

    the fold curve may contain a cusp bifurcation.

The script asks BifurcationKit to detect cusp, Bogdanov-Takens, and Zero-Hopf
points along the fold curve.

Important terminology
---------------------
The plotted curve is a saddle-node / fold bifurcation boundary.

It is not merely:
    - a trajectory threshold,
    - an endpoint cutoff,
    - or a heatmap contour.

It is the parameter locus where a stable and an unstable equilibrium meet.
