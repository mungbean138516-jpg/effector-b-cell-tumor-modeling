# June Update Figure Provenance

This audit maps the research figures in `presentations/post_meeting/june_update/June_update.pdf` (22 pages, created June 27, 2026) to the archived Julia code and repository outputs used to produce them.

## Figure map

| PDF page | Figure or claim | Generating code | Repository output | Status |
| --- | --- | --- | --- | --- |
| 3 | Up/down continuation tumor branches | `code/post_meeting/b6_equilibrium_jacobian/b6_equilibrium_jacobian_nextstep.jl`; cleaned by `code/post_meeting/up_down_continuation/phase1_clean_branch_plots.jl` | `results/post_meeting/phase1_clean_branches/phase1_clean_tumor_branches.png` | Matched |
| 4 | Hysteresis gap near the `b6` threshold | Same as page 3 | `results/post_meeting/phase1_clean_branches/phase1_clean_hysteresis_gap.png` | Matched |
| 5 | Stability and equilibrium residual checks | Same as page 3 | `results/post_meeting/phase1_clean_branches/phase1_clean_stability_screen.png`; `phase1_clean_residual_check.png` | Matched |
| 6 | Formal high- and low-tumor equilibrium branches | `code/post_meeting/formal_b6_continuation/b6_formal_palc_continuation_fixed.jl` | `results/post_meeting/formal_b6_continuation/equilibrium_branches_across_b6.png` | Matched |
| 7 | Stable and unstable equilibria meet | `code/post_meeting/zero_eigenvalue_classification/classify_zero_eigenvalue_special_point.jl` | `results/post_meeting/zero_eigenvalue_classification/stable_and_unstable_equilibria_meet.png` | Matched |
| 8 | Leading real eigenvalue crosses zero | Same as page 7 | `results/post_meeting/zero_eigenvalue_classification/real_eigenvalue_crosses_zero.png` | Matched |
| 9 | One isolated zero singular value | Same as page 7 | `results/post_meeting/zero_eigenvalue_classification/jacobian_has_one_zero_mode.png` | Matched |
| 11 | Tumor(0)-NK(0) basin maps across three `b6` values | `code/post_meeting/basin_of_attraction/b6_basin_of_attraction_analysis.jl`; summarized by `basin_clean_summary_figures.jl` | `results/post_meeting/basin_of_attraction/summary/basin_tumor_nk_across_b6.png` | Matched |
| 12 | Controlled basin fraction versus `b6` | `code/post_meeting/basin_of_attraction/basin_clean_summary_figures.jl` | `results/post_meeting/basin_of_attraction/summary/controlled_basin_fraction_clean.png` | Matched |
| 13 | Same `b6`, different tumor trajectories | `code/post_meeting/basin_of_attraction/basin_direction1_representative_trajectories.jl` | `results/post_meeting/basin_of_attraction/representative_trajectories/same_b6_different_tumor_trajectories.png` | Matched |
| 14 | Full-state trajectories for two initial states | Same as page 13 | `results/post_meeting/basin_of_attraction/representative_trajectories/different_initial_states_same_b6_full_trajectories.png` | Matched |
| 15 | Same initial condition, three different `b6` values | No exact standalone image or explicit generating routine was found in the supplied archive | Figure remains embedded in `June_update.pdf` | Missing original plot code |
| 16 | Three initial-condition basin maps at `b6 = 1.5e-6` | `code/post_meeting/basin_of_attraction/b6_basin_of_attraction_analysis.jl`; summarized by `basin_clean_summary_figures.jl` | `results/post_meeting/basin_of_attraction/summary/basin_three_maps_b6_1p5e-6.png` | Matched |
| 17 | Refined Tumor(0)-NK(0) phase map | `code/post_meeting/basin_of_attraction/basin_direction2_refined_tumor_nk_boundary.jl` | `results/post_meeting/basin_of_attraction/refined_boundary/refined_tumor0_nk0_basin_phase_map.png` | Matched |
| 18 | Minimum NK(0) required for control | Same as page 17 | `results/post_meeting/basin_of_attraction/refined_boundary/refined_tumor0_nk0_control_boundary.png` | Matched |
| 19, left | Formal saddle-node curve in the `b5`-`b6` plane | `code/post_meeting/two_parameter_fold/continue_saddle_node_in_b5_b6_fixed_v2.jl` | `results/post_meeting/two_parameter_fold/saddle_node_threshold_in_b5_b6_plane.png` | Matched |
| 19, right | Finite-time `b5`-`b6` heatmap reference | `code/post_meeting/step3_heatmap/step3_fix_heatmap_plotting.jl` | `results/post_meeting/step3_heatmap/fixed/step3_fixed_log_heatmap.png` | Matched |
| 20 | Three-panel `critical b6`, `B*`, and `critical b6 x B*` trade-off summary | No exact standalone image or generating routine was found. The archived two-parameter script exports the fold curve and its numerical CSV, but not this three-panel composition. | Figure remains embedded in `June_update.pdf`; source data are available in `results/post_meeting/two_parameter_fold/saddle_node_curve_b5_b6.csv` | Missing original plot code |
| 21 | Fold-curve equilibrium residuals | `code/post_meeting/two_parameter_fold/continue_saddle_node_in_b5_b6_fixed_v2.jl` | `results/post_meeting/two_parameter_fold/saddle_node_curve_residuals.png` | Matched |

Pages 1, 2, 10, and 22 contain titles, workflow summaries, classification definitions, or conclusions rather than research plots requiring separate provenance.

## Audit conclusion

Seventeen of the nineteen research figures were matched to archived source code and outputs. The two unmatched items are presentation-only composites on pages 15 and 20. Related simulation or continuation data exist, but no exact script producing those layouts was present in the supplied archive, so this repository does not claim reconstructed code as original provenance.
