using CSV
using DataFrames
using Plots
using Statistics

# ------------------------------------------------------------
# Phase 1 / Step 1: clean bifurcation-style plots for b6
# ------------------------------------------------------------
# Purpose:
#   This script does NOT rerun the ODE model.
#   It reads the existing up-sweep / down-sweep continuation outputs
#   and makes cleaner presentation-ready figures.
#
# Why this script exists:
#   Prof. MacLean suggested clarifying whether the b6 threshold is
#   switch-like or bifurcation-like. The up/down continuation already
#   gives a first hysteresis screen. This script makes the visualization
#   clearer by:
#     - using a plotting floor for very small tumor values
#     - separating branch, hysteresis-gap, and stability views
#     - adding threshold guide lines
#     - saving a compact numerical summary
# ------------------------------------------------------------

outdir = "phase1_clean_branch_outputs"
mkpath(outdir)

# Plotting constants
T_FLOOR = 1e-6          # plotting only; does not change data
T_CONTROL = 1e3
T_DOMINANT = 1e5
GAP_TOL = 0.5           # log10-scale gap; >0.5 means clearly separated branches

# ------------------------------------------------------------
# Find input files flexibly
# ------------------------------------------------------------
function first_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    error("Could not find any of the expected files:\n" * join(paths, "\n"))
end

branches_file = first_existing([
    "b6_equilibrium_branches.csv",
    joinpath("b6_equilibrium_jacobian_outputs", "b6_equilibrium_branches.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian", "b6_equilibrium_branches.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian", "b6_equilibrium_jacobian_outputs", "b6_equilibrium_branches.csv")
])

comparison_file = first_existing([
    "b6_equilibrium_up_down_comparison.csv",
    joinpath("b6_equilibrium_jacobian_outputs", "b6_equilibrium_up_down_comparison.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian", "b6_equilibrium_up_down_comparison.csv"),
    joinpath("results", "post_meeting", "b6_equilibrium_jacobian", "b6_equilibrium_jacobian_outputs", "b6_equilibrium_up_down_comparison.csv")
])

println("Using branches file: ", branches_file)
println("Using comparison file: ", comparison_file)

branches = CSV.read(branches_file, DataFrame)
comp = CSV.read(comparison_file, DataFrame)

# Expected columns in branches:
# direction, b6, T, max_real_eig, residual_scaled, phase
#
# Expected columns in comparison:
# b6, T_up, T_down, logT_up, logT_down, hysteresis_gap

up = sort(branches[branches.direction .== "up", :], :b6)
down = sort(branches[branches.direction .== "down", :], :b6)
comp = sort(comp, :b6)

# Helper for plotting floor
plot_T(x) = max.(x, T_FLOOR)

# ------------------------------------------------------------
# Estimate simple transition / hysteresis summaries
# ------------------------------------------------------------
function first_b6_below(df, threshold)
    idx = findfirst(df.T .< threshold)
    return isnothing(idx) ? NaN : df.b6[idx]
end

# Where branches first become clearly separated
idx_gap = findfirst(comp.hysteresis_gap .> GAP_TOL)
b6_gap_start = isnothing(idx_gap) ? NaN : comp.b6[idx_gap]

summary = DataFrame(
    quantity = [
        "gap_tolerance_log10",
        "first_b6_with_hysteresis_gap_gt_tol",
        "up_branch_first_T_below_1e3",
        "down_branch_first_T_below_1e3",
        "up_branch_first_T_below_1e5",
        "down_branch_first_T_below_1e5"
    ],
    value = [
        GAP_TOL,
        b6_gap_start,
        first_b6_below(up, T_CONTROL),
        first_b6_below(down, T_CONTROL),
        first_b6_below(up, T_DOMINANT),
        first_b6_below(down, T_DOMINANT)
    ]
)
CSV.write(joinpath(outdir, "phase1_branch_summary.csv"), summary)

# ------------------------------------------------------------
# Figure 1: clean branch plot
# ------------------------------------------------------------
p_branch = plot(
    up.b6, plot_T(up.T),
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    lw=2.8,
    color=:steelblue,
    label="up-sweep",
    xlabel="b6",
    ylabel="Tumor equilibrium candidate",
    title="Phase 1: up/down continuation reveals hysteresis-like branches",
    size=(950, 650),
    legend=:bottomleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=6Plots.mm,
    right_margin=8Plots.mm,
)
plot!(p_branch, down.b6, plot_T(down.T), marker=:diamond, lw=2.8, color=:darkorange, label="down-sweep")
hline!(p_branch, [T_CONTROL], linestyle=:dash, color=:black, label="T = 1e3")
hline!(p_branch, [T_DOMINANT], linestyle=:dot, color=:gray, label="T = 1e5")
if !isnan(b6_gap_start)
    vline!(p_branch, [b6_gap_start], linestyle=:dashdot, color=:red, label="gap starts")
end
savefig(p_branch, joinpath(outdir, "phase1_clean_tumor_branches.png"))

# ------------------------------------------------------------
# Figure 2: hysteresis gap plot
# ------------------------------------------------------------
p_gap = plot(
    comp.b6, comp.hysteresis_gap,
    xscale=:log10,
    marker=:circle,
    lw=2.8,
    color=:purple,
    xlabel="b6",
    ylabel="|log10(T_up + 1) - log10(T_down + 1)|",
    title="Phase 1: hysteresis gap between up/down branches",
    legend=false,
    size=(950, 550),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=6Plots.mm,
)
hline!(p_gap, [GAP_TOL], linestyle=:dash, color=:black, label=false)
if !isnan(b6_gap_start)
    vline!(p_gap, [b6_gap_start], linestyle=:dashdot, color=:red, label=false)
end
savefig(p_gap, joinpath(outdir, "phase1_clean_hysteresis_gap.png"))

# ------------------------------------------------------------
# Figure 3: stability screen
# ------------------------------------------------------------
if (:max_real_eig in propertynames(up)) && (:max_real_eig in propertynames(down))
    p_stab = plot(
        up.b6, up.max_real_eig,
        xscale=:log10,
        marker=:circle,
        lw=2.6,
        color=:steelblue,
        label="up-sweep",
        xlabel="b6",
        ylabel="max real eigenvalue",
        title="Phase 1: local stability screen at branch candidates",
        size=(950, 550),
        legend=:bottomright,
        left_margin=8Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=6Plots.mm,
    )
    plot!(p_stab, down.b6, down.max_real_eig, marker=:diamond, lw=2.6, color=:darkorange, label="down-sweep")
    hline!(p_stab, [0.0], linestyle=:dash, color=:black, label="stability boundary")
    savefig(p_stab, joinpath(outdir, "phase1_clean_stability_screen.png"))
end

# ------------------------------------------------------------
# Figure 4: residual screen if residuals are available
# ------------------------------------------------------------
if (:residual_scaled in propertynames(up)) && (:residual_scaled in propertynames(down))
    safe_res(x) = max.(x, 1e-16)
    p_res = plot(
        up.b6, safe_res(up.residual_scaled),
        xscale=:log10,
        yscale=:log10,
        marker=:circle,
        lw=2.6,
        color=:steelblue,
        label="up-sweep",
        xlabel="b6",
        ylabel="scaled residual ||f(x*)||",
        title="Phase 1: equilibrium-candidate residual check",
        size=(950, 550),
        legend=:bottomleft,
        left_margin=8Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=6Plots.mm,
    )
    plot!(p_res, down.b6, safe_res(down.residual_scaled), marker=:diamond, lw=2.6, color=:darkorange, label="down-sweep")
    hline!(p_res, [1e-6, 1e-8], linestyle=:dash, color=[:gray :black], label=["1e-6" "1e-8"])
    savefig(p_res, joinpath(outdir, "phase1_clean_residual_check.png"))
end

println("Phase 1 clean plotting complete.")
println("Outputs saved to: ", outdir)
println("Main figure: phase1_clean_tumor_branches.png")
println("Support figures: phase1_clean_hysteresis_gap.png, phase1_clean_stability_screen.png, phase1_clean_residual_check.png")
println("Summary CSV: phase1_branch_summary.csv")
