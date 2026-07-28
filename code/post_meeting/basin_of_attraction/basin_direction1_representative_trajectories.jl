using DifferentialEquations
using Plots
using CSV
using DataFrames

# ------------------------------------------------------------
# Direction 1: representative full-state trajectories
# ------------------------------------------------------------
# Goal:
#   Show that at the SAME b6 value, two different initial states can lead to
#   different long-time tumor regimes.
#
# This is the trajectory-based companion to the basin map.
# Prof. MacLean likes trajectory plots, so this figure is designed to be
# presentation-ready and easy to interpret.
#
# Figure-title rule:
#   Use descriptive result-based titles only.
# ------------------------------------------------------------

outdir = "basin_direction1_trajectories"
mkpath(outdir)

B6_TARGET = 1.5e-6
B5_FIXED = 1.0e-4
TMAX = 10000.0

T_CONTROL = 1.0e3
T_DOMINANT = 1.0e5

function first_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    error("Could not find basin_all_runs.csv")
end

basin_file = first_existing([
    joinpath("basin_of_attraction_outputs", "basin_all_runs.csv"),
    "basin_all_runs.csv",
    joinpath("results", "post_meeting", "basin_of_attraction", "basin_all_runs.csv")
])

basin_rows = CSV.read(basin_file, DataFrame)

# ----------------------------
# Model definition
# ----------------------------
function base_p(; b5=B5_FIXED, b6=B6_TARGET)
    return (
        a1 = 1.0e-1,
        h  = 1.0e7,
        b1 = 3.5e-6,
        b2 = 1.1e-7,
        z1 = 0.0,

        a2 = 1.0e2,
        a3 = 1.0e8,
        z2 = 2.0e-1,

        a4 = 1.4e4,
        a5 = 2.5e-2,
        b3 = 4.0e-5,
        z3 = 4.12e-2,

        a6 = 1.1e-7,
        a7 = 1.0e-1,
        b4 = 1.0e-4,
        z4 = 2.0e-2,

        g1 = 1.0e10,
        g2 = 2.02e7,
        g3 = 2.02e7,

        a8  = 1.0e4,
        a9  = 5.0e-2,
        b5  = b5,
        z5  = 2.0e-2,
        g4  = 2.02e7,
        g5  = 1.0e3,
        a10 = 2.0,
        a11 = 1.0e-1,
        b6  = b6
    )
end

function baseline_state(p; tumor_cells=2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function model!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)

    b_to_ctl = p.a11 * xB / (p.g5 + xB)
    growth = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0))

    du[1] = growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT
    du[2] = p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM

    nk_recruit = p.a5 * xT^2/(p.g2 + xT^2) * (1.0 + p.a10*xB/(p.g5 + xB))
    du[3] = p.a4 + nk_recruit - p.b3*xM*xN - p.z3*xN

    du[4] = p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) + b_to_ctl - p.b4*xM*xC - p.z4*xC
    du[5] = p.a8 + p.a9*xT^2/(p.g4 + xT^2) - p.b5*xM*xB - p.z5*xB
end

function bad_state(u, t, integrator)
    any(x -> !isfinite(x) || x < -1e-6, u)
end
cb_bad = DiscreteCallback(bad_state, terminate!)

function solve_case(; b6=B6_TARGET, tumor0=2.0, nk0=nothing, b0=nothing, saveat=10.0)
    p = base_p(b6=b6)
    u0 = baseline_state(p; tumor_cells=tumor0)

    if nk0 !== nothing
        u0[3] = nk0
    end
    if b0 !== nothing
        u0[5] = b0
    end

    prob = ODEProblem(model!, u0, (0.0, TMAX), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x -> x < -1e-6, u))
    return sol, p, u0
end

safe(v; eps=1e-8) = max.(collect(v), eps)

# ----------------------------
# Pick representative cases from Tumor(0) × NK(0), b6 closest to target
# ----------------------------
sub = basin_rows[basin_rows.map_name .== "Tumor0_NK0", :]
b6_values = sort(unique(sub.b6))
b6_use = b6_values[argmin(abs.(b6_values .- B6_TARGET))]
sub = sub[sub.b6 .== b6_use, :]

controlled = sub[sub.final_regime .== "controlled", :]
dominant = sub[sub.final_regime .== "tumor_dominant", :]

if nrow(controlled) == 0 || nrow(dominant) == 0
    error("Need both controlled and tumor-dominant cases at b6=$(b6_use). Try a different B6_TARGET.")
end

# Choose representative points not at extreme corners if possible.
controlled = sort(controlled, [:tumor0, :nk0])
dominant = sort(dominant, [:tumor0, :nk0])

c_case = controlled[round(Int, nrow(controlled)/2), :]
d_case = dominant[round(Int, nrow(dominant)/2), :]

sol_c, _, _ = solve_case(b6=b6_use, tumor0=c_case.tumor0, nk0=c_case.nk0, saveat=10.0)
sol_d, _, _ = solve_case(b6=b6_use, tumor0=d_case.tumor0, nk0=d_case.nk0, saveat=10.0)

cell_names = ["Tumor", "MDSC", "NK", "CTL", "B"]
cell_cols = [:firebrick, :orange, :forestgreen, :royalblue, :purple]

function make_full_state_plot(sol, title_text)
    p = plot(
        sol.t, safe(sol[1,:]),
        yscale=:log10,
        lw=2.8,
        color=cell_cols[1],
        label=cell_names[1],
        xlabel="Time (days)",
        ylabel="Cells",
        title=title_text,
        size=(850, 520),
        left_margin=8Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=7Plots.mm
    )
    for j in 2:5
        plot!(p, sol.t, safe(sol[j,:]), yscale=:log10, lw=2.5,
              color=cell_cols[j], label=cell_names[j])
    end
    hline!(p, [T_CONTROL], linestyle=:dash, color=:black, label="T = 1e3")
    hline!(p, [T_DOMINANT], linestyle=:dot, color=:gray, label="T = 1e5")
    return p
end

p_control = make_full_state_plot(
    sol_c,
    "Controlled outcome from one initial state"
)

p_dominant = make_full_state_plot(
    sol_d,
    "Tumor-dominant outcome from another initial state"
)

p_combined = plot(
    p_control, p_dominant,
    layout=(1,2),
    size=(1600, 560),
    plot_title="Different initial states at the same b6 produce different tumor trajectories"
)

savefig(p_combined, joinpath(outdir, "different_initial_states_same_b6_full_trajectories.png"))

# Tumor-only version, easier for talks.
p_tumor = plot(
    sol_d.t, safe(sol_d[1,:]),
    yscale=:log10,
    lw=3,
    color=:firebrick,
    label="tumor-dominant example",
    xlabel="Time (days)",
    ylabel="Tumor cells",
    title="Same b6 can lead to different long-time tumor regimes",
    size=(900, 550),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)
plot!(p_tumor, sol_c.t, safe(sol_c[1,:]), yscale=:log10, lw=3,
      color=:forestgreen, label="controlled example")
hline!(p_tumor, [T_CONTROL], linestyle=:dash, color=:black, label="T = 1e3")
hline!(p_tumor, [T_DOMINANT], linestyle=:dot, color=:gray, label="T = 1e5")

savefig(p_tumor, joinpath(outdir, "same_b6_different_tumor_trajectories.png"))

summary = DataFrame(
    outcome = ["controlled", "tumor_dominant"],
    b6 = [b6_use, b6_use],
    tumor0 = [c_case.tumor0, d_case.tumor0],
    nk0 = [c_case.nk0, d_case.nk0],
    final_tumor = [sol_c[1,end], sol_d[1,end]]
)

CSV.write(joinpath(outdir, "representative_trajectory_cases.csv"), summary)

println("Representative trajectory figures saved in: ", outdir)
println("Main figure: different_initial_states_same_b6_full_trajectories.png")
println("Talk-friendly figure: same_b6_different_tumor_trajectories.png")
