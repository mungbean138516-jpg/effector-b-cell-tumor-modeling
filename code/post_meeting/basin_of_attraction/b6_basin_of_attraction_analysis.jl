using DifferentialEquations
using Plots
using DataFrames
using CSV
using Statistics

# ------------------------------------------------------------
# Basin-of-attraction analysis near the b6 hysteresis region
# ------------------------------------------------------------
# This script answers a different question from the earlier IC robustness test.
#
# Earlier Step 2 question:
#   Does the long-time b6 threshold move when one initial condition is changed?
#
# This basin test question:
#   For the SAME fixed b6, can different starting tumor/immune states end in
#   different long-time regimes?
#
# If yes, this supports possible bistability / multiple basins of attraction.
#
# Figure-title rule:
#   Titles are descriptive only. No "Phase 1", "Step 4", etc.
# ------------------------------------------------------------

outdir = "basin_of_attraction_outputs"
mkpath(outdir)

# ----------------------------
# User-tunable settings
# ----------------------------
GRID_N = 19
TMAX = 10000.0
SAVEAT_SUMMARY = [365.0, TMAX]

# b6 values inside / near the hysteresis region.
B6_VALUES = [8.0e-7, 1.0e-6, 1.5e-6]

# fixed b5 baseline from previous analyses
B5_FIXED = 1.0e-4

# for the NK(0) × B(0) immune-context map, use a moderate tumor burden
# so that the map is not trivially controlled just because Tumor(0)=2.
FIXED_TUMOR_FOR_IMMUNE_MAP = 1.0e4

# thresholds for long-time regimes
T_CONTROL = 1.0e3
T_DOMINANT = 1.0e5

# initial-condition grids
tumor0_grid = 10 .^ range(0, 7, length=GRID_N)    # 1 to 1e7
nk0_grid    = 10 .^ range(2, 6, length=GRID_N)    # 1e2 to 1e6
b0_grid     = 10 .^ range(0, 6, length=GRID_N)    # 1 to 1e6, includes rare-B values

# ----------------------------
# Model parameters
# ----------------------------
function base_p(; b5=B5_FIXED, b6=1.0e-6)
    return (
        # original CIR-style parameters
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

        # effector-B extension
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

function solve_case(; b6=1.0e-6, tumor0=2.0, nk0=nothing, b0=nothing,
                    m0=nothing, ctl0=nothing, tmax=TMAX, saveat=SAVEAT_SUMMARY)
    p = base_p(b6=b6)
    u0 = baseline_state(p; tumor_cells=tumor0)

    if nk0 !== nothing
        u0[3] = nk0
    end
    if b0 !== nothing
        u0[5] = b0
    end
    if m0 !== nothing
        u0[2] = m0
    end
    if ctl0 !== nothing
        u0[4] = ctl0
    end

    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x -> x < -1e-6, u))
    return sol, p, u0
end

# ----------------------------
# Metrics
# ----------------------------
function final_regime(T)
    if T < T_CONTROL
        return "controlled"
    elseif T < T_DOMINANT
        return "intermediate"
    else
        return "tumor_dominant"
    end
end

function regime_code(T)
    if T < T_CONTROL
        return 0
    elseif T < T_DOMINANT
        return 1
    else
        return 2
    end
end

function value_at(sol, state_index, day)
    idx = argmin(abs.(sol.t .- day))
    return sol[state_index, idx]
end

function make_metrics(sol; map_name, b6, tumor0, nk0, b0)
    T365 = value_at(sol, 1, 365.0)
    Tfinal = value_at(sol, 1, TMAX)

    return (
        map_name = map_name,
        b6 = b6,
        tumor0 = tumor0,
        nk0 = nk0,
        b0 = b0,
        T365 = T365,
        Tfinal = Tfinal,
        logTfinal = log10(Tfinal + 1.0),
        final_regime = final_regime(Tfinal),
        regime_code = regime_code(Tfinal)
    )
end

# ----------------------------
# Run a basin map
# ----------------------------
function run_basin_map(map_name, x_grid, y_grid; b6=1.0e-6)
    rows = DataFrame()

    println("Running basin map: ", map_name, " at b6=", b6)

    for yval in y_grid
        for xval in x_grid
            if map_name == "Tumor0_NK0"
                tumor0 = xval
                nk0 = yval
                b0 = nothing
                sol, _, _ = solve_case(b6=b6, tumor0=tumor0, nk0=nk0)

            elseif map_name == "Tumor0_B0"
                tumor0 = xval
                nk0 = nothing
                b0 = yval
                sol, _, _ = solve_case(b6=b6, tumor0=tumor0, b0=b0)

            elseif map_name == "NK0_B0"
                tumor0 = FIXED_TUMOR_FOR_IMMUNE_MAP
                nk0 = yval
                b0 = xval
                sol, _, _ = solve_case(b6=b6, tumor0=tumor0, nk0=nk0, b0=b0)

            else
                error("Unknown map_name: $(map_name)")
            end

            push!(rows, make_metrics(sol;
                map_name=map_name,
                b6=b6,
                tumor0=(map_name == "NK0_B0" ? FIXED_TUMOR_FOR_IMMUNE_MAP : xval),
                nk0=(map_name == "Tumor0_NK0" ? yval : (map_name == "NK0_B0" ? yval : NaN)),
                b0=(map_name == "Tumor0_B0" ? yval : (map_name == "NK0_B0" ? xval : NaN))
            ); cols=:union)
        end
    end

    return rows
end

function build_matrix(rows, x_grid, y_grid, value_col, map_name)
    Z = fill(NaN, length(y_grid), length(x_grid))
    for (i, yval) in enumerate(y_grid)
        for (j, xval) in enumerate(x_grid)
            if map_name == "Tumor0_NK0"
                sub = rows[(rows.tumor0 .== xval) .& (rows.nk0 .== yval), :]
            elseif map_name == "Tumor0_B0"
                sub = rows[(rows.tumor0 .== xval) .& (rows.b0 .== yval), :]
            else
                sub = rows[(rows.b0 .== xval) .& (rows.nk0 .== yval), :]
            end
            if nrow(sub) == 1
                Z[i,j] = sub[1, value_col]
            end
        end
    end
    return Z
end

function plot_basin_map(rows, x_grid, y_grid; map_name, b6)
    if map_name == "Tumor0_NK0"
        xlab = "log10(Tumor(0))"
        ylab = "log10(NK(0))"
        title = "Initial tumor and NK levels determine final tumor regime"
        baseline_x = log10(2.0)
        baseline_y = log10(baseline_state(base_p())[3])
    elseif map_name == "Tumor0_B0"
        xlab = "log10(Tumor(0))"
        ylab = "log10(B(0))"
        title = "Initial tumor and B-cell levels determine final tumor regime"
        baseline_x = log10(2.0)
        baseline_y = log10(baseline_state(base_p())[5])
    else
        xlab = "log10(B(0))"
        ylab = "log10(NK(0))"
        title = "Initial NK and B-cell levels determine final tumor regime"
        baseline_x = log10(baseline_state(base_p())[5])
        baseline_y = log10(baseline_state(base_p())[3])
    end

    Z = build_matrix(rows, x_grid, y_grid, :regime_code, map_name)
    phase_colors = cgrad([:forestgreen, :khaki, :firebrick], categorical=true)

    p = heatmap(
        log10.(x_grid), log10.(y_grid), Z,
        color=phase_colors,
        clims=(0,2),
        xlabel=xlab,
        ylabel=ylab,
        title=title,
        colorbar_ticks=([0,1,2], ["controlled", "intermediate", "dominant"]),
        size=(900, 700),
        right_margin=10Plots.mm,
        bottom_margin=8Plots.mm,
        left_margin=8Plots.mm,
        top_margin=8Plots.mm
    )

    # Baseline marker
    scatter!(p, [baseline_x], [baseline_y], marker=:star5, ms=9,
             color=:black, label="baseline")

    # Put the b6 value in the plot without making title too long.
    annotate!(p, minimum(log10.(x_grid)) + 0.05*(maximum(log10.(x_grid))-minimum(log10.(x_grid))),
                 maximum(log10.(y_grid)) - 0.08*(maximum(log10.(y_grid))-minimum(log10.(y_grid))),
                 text("b6 = $(b6)", 10, :left))

    return p
end

# ----------------------------
# Main run
# ----------------------------
all_rows = DataFrame()

map_specs = [
    ("Tumor0_NK0", tumor0_grid, nk0_grid),
    ("Tumor0_B0", tumor0_grid, b0_grid),
    ("NK0_B0", b0_grid, nk0_grid)
]

for b6 in B6_VALUES
    for (map_name, x_grid, y_grid) in map_specs
        rows = run_basin_map(map_name, x_grid, y_grid; b6=b6)
        append!(all_rows, rows; cols=:union)

        CSV.write(joinpath(outdir, "basin_$(map_name)_b6_$(round(b6, sigdigits=3)).csv"), rows)

        p = plot_basin_map(rows, x_grid, y_grid; map_name=map_name, b6=b6)
        savefig(p, joinpath(outdir, "basin_$(map_name)_b6_$(round(b6, sigdigits=3)).png"))
    end
end

CSV.write(joinpath(outdir, "basin_all_runs.csv"), all_rows)

# ----------------------------
# Basin fraction summary
# ----------------------------
frac_rows = DataFrame(
    map_name=String[],
    b6=Float64[],
    frac_controlled=Float64[],
    frac_intermediate=Float64[],
    frac_dominant=Float64[]
)

for b6 in B6_VALUES
    for (map_name, _, _) in map_specs
        sub = all_rows[(all_rows.b6 .== b6) .& (all_rows.map_name .== map_name), :]
        n = nrow(sub)
        push!(frac_rows, (
            map_name,
            b6,
            sum(sub.final_regime .== "controlled") / n,
            sum(sub.final_regime .== "intermediate") / n,
            sum(sub.final_regime .== "tumor_dominant") / n
        ))
    end
end

CSV.write(joinpath(outdir, "basin_fraction_summary.csv"), frac_rows)

p_frac = plot(
    xlabel="b6",
    ylabel="fraction of IC grid ending controlled",
    xscale=:log10,
    title="Controlled basin fraction increases with b6",
    size=(850, 550),
    legend=:topleft,
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)

for map_name in unique(frac_rows.map_name)
    sub = frac_rows[frac_rows.map_name .== map_name, :]
    plot!(p_frac, sub.b6, sub.frac_controlled,
          marker=:circle, lw=2.5, label=replace(map_name, "_" => " × "))
end

savefig(p_frac, joinpath(outdir, "controlled_basin_fraction_vs_b6.png"))

# ----------------------------
# Representative trajectories from Tumor(0) × NK(0) map at b6 = 1e-6
# ----------------------------
function plot_representative_trajectories()
    target_b6 = 1.0e-6
    sub = all_rows[(all_rows.b6 .== target_b6) .& (all_rows.map_name .== "Tumor0_NK0"), :]

    controlled = sub[sub.final_regime .== "controlled", :]
    dominant = sub[sub.final_regime .== "tumor_dominant", :]

    if nrow(controlled) == 0 || nrow(dominant) == 0
        println("Could not make representative trajectory plot: missing controlled or dominant cases.")
        return
    end

    # Choose simple representative examples.
    c = controlled[round(Int, nrow(controlled)/2), :]
    d = dominant[round(Int, nrow(dominant)/2), :]

    sol_c, _, _ = solve_case(b6=target_b6, tumor0=c.tumor0, nk0=c.nk0, saveat=10.0)
    sol_d, _, _ = solve_case(b6=target_b6, tumor0=d.tumor0, nk0=d.nk0, saveat=10.0)

    p = plot(
        sol_d.t, max.(sol_d[1,:], 1e-8),
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
        top_margin=7Plots.mm,
    )
    plot!(p, sol_c.t, max.(sol_c[1,:], 1e-8),
          yscale=:log10, lw=3, color=:forestgreen,
          label="controlled example")
    hline!(p, [T_CONTROL], linestyle=:dash, color=:black, label="T = 1e3")
    hline!(p, [T_DOMINANT], linestyle=:dot, color=:gray, label="T = 1e5")

    savefig(p, joinpath(outdir, "representative_basin_trajectories_b6_1e-6.png"))
end

plot_representative_trajectories()

println("Basin-of-attraction analysis complete.")
println("Outputs saved to: ", outdir)
println("Key figures:")
println("  basin_Tumor0_NK0_b6_*.png")
println("  basin_Tumor0_B0_b6_*.png")
println("  basin_NK0_B0_b6_*.png")
println("  controlled_basin_fraction_vs_b6.png")
println("  representative_basin_trajectories_b6_1e-6.png")
