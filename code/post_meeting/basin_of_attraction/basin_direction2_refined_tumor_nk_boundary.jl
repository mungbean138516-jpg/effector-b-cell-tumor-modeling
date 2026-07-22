using DifferentialEquations
using Plots
using DataFrames
using CSV

# ------------------------------------------------------------
# Direction 2: refined Tumor(0) × NK(0) basin boundary
# ------------------------------------------------------------
# Goal:
#   Refine the most important basin map:
#      Tumor(0) × NK(0) at b6 = 1.5e-6
#
# Why:
#   The first basin map showed that tumor/NK starting balance matters.
#   This script increases resolution near that boundary and makes cleaner,
#   presentation-ready figures.
#
# Figure-title rule:
#   Use descriptive result-based titles only.
# ------------------------------------------------------------

outdir = "basin_direction2_refined_boundary"
mkpath(outdir)

B6_TARGET = 1.5e-6
B5_FIXED = 1.0e-4
TMAX = 10000.0

# First pass grid. If this runs too slowly, reduce GRID_N to 25.
GRID_N = 41

# Refined ranges based on previous basin maps:
# Tumor(0): wide, but dense enough for boundary
# NK(0): include low, baseline, and high regions
tumor0_grid = 10 .^ range(0, 7, length=GRID_N)      # 1 to 1e7
nk0_grid = 10 .^ range(2, 6, length=GRID_N)         # 1e2 to 1e6

T_CONTROL = 1.0e3
T_DOMINANT = 1.0e5

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

function solve_case(tumor0, nk0; b6=B6_TARGET, tmax=TMAX)
    p = base_p(b6=b6)
    u0 = baseline_state(p; tumor_cells=tumor0)
    u0[3] = nk0

    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=[365.0, tmax],
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x -> x < -1e-6, u))
    return sol
end

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

rows = DataFrame()

println("Running refined Tumor(0) × NK(0) basin boundary at b6=$(B6_TARGET)")
println("Grid: ", length(tumor0_grid), " × ", length(nk0_grid), " = ", length(tumor0_grid)*length(nk0_grid), " ODE runs")

for nk0 in nk0_grid
    for tumor0 in tumor0_grid
        sol = solve_case(tumor0, nk0)
        T365 = sol[1,1]
        Tfinal = sol[1,end]

        push!(rows, (
            b6 = B6_TARGET,
            tumor0 = tumor0,
            nk0 = nk0,
            T365 = T365,
            Tfinal = Tfinal,
            logTfinal = log10(Tfinal + 1.0),
            final_regime = final_regime(Tfinal),
            regime_code = regime_code(Tfinal)
        ); cols=:union)
    end
end

CSV.write(joinpath(outdir, "refined_tumor0_nk0_basin_b6_1p5e-6.csv"), rows)

# Build matrix
Z = fill(NaN, length(nk0_grid), length(tumor0_grid))
Zlog = fill(NaN, length(nk0_grid), length(tumor0_grid))

for (i, nk0) in enumerate(nk0_grid)
    for (j, tumor0) in enumerate(tumor0_grid)
        sub = rows[(rows.nk0 .== nk0) .& (rows.tumor0 .== tumor0), :]
        if nrow(sub) == 1
            Z[i,j] = sub.regime_code[1]
            Zlog[i,j] = sub.logTfinal[1]
        end
    end
end

phase_colors = cgrad([:forestgreen, :khaki, :firebrick], categorical=true)

# Discrete basin map
p_phase = heatmap(
    log10.(tumor0_grid), log10.(nk0_grid), Z,
    color=phase_colors,
    clims=(0,2),
    xlabel="log10(Tumor(0))",
    ylabel="log10(NK(0))",
    title="Tumor and NK initial levels separate long-time tumor regimes",
    colorbar_ticks=([0,1,2], ["controlled", "intermediate", "dominant"]),
    size=(900, 700),
    right_margin=10Plots.mm,
    bottom_margin=8Plots.mm,
    left_margin=8Plots.mm,
    top_margin=8Plots.mm
)

# baseline marker
baseline = baseline_state(base_p())
scatter!(p_phase, [log10(2.0)], [log10(baseline[3])], marker=:star5, ms=10, color=:black, label="baseline")

savefig(p_phase, joinpath(outdir, "refined_tumor0_nk0_basin_phase_map.png"))

# Continuous logT map
p_log = heatmap(
    log10.(tumor0_grid), log10.(nk0_grid), Zlog,
    color=cgrad(:Blues),
    xlabel="log10(Tumor(0))",
    ylabel="log10(NK(0))",
    title="Final tumor burden across Tumor(0) and NK(0)",
    colorbar_title="log10(T(10000)+1)",
    size=(900, 700),
    right_margin=10Plots.mm,
    bottom_margin=8Plots.mm,
    left_margin=8Plots.mm,
    top_margin=8Plots.mm
)
scatter!(p_log, [log10(2.0)], [log10(baseline[3])], marker=:star5, ms=10, color=:black, label="baseline")

savefig(p_log, joinpath(outdir, "refined_tumor0_nk0_basin_logT_map.png"))

# Boundary summary by tumor0:
# for each tumor0, find the lowest NK(0) producing controlled final state
boundary = DataFrame(tumor0=Float64[], nk0_control_boundary=Float64[])

for tumor0 in tumor0_grid
    sub = rows[rows.tumor0 .== tumor0, :]
    sub = sort(sub, :nk0)
    idx = findfirst(sub.final_regime .== "controlled")
    nk_boundary = isnothing(idx) ? NaN : sub.nk0[idx]
    push!(boundary, (tumor0, nk_boundary))
end

CSV.write(joinpath(outdir, "refined_tumor0_nk0_control_boundary.csv"), boundary)

p_boundary = plot(
    boundary.tumor0, boundary.nk0_control_boundary,
    xscale=:log10,
    yscale=:log10,
    marker=:circle,
    lw=2.8,
    xlabel="Tumor(0)",
    ylabel="minimum NK(0) for controlled outcome",
    title="More initial tumor requires more NK for control",
    legend=false,
    size=(850, 550),
    left_margin=8Plots.mm,
    bottom_margin=7Plots.mm,
    top_margin=7Plots.mm
)
savefig(p_boundary, joinpath(outdir, "refined_tumor0_nk0_control_boundary.png"))

println("Refined Tumor(0) × NK(0) basin analysis complete.")
println("Outputs saved to: ", outdir)
println("Main figure: refined_tumor0_nk0_basin_phase_map.png")
println("Support figures: refined_tumor0_nk0_basin_logT_map.png, refined_tumor0_nk0_control_boundary.png")
