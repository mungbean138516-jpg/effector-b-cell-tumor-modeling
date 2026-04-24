using DifferentialEquations
using Plots
using DataFrames
using CSV
using Statistics

# Step 3: refined b5 x b6 phase diagram
# Goal:
#   turn the b5 effect into a clearer long-time phase-style result.

outdir = "step3_b5_b6_outputs"
mkpath(outdir)

t_control = 1.0e3
t_intermediate = 1.0e5

function logspace10(lo, hi, n)
    return 10 .^ range(log10(lo), log10(hi), length=n)
end

safe(v; eps=1e-8) = max.(collect(v), eps)

function nearest_index(times, day)
    return argmin(abs.(times .- day))
end

function value_at(sol, state_index, day)
    idx = nearest_index(sol.t, day)
    return sol[state_index, idx]
end

function phase_label(T)
    if T < t_control
        return "controlled"
    elseif T < t_intermediate
        return "intermediate"
    else
        return "tumor_dominant"
    end
end

function phase_code(T)
    if T < t_control
        return 0
    elseif T < t_intermediate
        return 1
    else
        return 2
    end
end

function base_p(; b5=1.0e-4, b6=1.0e-7)
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

# Baseline initial condition used in Step 1 and Step 2.
function baseline_state(p; tumor_cells=2.0)
    xM0 = p.a2 / p.z2
    xN0 = p.z2 * p.a4 / (p.a2 * p.b3 + p.z2 * p.z3)
    xC0 = 0.0
    xB0 = p.a8 / (p.z5 + p.b5 * xM0)
    return [tumor_cells, xM0, xN0, xC0, xB0]
end

function model!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)

    # direct B -> CTL term fixed for this round
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

function solve_one_case(b5, b6; tmax=10000.0, saveat=[365.0, 10000.0])
    p = base_p(b5=b5, b6=b6)
    u0 = baseline_state(p)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x->x < -1e-6, u))
    return sol, p, u0
end

function metrics(sol, b5, b6)
    T365    = value_at(sol, 1, 365.0)
    T10000  = value_at(sol, 1, 10000.0)
    M10000  = value_at(sol, 2, 10000.0)
    NK10000 = value_at(sol, 3, 10000.0)
    CTL10000= value_at(sol, 4, 10000.0)
    B10000  = value_at(sol, 5, 10000.0)

    return (
        b5=b5,
        b6=b6,
        T365=T365,
        T10000=T10000,
        logT10000=log10(T10000 + 1.0),
        M10000=M10000,
        NK10000=NK10000,
        CTL10000=CTL10000,
        B10000=B10000,
        phase=phase_label(T10000),
        phase_code=phase_code(T10000)
    )
end

# -------------------------
# Step 3 grid
# -------------------------
b5_grid = logspace10(1e-5, 1e-3, 35)
b6_grid = logspace10(2e-7, 3e-6, 35)

println("Step 3: refined b5 x b6 phase diagram")
println("  CTL formulation: direct")
println("  baseline initial condition")
println("  tmax = 10,000 days")
println("  b5 grid points: ", length(b5_grid))
println("  b6 grid points: ", length(b6_grid))
println("  total ODE runs: ", length(b5_grid) * length(b6_grid))

# save the baseline IC values for the summary slide / notes
p_base = base_p()
u0_base = baseline_state(p_base)
baseline_ic = DataFrame(
    variable = ["Tumor", "MDSC", "NK", "CTL", "B"],
    baseline_ic = u0_base
)
CSV.write(joinpath(outdir, "step3_baseline_ic.csv"), baseline_ic)


# run the 2D grid

rows = DataFrame()
for (i, b5) in enumerate(b5_grid)
    println("Running b5 = ", b5, " (", i, "/", length(b5_grid), ")")
    for b6 in b6_grid
        sol, _, _ = solve_one_case(b5, b6)
        push!(rows, metrics(sol, b5, b6); cols=:union)
    end
end

CSV.write(joinpath(outdir, "step3_b5_b6_grid.csv"), rows)
println("Saved CSV: ", joinpath(outdir, "step3_b5_b6_grid.csv"))


# reshape helpers for heatmaps

function build_matrix(df, colname, b5_values, b6_values)
    Z = fill(NaN, length(b5_values), length(b6_values))
    for (i, b5) in enumerate(b5_values)
        for (j, b6) in enumerate(b6_values)
            sub = df[(df.b5 .== b5) .& (df.b6 .== b6), :]
            if nrow(sub) == 1
                Z[i, j] = sub[1, colname]
            end
        end
    end
    return Z
end

z_log = build_matrix(rows, :logT10000, b5_grid, b6_grid)
z_phase = build_matrix(rows, :phase_code, b5_grid, b6_grid)
z_t = build_matrix(rows, :T10000, b5_grid, b6_grid)


# clean continuous heatmap

p_heat = heatmap(
    b6_grid, b5_grid, z_log,
    xscale=:log10,
    yscale=:log10,
    color=cgrad(:Blues),
    xlabel="b6",
    ylabel="b5",
    title="Step 3: b5 × b6 phase diagram (log10(T(10000)+1))",
    colorbar_title="log10(T(10000)+1)",
    size=(900, 700),
    right_margin=8Plots.mm,
    bottom_margin=6Plots.mm,
    left_margin=8Plots.mm,
    top_margin=6Plots.mm,
)

# add a contour at T(10000) = 1e3 to show the control boundary more clearly
contour!(
    p_heat,
    b6_grid, b5_grid, z_t,
    levels=[t_control],
    color=:black,
    linewidth=2.2,
    label="T(10000)=1e3"
)

savefig(p_heat, joinpath(outdir, "step3_heatmap_logT10000.png"))


# discrete phase heatmap
# 0 = controlled, 1 = intermediate, 2 = tumor-dominant
phase_colors = cgrad([:white, :lightskyblue, :navy])
p_phase = heatmap(
    b6_grid, b5_grid, z_phase,
    xscale=:log10,
    yscale=:log10,
    color=phase_colors,
    clims=(0,2),
    xlabel="b6",
    ylabel="b5",
    title="Step 3: b5 × b6 discrete phase map",
    colorbar_ticks=([0,1,2], ["controlled", "intermediate", "dominant"]),
    size=(900, 700),
    right_margin=12Plots.mm,
    bottom_margin=6Plots.mm,
    left_margin=8Plots.mm,
    top_margin=6Plots.mm,
)
savefig(p_phase, joinpath(outdir, "step3_heatmap_phase.png"))


# estimate the control boundary for each b5

boundary_rows = DataFrame(b5=Float64[], b6_control_boundary=Float64[])
for b5 in b5_grid
    sub = rows[rows.b5 .== b5, :]
    sub = sort(sub, :b6)
    idx = findfirst(sub.T10000 .< t_control)
    boundary_b6 = isnothing(idx) ? NaN : sub.b6[idx]
    push!(boundary_rows, (b5, boundary_b6))
end
CSV.write(joinpath(outdir, "step3_control_boundary.csv"), boundary_rows)

p_boundary = plot(
    boundary_rows.b5, boundary_rows.b6_control_boundary,
    marker=:circle,
    lw=2.5,
    xscale=:log10,
    yscale=:log10,
    xlabel="b5",
    ylabel="control boundary b6",
    title="Step 3: b5 shifts the b6 control boundary",
    legend=false,
    size=(800, 550)
)
savefig(p_boundary, joinpath(outdir, "step3_b6_boundary_vs_b5.png"))


# selected b5 slices: easier to present than the full heatmap alone

selected_b5 = [1e-5, 1e-4, 1e-3]
p_slices = plot(
    xlabel="b6",
    ylabel="T(10000)",
    xscale=:log10,
    yscale=:log10,
    title="Selected b5 slices through the Step 3 phase diagram",
    size=(850, 550)
)

for b5 in selected_b5
    sub = rows[isapprox.(rows.b5, b5; atol=0.0, rtol=1e-12), :]
    sub = sort(sub, :b6)
    plot!(p_slices, sub.b6, safe(sub.T10000), lw=2.5, marker=:circle, label="b5=$(b5)")
end
hline!(p_slices, [t_control], linestyle=:dash, color=:black, label="T(10000)=1e3")
savefig(p_slices, joinpath(outdir, "step3_selected_b5_slices.png"))

println("Step 3 complete.")
println("Files saved in: ", outdir)
println("Suggested main figure for slides: step3_heatmap_logT10000.png")
println("Suggested support figures: step3_heatmap_phase.png, step3_b6_boundary_vs_b5.png, step3_selected_b5_slices.png")
