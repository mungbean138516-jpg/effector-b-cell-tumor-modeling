using DifferentialEquations
using Plots
using DataFrames
using CSV
using Statistics

# Step 2 updated version
# initial-condition robustness for the long-time b6 threshold

outdir = "b6_ic_robustness_outputs_v3"
mkpath(outdir)

default(fontfamily="Arial", size=(1400,900), dpi=220, legendfontsize=10,
        guidefontsize=13, tickfontsize=10, titlefontsize=18)

# Helper functions
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

function time_to_threshold(sol, threshold)
    for i in eachindex(sol.t)
        if sol[1,i] >= threshold
            return sol.t[i]
        end
    end
    return NaN
end

function auc_until(sol, tmax)
    keep = findall(t -> t <= tmax, sol.t)
    if length(keep) <= 1
        return NaN
    end
    t = sol.t[keep]
    y = sol[1, keep]
    total = 0.0
    for i in 1:(length(t)-1)
        total += 0.5 * (y[i] + y[i+1]) * (t[i+1] - t[i])
    end
    return total
end

function state_label(T)
    if T < 1e3
        return "controlled"
    elseif T < 1e5
        return "intermediate"
    else
        return "tumor_dominant"
    end
end

# Parameters
function base_p(; b6=1e-7)
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
        b5  = 1.0e-4,
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

# save baseline IC summary
function save_baseline_ic_summary(p)
    u0 = baseline_state(p)
    df = DataFrame(
        variable = ["Tumor(0)", "MDSC(0)", "NK(0)", "CTL(0)", "B(0)"],
        baseline_value = u0,
        note = [
            "tumor seed used in Step 2",
            "tumor-free equilibrium",
            "tumor-free equilibrium",
            "tumor-free equilibrium (0 in baseline)",
            "B steady state from the added B equation"
        ]
    )
    CSV.write(joinpath(outdir, "baseline_initial_conditions.csv"), df)

    open(joinpath(outdir, "baseline_initial_conditions.txt"), "w") do io
        println(io, "Baseline initial conditions used in Step 2")
        println(io, "----------------------------------------")
        println(io, "Tumor(0) = ", u0[1])
        println(io, "MDSC(0)  = ", u0[2])
        println(io, "NK(0)    = ", u0[3])
        println(io, "CTL(0)   = ", u0[4], "  (when CTL is perturbed, a positive reference level is used)")
        println(io, "B(0)     = ", u0[5])
        println(io)
        println(io, "Current Step 2 multiplier set:")
        println(io, "0.1x, 0.2x, 0.5x, 1x, 2x, 5x, 10x")
    end
end

function start_state_with_multiplier(p, changed_comp::String, mult::Float64; tumor_cells=2.0)
    u0 = baseline_state(p; tumor_cells=tumor_cells)
    comp_index = Dict("Tumor"=>1, "MDSC"=>2, "NK"=>3, "CTL"=>4, "B"=>5)
    idx = comp_index[changed_comp]

    if changed_comp == "CTL" && u0[idx] == 0.0
        ctl_ref = 1654.4182572290226
        u0[idx] = mult * ctl_ref
    else
        u0[idx] *= mult
    end
    return u0
end

# ODE model
function model!(du, u, p, t)
    xT, xM, xN, xC, xB = max.(u, 0.0)
    B_to_CTL = p.a11 * xB / (p.g5 + xB)

    growth = p.a1 * xT * log(max(p.h / max(xT, 1e-12), 1.0))
    du[1] = growth - p.b1*xT*xN - p.b2*xT*xC - p.b6*xT*xB - p.z1*xT
    du[2] = p.a2 + p.a3*xT/(p.g1 + xT) - p.z2*xM
    nk_recruit = p.a5 * xT^2/(p.g2 + xT^2) * (1.0 + p.a10*xB/(p.g5 + xB))
    du[3] = p.a4 + nk_recruit - p.b3*xM*xN - p.z3*xN
    du[4] = p.a6*xT*xN + p.a7*xT^2/(p.g3 + xT^2) + B_to_CTL - p.b4*xM*xC - p.z4*xC
    du[5] = p.a8 + p.a9*xT^2/(p.g4 + xT^2) - p.b5*xM*xB - p.z5*xB
end

function bad_state(u, t, integrator)
    any(x -> !isfinite(x) || x < -1e-6, u)
end
cb_bad = DiscreteCallback(bad_state, terminate!)

function solve_case(b6, changed_comp, mult; tmax=10000.0, saveat=5.0)
    p = base_p(b6=b6)
    u0 = start_state_with_multiplier(p, changed_comp, mult)
    prob = ODEProblem(model!, u0, (0.0, tmax), p)
    sol = solve(prob, Tsit5();
                saveat=saveat,
                reltol=1e-7,
                abstol=1e-9,
                callback=cb_bad,
                isoutofdomain=(u,p,t)->any(x->x < -1e-6, u))
    return sol
end

function metrics(sol, b6, changed_comp, mult)
    T100    = value_at(sol, 1, 100.0)
    T365    = value_at(sol, 1, 365.0)
    T10000  = value_at(sol, 1, 10000.0)
    t1e5 = time_to_threshold(sol, 1e5)
    t1e6 = time_to_threshold(sol, 1e6)
    return (
        changed_comp = changed_comp,
        multiplier = mult,
        b6 = b6,
        T100 = T100,
        T365 = T365,
        T10000 = T10000,
        logT10000 = log10(T10000 + 1.0),
        M10000 = value_at(sol, 2, 10000.0),
        NK10000 = value_at(sol, 3, 10000.0),
        CTL10000 = value_at(sol, 4, 10000.0),
        B10000 = value_at(sol, 5, 10000.0),
        t_to_1e5 = t1e5,
        t_to_1e6 = t1e6,
        t_to_1e5_plot = isnan(t1e5) ? 10000.0 : t1e5,
        t_to_1e6_plot = isnan(t1e6) ? 10000.0 : t1e6,
        auc_0_365 = auc_until(sol, 365.0),
        auc_0_1000 = auc_until(sol, 1000.0),
        final_state = state_label(T10000)
    )
end

# Step 2 design
changed_comps = ["Tumor", "MDSC", "NK", "CTL", "B"]
mults = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
b6_grid = logspace10(3e-7, 2e-6, 40)
traj_b6 = 8.0e-7

println("Step 2 (updated): initial-condition robustness")
println("  CTL formulation: direct")
println("  b5 fixed at 1e-4")
println("  changed compartments: ", changed_comps)
println("  multipliers: ", mults)
println("  b6 grid points: ", length(b6_grid))
println("  total ODE runs: ", length(changed_comps)*length(mults)*length(b6_grid))

save_baseline_ic_summary(base_p())

# Run all cases
all_rows = DataFrame()
for comp in changed_comps
    for mult in mults
        for (k, b6) in enumerate(b6_grid)
            println("Running comp=", comp, " mult=", mult, " b6=", b6, " (", k, "/", length(b6_grid), ")")
            sol = solve_case(b6, comp, mult)
            push!(all_rows, metrics(sol, b6, comp, mult); cols=:union)
        end
    end
end
CSV.write(joinpath(outdir, "b6_ic_robustness_all_runs.csv"), all_rows)

# Threshold summary
threshold_rows = DataFrame()
for comp in changed_comps
    for mult in mults
        sub = all_rows[(all_rows.changed_comp .== comp) .& (all_rows.multiplier .== mult), :]
        sort!(sub, :b6)
        th_1e5 = NaN; th_1e4 = NaN; th_1e3 = NaN
        for r in eachrow(sub)
            if isnan(th_1e5) && r.T10000 < 1e5
                th_1e5 = r.b6
            end
            if isnan(th_1e4) && r.T10000 < 1e4
                th_1e4 = r.b6
            end
            if isnan(th_1e3) && r.T10000 < 1e3
                th_1e3 = r.b6
            end
        end
        push!(threshold_rows, (changed_comp=comp, multiplier=mult,
                              threshold_T1e5=th_1e5,
                              threshold_T1e4=th_1e4,
                              threshold_T1e3=th_1e3); cols=:union)
    end
end
CSV.write(joinpath(outdir, "b6_ic_robustness_thresholds.csv"), threshold_rows)

# Plot 1 relative threshold view
function plot_thresholds_relative()
    comps_colors = Dict("Tumor"=>:red, "MDSC"=>:orange, "NK"=>:green, "CTL"=>:blue, "B"=>:purple)
    p1 = plot(xlabel="Initial-condition multiplier", ylabel="relative threshold",
              yscale=:identity, title="Relative threshold for T(10000) < 1e5",
              legend=:outerright)
    p2 = plot(xlabel="Initial-condition multiplier", ylabel="relative threshold",
              yscale=:identity, title="Relative threshold for T(10000) < 1e3",
              legend=false)

    for comp in changed_comps
        sub = threshold_rows[threshold_rows.changed_comp .== comp, :]
        sort!(sub, :multiplier)
        base1 = sub[sub.multiplier .== 1.0, :threshold_T1e5][1]
        base2 = sub[sub.multiplier .== 1.0, :threshold_T1e3][1]
        y1 = sub.threshold_T1e5 ./ base1
        y2 = sub.threshold_T1e3 ./ base2
        plot!(p1, sub.multiplier, y1, marker=:circle, lw=2.6, color=comps_colors[comp], label=comp)
        plot!(p2, sub.multiplier, y2, marker=:circle, lw=2.6, color=comps_colors[comp], label=comp)
    end
    hline!(p1, [1.0], ls=:dash, color=:gray50, lw=1.5, label=false)
    hline!(p2, [1.0], ls=:dash, color=:gray50, lw=1.5, label=false)
    fig = plot(p1, p2, layout=(1,2), size=(1250,450), margin=5Plots.mm)
    savefig(fig, joinpath(outdir, "ic_thresholds_relative_better.png"))
end
plot_thresholds_relative()

# Plot 2: raw threshold view
function plot_thresholds_raw()
    comps_colors = Dict("Tumor"=>:red, "MDSC"=>:orange, "NK"=>:green, "CTL"=>:blue, "B"=>:purple)
    p1 = plot(xlabel="Initial-condition multiplier", ylabel="b6 threshold",
              yscale=:log10, title="Raw threshold for T(10000) < 1e5",
              legend=:outerright)
    p2 = plot(xlabel="Initial-condition multiplier", ylabel="b6 threshold",
              yscale=:log10, title="Raw threshold for T(10000) < 1e3",
              legend=false)
    for comp in changed_comps
        sub = threshold_rows[threshold_rows.changed_comp .== comp, :]
        sort!(sub, :multiplier)
        plot!(p1, sub.multiplier, sub.threshold_T1e5, marker=:circle, lw=2.6, color=comps_colors[comp], label=comp)
        plot!(p2, sub.multiplier, sub.threshold_T1e3, marker=:circle, lw=2.6, color=comps_colors[comp], label=comp)
    end
    fig = plot(p1, p2, layout=(1,2), size=(1250,450), margin=5Plots.mm)
    savefig(fig, joinpath(outdir, "ic_thresholds_raw_better.png"))
end
plot_thresholds_raw()

# Plot 3 discrete phase heatmap
function phase_code(label)
    label == "controlled" && return 0
    label == "intermediate" && return 1
    return 2
end

function phase_matrix_for(comp)
    Z = fill(NaN, length(mults), length(b6_grid))
    for (i, mult) in enumerate(mults)
        for (j, b6) in enumerate(b6_grid)
            row = all_rows[(all_rows.changed_comp .== comp) .&
                           (all_rows.multiplier .== mult) .&
                           (abs.(all_rows.b6 .- b6) .< 1e-20), :]
            if nrow(row) > 0
                Z[i,j] = phase_code(row.final_state[1])
            end
        end
    end
    return Z
end

function plot_discrete_phase_heatmaps()
    hm_plots = Any[]
    phase_colors = cgrad([:lightskyblue1, :deepskyblue3, :navy], 3, categorical=true)
    for comp in changed_comps
        Z = phase_matrix_for(comp)
        p = heatmap(log10.(b6_grid), mults, Z,
                    color=phase_colors, clims=(0,2),
                    xlabel="log10(b6)", ylabel="IC multiplier",
                    title="Change $(comp)(0)",
                    colorbar_title="regime")
        push!(hm_plots, p)
    end
    fig = plot(hm_plots..., layout=(3,2), size=(1250,1100), margin=4Plots.mm)
    savefig(fig, joinpath(outdir, "ic_discrete_phase_heatmap.png"))
end
plot_discrete_phase_heatmaps()

# Plot 4 supplementary full-state plots (NK and B focus)
function plot_full_state_focus(changed_comp::String)
    chosen_mults = [0.1, 1.0, 10.0]
    cols = Dict(0.1=>:black, 1.0=>:dodgerblue3, 10.0=>:red)
    labels = ["Tumor", "MDSC", "NK", "CTL", "B"]
    subplots = Any[]
    for state_idx in 1:5
        p = plot(xlabel=state_idx==5 ? "Time (days)" : "",
                 ylabel=labels[state_idx], yscale=:log10,
                 legend=state_idx==1 ? :outerright : false)
        for mult in chosen_mults
            sol = solve_case(traj_b6, changed_comp, mult)
            plot!(p, sol.t, safe(sol[state_idx,:]), lw=2.3, color=cols[mult], label="$(mult)x")
        end
        push!(subplots, p)
    end
    fig = plot(subplots..., layout=(5,1), size=(900,1000), margin=4Plots.mm,
               plot_title="Supplementary state trajectories: vary $(changed_comp)(0), b6=$(traj_b6)")
    savefig(fig, joinpath(outdir, "supp_state_$(lowercase(changed_comp))0_b6_8e7.png"))
end
plot_full_state_focus("NK")
plot_full_state_focus("B")

println("Done. Outputs saved to: ", outdir)
